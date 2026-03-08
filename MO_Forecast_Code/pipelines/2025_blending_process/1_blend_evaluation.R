# ==============================================================================
# File: 1_blend_evaluation.R
# ==============================================================================
# Purpose
#   Cross-validate weekly-bin multinomial onset models (via nnet) using a YAML
#   spec to control included models and blocks. Saves metrics, reliability
#   plots, MME (multi-model ensemble) weights, and map_inputs RDS for downstream figure generation.
#
# Inputs
#   - specs/2025_blend/<spec_id>.yml (default: cv_models.yml)
#   - Monsoon_Data/Processed_Data/2025_pipeline_input/cv_data_*.rds (from script 0;
#     exact file depends on spec$run$cutoff_mode: mok, clim_mok_date, or no_mok_filter)
#   - Monsoon_Data/dissemination_cells.csv
#
# Outputs
#   - Monsoon_Data/results/2025_model_evaluation/* (metrics, reliability, summaries, MME weights)
#   - Monsoon_Data/results/2025_model_evaluation/map_inputs{output_tag}.rds
#
# Usage
#   Rscript pipelines/2025_blending_process/1_blend_evaluation.R [--spec_id cv_models] [--cores N]
#
# Dependencies
#   - R/2025_blending_process/blend_evaluation_utils.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tictoc)
  library(pROC)
  library(ModelMetrics)
  library(ggplot2)
  library(purrr)
  library(furrr)
  library(maps)
  library(readr)
  library(lubridate)
  library(colorspace)
  library(yaml)
  library(parallel)
  library(nnet)
  library(tibble)
})

path_box <- "Monsoon_Data/results/2025_model_evaluation"

# Project functions and utilities
source("R/2025_blending_process/blend_evaluation_utils.R")

# ---- Parse CLI args ----
if (!interactive()) {
  option_list <- list(
    optparse::make_option("--spec_id", type = "character", default = "cv_models",
                          help = "Spec file name (without .yml) in specs/2025_blend/"),
    optparse::make_option("--cores", type = "integer", default = NULL,
                          help = "Override number of parallel cores (for running variants concurrently)")
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  spec_id <- opt$spec_id
  if (!is.null(opt$cores)) options(mc.cores = opt$cores)
} else {
  spec_id <- "hindcast_1965_1978_no_mok_filter"
}

# ---- Read spec ----
spec_path <- file.path("specs", "2025_blend", paste0(spec_id, ".yml"))
spec <- yaml::read_yaml(spec_path)

cutoff_mode <- spec$run$cutoff_mode %||% "mok"
MR <- isTRUE(spec$run$MR %||% TRUE)
training_years <- as.integer(unlist(spec$run$training_years %||% (2000:2024)))
true_holdout_years <- as.integer(unlist(spec$run$true_holdout_years %||% integer(0)))
cv_holdout_years <- as.integer(unlist(spec$run$cv_holdout_years %||% integer(0)))
holdout_years <- c(true_holdout_years, cv_holdout_years)

cv_methods <- unlist(spec$cv$methods %||% c("global"))

cutoff_tag <- make_cutoff_tag(cutoff_mode)
year_tag <- make_year_tag(holdout_years)
output_tag <- paste0(cutoff_tag, year_tag)

# ---- dissemination cells ----
dissemination_cells <- readr::read_csv("Monsoon_Data/dissemination_cells.csv", show_col_types = FALSE) %>%
  dplyr::transmute(lat = as.numeric(lat), lon = as.numeric(lon)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(id = paste0(lat, "_", lon))

# ---- India boundary ----
india_boundary <- read.csv("Monsoon_Data/india_boundary_for_ggplot.csv") %>%
  dplyr::filter(id == 253)

# ---- Load data ----
input_rds_file <- input_rds_from_cutoff(cutoff_mode)

tictoc::tic()
wide_df <- readr::read_rds(file.path("Monsoon_Data/Processed_Data/2025_pipeline_input", input_rds_file)) %>%
  dplyr::filter(!is.na(outcome))

# ---- Build formula models from spec ----
formulas_list <- build_formulas_from_spec(spec, cutoff_mode)
model_names <- names(formulas_list)

# ---- Section 1: Main CV loop ----
# Outer loop over CV methods (e.g. "global"), inner loop over nnet formula models
# from the spec.  Each (method, model) pair produces leave-year-out CV predictions
# via compute_cv_*().  save_reliability_and_metrics() computes per-cell and pooled
# metrics (Brier, RPS, AUC) and writes reliability-diagram CSVs; it returns updated
# accumulators cv_metrics_list and yearly_metrics.
cv_metrics_list <- list()
yearly_metrics <- list()
clim_cv_preds_list <- list()

for (method in cv_methods) {
  cv_metrics_list[[method]] <- list()
  yearly_metrics[[method]] <- list()
  
  for (mod in model_names) {
    cat("Processing model", mod, "with", method, "CV...\n")
    
    # Dispatch to the appropriate CV strategy
    cv_preds <- switch(
      method,
      global = compute_cv_global(
        formulas_list[[mod]],
        data_train = restrict_to_allowed(wide_df, dissemination_cells),
        data_pred  = wide_df,
        holdout_years,
        true_holdout_years = true_holdout_years
      ),
      local = compute_cv_local(formulas_list[[mod]], wide_df, holdout_years),
      neighbors = compute_cv_neighbors(formulas_list[[mod]], wide_df, holdout_years),
      cluster2 = compute_cv_clusters(formulas_list[[mod]], wide_df, holdout_years, "cluster2"),
      cluster3 = compute_cv_clusters(formulas_list[[mod]], wide_df, holdout_years, "cluster3")
    )
    
    if (is.null(cv_preds)) next
    
    # Compute and save metrics; returns updated accumulator lists
    res <- save_reliability_and_metrics(
      cv_preds    = cv_preds,
      model_name  = mod,
      method      = method,
      path_box    = path_box,
      cv_metrics_list = cv_metrics_list,
      yearly_metrics  = yearly_metrics,
      allowed_cells = dissemination_cells,
      run_tag     = output_tag
    )
    cv_metrics_list <- res$cv_metrics_list
    yearly_metrics  <- res$yearly_metrics
  }
  # ---- Section 2: Climatology logit tasks ----
  # Processes climatology-based models defined in spec$extras$clim_logits.
  # Two sub-types:
  #   - "unc_clim_raw" (static climatology): converts logit columns to 5-bin + earlier
  #   - "clim_raw" (evolving expectations): 5-bin only, with optional training-window
  #     variants (e.g. clim_raw_1960_1999) via window_start_years / window_end_year.
  # Each variant is scored and stored in clim_cv_preds_list for later blending.
  # handle_one() uses <<- to write back into the enclosing for-loop scope.
  clim_logit_tasks <- spec$extras$clim_logits %||% list()
  if (length(clim_logit_tasks) > 0) {
    for (task in clim_logit_tasks) {
      nm <- task$name %||% stop("extras.clim_logits task missing 'name'")
      base_prefix <- task$base_col_prefix %||% stop("extras.clim_logits task missing 'base_col_prefix'")
      
      # Optional window variants for clim_raw
      win_starts <- as.integer(unlist(task$window_start_years %||% integer(0)))
      win_end <- task$window_end_year %||% NA_integer_
      
      # Helper to generate ONE model + store it for blending + score it
      handle_one <- function(model_nm, start_year = NULL, end_year = NULL) {
        
        # unc_clim_raw: include cv_earlier (ONLY here)
        if (identical(model_nm, "unc_clim_raw")) {
          earlier_col <- task$earlier_col %||% stop("unc_clim_raw requires extras.clim_logits.earlier_col in YAML")
          cvp <- make_raw_preds_from_wide_logit(
            wide_df = wide_df,
            base_col_prefix = base_prefix,
            holdout_years = holdout_years,
            earlier_col = earlier_col,
            earlier_is_logit = isTRUE(task$earlier_is_logit %||% TRUE),
            add_cv_earlier = TRUE,
            renormalize_6 = TRUE
          )
          
        } else {
          # clim_raw (and its window variants): 5-bin only
          cvp <- make_raw_preds_from_wide_logit_window(
            wide_df = wide_df,
            base_col_prefix = base_prefix,
            holdout_years = holdout_years,
            start_year = start_year, end_year = end_year
          )
        }
        
        if (is.null(cvp)) return(invisible(NULL))
        
        clim_cv_preds_list[[model_nm]] <<- cvp  # <<- needed: write back to enclosing for-loop scope from handle_one()
        
        res <- save_reliability_and_metrics(
          cvp, model_nm, method, path_box,
          cv_metrics_list, yearly_metrics,
          allowed_cells = dissemination_cells, run_tag = output_tag
        )
        cv_metrics_list <<- res$cv_metrics_list  # <<- needed: accumulate results in enclosing scope
        yearly_metrics  <<- res$yearly_metrics
        
        invisible(NULL)
      }
      
      if (length(win_starts) == 0 || is.na(win_end) || !identical(nm, "clim_raw")) {
        # no windows OR not clim_raw -> just do the base model name
        handle_one(nm, start_year = NULL, end_year = NULL)
      } else {
        # clim_raw window variants
        for (sy in win_starts) {
          model_nm <- if (sy == 1900) nm else paste0(nm, "_", sy, "_", win_end)
          handle_one(model_nm, start_year = sy, end_year = win_end)
        }
      }
    }
  }
  
  # ---- Section 3: Forecast tasks (raw + calibrated) ----
  # Processes NWP forecast sources (NGCM, AIFS) defined in spec$extras$forecasts.
  # Two modes per forecast:
  #   - raw: Direct model probabilities, no calibration.  Optionally computes fair
  #     Brier score (finite-sample corrected).
  #   - calibrated: Platt-calibrated (one-vs-rest logistic) then renormalized to sum to 1.
  #     Calibrated preds are stored by variant name (mok / clim_mok_date / no_mok_filter)
  #     for use in the blending step below.  Platt weights can optionally be exported to
  #     RDS for reuse in the 2025 evaluation script.
  forecast_tasks <- spec$extras$forecasts %||% list()
  cal_preds_store <- list()

  if (length(forecast_tasks) > 0) {
    for (task in forecast_tasks) {
      f_name   <- task$name   %||% stop("extras.forecasts task missing 'name'")
      variant  <- task$variant %||% "base"
      
      do_raw      <- isTRUE(task$raw %||% FALSE)
      do_calibrated <- isTRUE(task$calibrated %||% FALSE)
      
      base_label <- forecast_label(f_name, variant)
      
      # -----------------
      # RAW
      # -----------------
      if (do_raw) {
        cv_preds_raw <- make_raw_preds_from_wide(
          wide_df        = wide_df,
          forecast_name  = f_name,
          variant        = variant,
          holdout_years  = holdout_years,
          spec           = spec
        )
        
        if (!is.null(cv_preds_raw)) {
          model_nm <- paste0(base_label, "_raw")
          
          res <- save_reliability_and_metrics(
            cv_preds_raw, model_nm, method, path_box,
            cv_metrics_list, yearly_metrics,
            allowed_cells = dissemination_cells, run_tag = output_tag
          )
          cv_metrics_list <- res$cv_metrics_list
          yearly_metrics  <- res$yearly_metrics
          
          if (isTRUE(task$fair_brier %||% FALSE)) {
            fb <- compute_fair_brier5(cv_preds_raw, dissemination_cells)
            cat("Raw fair Brier (", model_nm, ", ", method, ", pooled): ", fb, "\n", sep = "")
            
            fair_file <- file.path(path_box, paste0("fair_brier_", f_name, output_tag, ".rds"))
            fair_row  <- tibble::tibble(model = model_nm, cv_method = method, fair_brier = fb)
            
            if (file.exists(fair_file)) {
              old <- readRDS(fair_file)
              fair_row <- dplyr::bind_rows(old, fair_row) %>%
                dplyr::distinct(model, cv_method, .keep_all = TRUE)
            }
            saveRDS(fair_row, fair_file)
          }
        }
      }
      
      # -----------------
      # JUST_CAL
      # -----------------
      if (do_calibrated) {
        cv_preds_cal <- make_calibrated_preds_from_wide(
          wide_df            = wide_df,
          forecast_name      = f_name,
          variant            = variant,
          training_years     = training_years,
          holdout_years      = holdout_years,
          true_holdout_years = true_holdout_years,
          allowed_cells      = dissemination_cells,
          spec               = spec
        )
        # Store calibrated preds for downstream MME blending
        key <- if (identical(variant, "base")) f_name else paste0(f_name, "_calibrated_", variant)
        cal_preds_store[[key]] <- cv_preds_cal
        
        if (!is.null(cv_preds_cal)) {
          model_nm <- if (identical(variant, "base")) paste0(f_name, "_calibrated") else paste0(f_name, "_calibrated_", variant)
          
          res <- save_reliability_and_metrics(
            cv_preds_cal, model_nm, method, path_box,
            cv_metrics_list, yearly_metrics,
            allowed_cells = dissemination_cells, run_tag = output_tag
          )
          cv_metrics_list <- res$cv_metrics_list
          yearly_metrics  <- res$yearly_metrics
        }
        
        # Export reusable platt weights (training years only)
        if (isTRUE(task$export_platt_weights %||% FALSE)) {
          # reuse the same df we used inside make_calibrated... but we need it here for export
          suf  <- get_forecast_variant_suffix(spec, variant)
          cols <- forecast_prob_cols(f_name, suf)
          
          missing <- setdiff(unname(cols), names(wide_df))
          if (length(missing) == 0) {
            df_export <- wide_df %>%
              dplyr::filter(.data$year %in% c(training_years, holdout_years)) %>%
              dplyr::transmute(
                outcome, year,
                week1 = .data[[cols["week1"]]],
                week2 = .data[[cols["week2"]]],
                week3 = .data[[cols["week3"]]],
                week4 = .data[[cols["week4"]]],
                later = pmax(0, 1 - (week1 + week2 + week3 + week4))
              )
            
            platt_fit <- fit_platt_weights_export(
              df = df_export,
              prob_cols = c("week1","week2","week3","week4","later"),
              training_years = training_years,
              outcome_col = "outcome",
              year_col = "year"
            )
            
            saveRDS(
              platt_fit$weights_list,
              file.path(path_box, paste0("platt_weights_", base_label, "_calibrated_list", output_tag, ".rds"))
            )
            saveRDS(
              platt_fit$weights_df,
              file.path(path_box, paste0("platt_weights_", base_label, "_calibrated_df", output_tag, ".rds"))
            )
          }
        }
      }
    }
  }
  
  # ---------------------------------------------------------
  # Section 4: MME optimization
  # ---------------------------------------------------------
  # Finds optimal convex-combination weights for an N-source multimodel ensemble:
  #   P_mix = sum_i w_i * P_i, where weights live on the N-simplex.
  #
  # Sources are specified via spec$mme$blend_models, which lists each model's
  # name and source type (clim or forecast with cal_variant).
  #
  # Strategy: coarse grid search (0.1-step) over the simplex to find a good
  # starting point, then refine with BFGS in unconstrained softmax space.
  # Objective is pooled RPS over dissemination cells.
  #
  # Write vs read path logic: when do_opt_this_run is TRUE (full CV run with
  # holdout_years = 2000:2024), we optimize and write weights.  Otherwise
  # (e.g. the 1965-1978 holdout period), we read previously saved weights via weights_year_tag.
  if (isTRUE(spec$mme$enabled %||% FALSE) && length(clim_cv_preds_list) > 0) {

    mme_variants <- unlist(spec$mme$variants %||% c("clim_mok_date"))
    mc_cores <- getOption("mc.cores", spec$mme$mc_cores %||% 25)
    blend_models <- spec$mme$blend_models %||% list()

    bins5 <- c("week1","week2","week3","week4","later")
    cols5 <- paste0("cv_", bins5)

    n_blend <- length(blend_models)
    blend_names <- vapply(blend_models, function(bm) bm$name, character(1))

    # Helper: generate all points on the N-simplex with given step size
    simplex_grid <- function(n, step = 0.1) {
      if (n == 1) return(matrix(1, nrow = 1, ncol = 1))
      vals <- seq(0, 1, by = step)
      grid <- expand.grid(replicate(n - 1, vals, simplify = FALSE))
      grid <- grid[rowSums(grid) <= 1 + 1e-9, , drop = FALSE]
      grid[[n]] <- 1 - rowSums(grid[, 1:(n-1), drop = FALSE])
      as.matrix(grid)
    }

    # Softmax for N dimensions
    softmax_n <- function(theta) {
      e <- exp(theta - max(theta))
      e / sum(e)
    }

    for (mme_variant in mme_variants) {

      cat("\n=== MME optimization (RPS-only) variant:", mme_variant, "===\n")

      # Resolve sources from blend_models
      mme_sources <- list()
      for (bm in blend_models) {
        if (identical(bm$source, "clim")) {
          src <- clim_cv_preds_list[[bm$name]]
          if (is.null(src)) stop("MME blend_model '", bm$name, "' (source: clim) not found in clim_cv_preds_list.")
        } else if (identical(bm$source, "forecast")) {
          key <- if (identical(bm$cal_variant, "base")) bm$name else paste0(bm$name, "_calibrated_", bm$cal_variant)
          src <- cal_preds_store[[key]]
          if (is.null(src)) stop("MME blend_model '", bm$name, "' (source: forecast, key: ", key, ") not found in cal_preds_store.")
        } else {
          stop("Unknown MME blend_model source: ", bm$source)
        }
        mme_sources[[bm$name]] <- src
      }

      id_vars <- c("time","id","lat","lon","year","outcome")

      # Join all sources into a single data frame
      join_list <- lapply(blend_names, function(nm) {
        mme_sources[[nm]] %>%
          dplyr::select(dplyr::all_of(id_vars), dplyr::starts_with("cv_")) %>%
          dplyr::rename_with(~ gsub("^cv_", paste0("cv_", nm, "_"), .x), dplyr::starts_with("cv_"))
      })

      mme_base_full <- join_list[[1]]
      if (n_blend > 1) {
        for (j in 2:n_blend) {
          mme_base_full <- dplyr::inner_join(mme_base_full, join_list[[j]], by = id_vars)
        }
      }

      mme_base_allowed <- restrict_to_allowed(mme_base_full, dissemination_cells)

      # Build P matrices for each source (full and allowed)
      P_full_list <- lapply(blend_names, function(nm) {
        sel_cols <- paste0("cv_", nm, "_", bins5)
        as.matrix(mme_base_full[, sel_cols, with = FALSE])
      })
      names(P_full_list) <- blend_names

      P_allowed_list <- lapply(blend_names, function(nm) {
        sel_cols <- paste0("cv_", nm, "_", bins5)
        as.matrix(mme_base_allowed[, sel_cols, with = FALSE])
      })
      names(P_allowed_list) <- blend_names

      Y_allowed <- vapply(
        bins5,
        function(b) as.integer(mme_base_allowed$outcome == b),
        numeric(nrow(mme_base_allowed))
      )
      Y_allowed <- matrix(Y_allowed, nrow = nrow(mme_base_allowed), ncol = length(bins5), byrow = FALSE)

      P_allowed_clean <- lapply(P_allowed_list, clean_probs5)

      fast_rps_from_w <- function(w) {
        P_mix <- Reduce("+", mapply(function(wi, Pi) wi * Pi, w, P_allowed_clean, SIMPLIFY = FALSE))
        P_mix <- clean_probs5(P_mix)
        pooled_rps5(P_mix, Y_allowed)
      }

      # Use first clim source name for file naming (backward compat)
      clim_model <- blend_names[vapply(blend_models, function(bm) identical(bm$source, "clim"), logical(1))][1]
      if (is.na(clim_model)) clim_model <- blend_names[1]
      mme_prefix <- paste0("mme_", mme_variant, "_", clim_model, "_opt")

      do_opt_this_run <- isTRUE(spec$mme$optimize_if_holdout_is_full_2000_2024 %||% TRUE) &&
        identical(holdout_years, 2000:2024)

      # Write path always uses the current run's year_tag
      mme_opt_file_write <- file.path(
        path_box,
        paste0("mme_optimized_weights", cutoff_tag, "_", clim_model, "_", mme_variant, year_tag, ".rds")
      )

      # Read path can be overridden via spec$mme$weights_year_tag
      # (e.g., hindcast reads weights saved by the CV run)
      weights_read_tag <- spec$mme$weights_year_tag %||% year_tag
      mme_opt_file_read <- file.path(
        path_box,
        paste0("mme_optimized_weights", cutoff_tag, "_", clim_model, "_", mme_variant, weights_read_tag, ".rds")
      )

      w_col_names <- paste0("w_", blend_names)

      get_saved_mme_weights_rps <- function(file_path) {
        if (!file.exists(file_path)) stop("MME weights file not found: ", file_path)
        tbl <- readRDS(file_path)
        row <- tbl %>% dplyr::filter(.data$objective == "rps") %>% dplyr::slice(1)
        if (nrow(row) == 0) stop("Saved weights file has no row for objective='rps'")
        w <- as.numeric(row[1, w_col_names])
        names(w) <- w_col_names
        w
      }

      if (do_opt_this_run) {
        # Coarse 0.1-step grid search over the N-simplex to initialize BFGS
        weight_mat <- simplex_grid(n_blend, step = 0.1)
        colnames(weight_mat) <- w_col_names

        mme_list <- parallel::mclapply(
          X = seq_len(nrow(weight_mat)),
          FUN = function(i) {
            w <- weight_mat[i, ]
            row <- as.list(w)
            row$rps <- fast_rps_from_w(w)
            tibble::as_tibble(row)
          },
          mc.cores = mc_cores
        )
        mme_results <- dplyr::bind_rows(mme_list)

        best_rps_grid <- mme_results %>%
          dplyr::filter(!is.na(rps)) %>%
          dplyr::arrange(rps) %>%
          dplyr::slice(1)

        w0 <- as.numeric(best_rps_grid[1, w_col_names])
        w0 <- pmax(w0, 1e-6); w0 <- w0 / sum(w0)
        theta0 <- log(w0)

        opt <- optim(
          par = theta0,
          fn  = function(theta) fast_rps_from_w(softmax_n(theta)),
          method = "BFGS",
          control = list(maxit = 200)
        )
        w_star <- softmax_n(opt$par)

        mme_opt_tbl <- tibble::tibble(objective = "rps")
        for (k in seq_along(w_col_names)) {
          mme_opt_tbl[[w_col_names[k]]] <- w_star[k]
        }
        mme_opt_tbl$rps <- fast_rps_from_w(w_star)

        saveRDS(mme_opt_tbl, mme_opt_file_write)

        opt_w <- w_star
        opt_rps <- mme_opt_tbl$rps
      } else {
        cat("\n=== MME comparison: skipping optimization (reading saved RPS weights) ===\n")
        opt_w <- get_saved_mme_weights_rps(mme_opt_file_read)
        opt_rps <- fast_rps_from_w(opt_w)
      }

      cat("MME weights: ", paste(round(opt_w, 4), collapse = ", "), "\n")
      cat("RPS (pooled, allowed): ", opt_rps, "\n")

      P_mix_full <- Reduce("+", mapply(
        function(wi, Pi) wi * clean_probs5(Pi),
        opt_w, P_full_list, SIMPLIFY = FALSE
      ))
      P_mix_full <- clean_probs5(P_mix_full)

      cv_preds_mme_opt_rps <- dplyr::mutate(
        mme_base_full[, ..id_vars],
        cv_week1 = P_mix_full[,1],
        cv_week2 = P_mix_full[,2],
        cv_week3 = P_mix_full[,3],
        cv_week4 = P_mix_full[,4],
        cv_later = P_mix_full[,5]
      )

      res <- save_reliability_and_metrics(
        cv_preds_mme_opt_rps,
        paste0(mme_prefix, "_rps"),
        method, path_box,
        cv_metrics_list, yearly_metrics,
        allowed_cells = dissemination_cells,
        run_tag = output_tag
      )
      cv_metrics_list <- res$cv_metrics_list
      yearly_metrics  <- res$yearly_metrics
    }
  }
}

# -------------------------------------------------------------------
# Section 5: Year-by-year metrics table
# -------------------------------------------------------------------
# Computes per-year Brier skill and RPSS relative to unc_clim_raw (static
# climatology baseline).  Only the global CV method is used.  Output is saved
# as both RDS and CSV for the yearly performance figure (3_produce_figures.R).
if ("global" %in% names(yearly_metrics)) {
  yearly_global <- yearly_metrics[["global"]] %>% purrr::compact() %>% dplyr::bind_rows()
  
  keep_models <- c(
    "blended_model",
    "unc_clim_raw",
    "clim_raw",
    "ngcm_raw",
    "ngcm_clim_mok_date_raw",
    "ngcm_calibrated_clim_mok_date",
    "ngcm_calibrated"
  )
  
  yearly_global <- yearly_global %>% dplyr::filter(model %in% keep_models)
  
  baseline_yearly <- yearly_global %>%
    dplyr::filter(model == "unc_clim_raw") %>%
    dplyr::select(year, brier_ref = brier, rps_ref = rps, auc_ref = auc)
  
  yearly_metrics_global <- yearly_global %>%
    dplyr::left_join(baseline_yearly, by = "year") %>%
    dplyr::mutate(
      brier_skill = 1 - (brier / brier_ref),
      rpss        = 1 - (rps   / rps_ref)
    ) %>%
    dplyr::select(year, model, brier, brier_skill, rps, rpss, auc) %>%
    dplyr::arrange(year, model)
  
  saveRDS(yearly_metrics_global, file.path(path_box, paste0("yearly_metrics_global", output_tag, ".rds")))
  readr::write_csv(yearly_metrics_global, file.path(path_box, paste0("yearly_metrics_global", output_tag, ".csv")))
}

# -------------------------------------------------------------------
# Section 6: Summaries & map inputs
# -------------------------------------------------------------------
# Aggregates per-cell metrics into pooled model summaries (with skill scores
# relative to unc_clim_raw) and per-cell comparison tables for map plotting.
# Builds polygon geometries for the grid cells and packages everything into
# a single map_inputs RDS consumed by 3_produce_figures.R.
summary_maps <- list()
summary_models <- list()

for (method in cv_methods) {
  all_cells <- dplyr::bind_rows(cv_metrics_list[[method]])
  
  model_avg <- summarize_models_pooled(all_cells, baseline_model = "unc_clim_raw")
  
  summary_models[[method]] <- model_avg
  
  cat("For CV method", method, "best model is",
      model_avg %>% dplyr::filter(brier == min(brier, na.rm = TRUE)) %>% dplyr::pull(model),
      "\n")
  
  summary_maps[[method]] <- summarize_maps_compare(
    all_cells = all_cells,
    method = method,
    clim_model = "clim_raw",
    final_model = "blended_model"
  )
  
  # Save per-cell metrics for mapping in 3_produce_figures.R
  cell_fn <- paste0("cell_metrics_", min(holdout_years), "_", max(holdout_years), cutoff_tag, ".rds")
  saveRDS(all_cells, file.path(path_box, cell_fn))
}

# Save summaries
if (MR) {
  if (length(holdout_years) == 1) {
    fn <- paste0("summary_models_", holdout_years, cutoff_tag, ".rds")
  } else {
    fn <- paste0("summary_models_", min(holdout_years), "_", max(holdout_years), cutoff_tag, ".rds")
  }
  saveRDS(summary_models$global, file.path(path_box, fn))
} else {
  saveRDS(summary_models$global, file.path(path_box, "summary_models_q.rds"))
}

# -------------------------------------------------------------------
# Build & save map inputs for separate mapping script
# -------------------------------------------------------------------
polys <- build_polygons_for_mapping(
  grid_centers = wide_df %>% dplyr::distinct(id, lat, lon),
  allowed_cells = dissemination_cells
)

map_inputs <- list(
  polygons_df = polys$polygons_df,
  allowed_polygons_df = polys$allowed_polygons_df,
  india_boundary = india_boundary,
  summary_maps = summary_maps,
  output_tag = output_tag,
  MR = MR
)

saveRDS(map_inputs, file.path(path_box, paste0("map_inputs", output_tag, ".rds")))

tictoc::toc()
cat("Saved map inputs: ", file.path(path_box, paste0("map_inputs", output_tag, ".rds")), "\n", sep = "")
