# ==============================================================================
# File: blend_evaluation_utils.R
# ==============================================================================
# Purpose
#   Utility functions for the 2025 weekly-bin blending evaluation pipeline.
#   Provides spec-driven formula building, raw/calibrated prediction helpers,
#   Platt scaling, blend weight optimization, cell polygon building, reliability
#   chart construction, per-cell and pooled metric computation, and mapping
#   utilities.
#
# Function index (key functions)
#   get_forecast_variant_suffix(spec, variant)
#     Map variant name ("base","mok","clim_mok_date") to column suffix string.
#
#   forecast_prob_cols(forecast_name, variant_suffix)
#     Build named vector of week1..week4 probability column names.
#
#   make_raw_preds_from_wide(wide_df, forecast_name, variant, holdout_years, spec)
#     Extract raw model predictions for holdout years (no calibration).
#
#   make_calibrated_preds_from_wide(wide_df, forecast_name, variant,
#       training_years, holdout_years, true_holdout_years, allowed_cells, spec)
#     Platt-calibrated predictions via platt_cv_multibin.
#
#   build_formulas_from_spec(spec, cutoff_mode)
#     Construct named list of R formulas from spec$models$formulas,
#     with optional windowed variants.
#
#   compute_cv_global(formula, data_train, holdout_years, ...)
#     Leave-one-year-out CV using nnet::multinom, pooling all cells.
#
#   compute_cv_local(formula, data, holdout_years, ...)
#     Per-cell leave-one-year-out CV.
#
#   compute_cv_neighbors(formula, data, holdout_years, ...)
#     Neighbor-cell CV (train on cells within +/-2 lat/lon).
#
#   compute_cv_clusters(formula, data, holdout_years, cluster_var, ...)
#     Cluster-based CV (train on cells in the same cluster).
#
#   compute_cell_metrics_fast(df, multinomial, allowed_cells, ...)
#     Compute per-cell and pooled Brier/RPS/AUC from cv_* probability columns.
#
#   save_reliability_and_metrics(cv_preds, model_name, method, path_box,
#       cv_metrics_list, yearly_metrics, ...)
#     Calibration plot + cell metrics + yearly metrics in one call.
#
#   platt_cv_multibin(df, prob_cols, holdout_years, ...)
#     Cross-validated per-bin Platt calibration (glm logistic).
#
#   summarize_models_pooled(all_cells, baseline_model)
#     Compute skill scores vs a baseline from the "ALL" row.
#
#   summarize_maps_compare(all_cells, method, clim_model, final_model)
#     Per-cell comparison of a model vs climatology.
#
#   clean_probs5(P)
#     Sanitize 5-bin probability matrix (non-finite -> 0, row-normalize).
#
#   pooled_rps5(P, y_onehot)
#     Pooled RPS from 5-bin probability and one-hot matrices.
#
# Dependencies
#   - R/_shared/misc.R
#   - dplyr, data.table, matrixStats, future.apply, purrr, tibble, tidyr,
#     ggplot2, nnet, parallel
# ==============================================================================

source("R/_shared/misc.R")

get_forecast_variant_suffix <- function(spec, variant) {
  m <- spec$extras$forecast_variants %||% list(base = "", mok = "_mok", clim_mok_date = "_clim_mok_date")
  if (!variant %in% names(m)) stop("Unknown forecast variant: ", variant)
  m[[variant]]
}

forecast_prob_cols <- function(forecast_name, variant_suffix) {
  # expects columns like: {name}_p_onset{suffix}_week1 ... week4
  base <- paste0(forecast_name, "_p_onset", variant_suffix)
  c(
    week1 = paste0(base, "_week1"),
    week2 = paste0(base, "_week2"),
    week3 = paste0(base, "_week3"),
    week4 = paste0(base, "_week4")
  )
}

make_raw_preds_from_wide <- function(wide_df, forecast_name, variant, holdout_years, spec) {
  suf <- get_forecast_variant_suffix(spec, variant)
  cols <- forecast_prob_cols(forecast_name, suf)
  
  missing <- setdiff(unname(cols), names(wide_df))
  if (length(missing) > 0) {
    warning("Skipping raw for ", forecast_name, " (", variant, "); missing cols: ", paste(missing, collapse = ", "))
    return(NULL)
  }
  
  wide_df %>%
    dplyr::filter(.data$year %in% holdout_years) %>%
    dplyr::transmute(
      outcome, time, id, lat, lon, year,
      cv_week1 = .data[[cols["week1"]]],
      cv_week2 = .data[[cols["week2"]]],
      cv_week3 = .data[[cols["week3"]]],
      cv_week4 = .data[[cols["week4"]]],
      cv_later = pmin(1, pmax(0, 1 - (
        .data[[cols["week1"]]] + .data[[cols["week2"]]] + .data[[cols["week3"]]] + .data[[cols["week4"]]]
      )))
    )
}

make_calibrated_preds_from_wide <- function(wide_df, forecast_name, variant,
                                          training_years, holdout_years, true_holdout_years,
                                          allowed_cells, spec) {
  suf <- get_forecast_variant_suffix(spec, variant)
  cols <- forecast_prob_cols(forecast_name, suf)
  
  missing <- setdiff(unname(cols), names(wide_df))
  if (length(missing) > 0) {
    warning("Skipping calibrated for ", forecast_name, " (", variant, "); missing cols: ", paste(missing, collapse = ", "))
    return(NULL)
  }
  
  df <- wide_df %>%
    dplyr::filter(.data$year %in% c(training_years, holdout_years)) %>%
    dplyr::transmute(
      outcome, time, id, lat, lon, year,
      week1 = .data[[cols["week1"]]],
      week2 = .data[[cols["week2"]]],
      week3 = .data[[cols["week3"]]],
      week4 = .data[[cols["week4"]]],
      later = pmax(0, 1 - (week1 + week2 + week3 + week4))
    )
  
  platt_cv_multibin(
    df                 = df,
    prob_cols          = c("week1","week2","week3","week4","later"),
    holdout_years      = holdout_years,
    true_holdout_years = true_holdout_years,
    outcome_col        = "outcome",
    year_col           = "year",
    cv_prefix          = "cv_",
    allowed_cells      = allowed_cells
  )
}

forecast_label <- function(forecast_name, variant) {
  # determines filenames/model names: ngcm_raw, ngcm_clim_mok_date_raw, ngcm_calibrated_clim_mok_date, etc.
  if (identical(variant, "base")) forecast_name else paste0(forecast_name, "_", variant)
}


# Year tag for filenames.
# years: integer vector; returns "", "_2024", or "_1965_1978".
make_year_tag <- function(years) {
  years <- sort(unique(years))
  if (length(years) == 0) return("")
  if (length(years) == 1) return(paste0("_", years))
  paste0("_", min(years), "_", max(years))
}

# Cutoff tag for filenames.
# cutoff_mode: "mok", "clim_mok_date", or "no_mok_filter"; returns "", "_clim_mok_date", or "_no_mok_filter".
make_cutoff_tag <- function(cutoff_mode) {
  dplyr::case_when(
    identical(cutoff_mode, "clim_mok_date")  ~ "_clim_mok_date",
    identical(cutoff_mode, "no_mok_filter")  ~ "_no_mok_filter",
    TRUE                                     ~ ""
  )
}

# Expand formula shortcuts in a formula string.
# - RHS terms containing "_qx" expand to "_week1"..."_week4".
# formula_str: string with optional LHS; returns an R formula object.
expand_formula <- function(formula_str) {
  parts <- strsplit(formula_str, "~", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    lhs <- trimws(parts[1])
    rhs <- trimws(parts[2])
  } else {
    lhs <- ""
    rhs <- trimws(parts[1])
  }

  terms <- strsplit(rhs, "\\+")[[1]]
  terms <- trimws(terms)

  expanded_terms <- sapply(terms, function(term) {
    expanded <- term

    # Expand _qx -> _week1.._week4
    new_expanded <- c()
    for (candidate in expanded) {
      if (grepl("_qx", candidate, fixed = TRUE)) {
        for (i in 1:4) {
          new_expanded <- c(new_expanded, gsub("_qx", paste0("_week", i), candidate, fixed = TRUE))
        }
      } else {
        new_expanded <- c(new_expanded, candidate)
      }
    }
    expanded <- new_expanded

    if (length(expanded) > 1) paste0("(", paste(expanded, collapse = " + "), ")") else expanded
  })
  
  new_rhs <- paste(expanded_terms, collapse = " + ")
  new_formula_str <- if (lhs != "") paste(lhs, "~", new_rhs) else new_rhs
  as.formula(new_formula_str)
}

# Window suffix string like "_2000_2024".
make_window_suffix <- function(start_year, end_year) paste0("_", start_year, "_", end_year)

make_raw_preds_from_wide_logit_window <- function(wide_df, base_col_prefix,
                                                  holdout_years,
                                                  start_year = NULL, end_year = NULL) {
  window <- if (!is.null(start_year) && !is.null(end_year) && start_year != 1900) {
    paste0("_", start_year, "_", end_year)
  } else {
    ""  # baseline/unwindowed (incl 1900 convention)
  }
  prefix <- paste0(base_col_prefix, window)
  
  wk_cols <- paste0(prefix, "_week", 1:4)
  missing <- setdiff(wk_cols, names(wide_df))
  if (length(missing) > 0) {
    warning("Skipping clim logit raw for prefix=", prefix, "; missing cols: ",
            paste(missing, collapse = ", "))
    return(NULL)
  }
  
  wide_df %>%
    dplyr::filter(.data$year %in% holdout_years) %>%
    dplyr::mutate(
      p1 = plogis(.data[[wk_cols[1]]]),
      p2 = plogis(.data[[wk_cols[2]]]),
      p3 = plogis(.data[[wk_cols[3]]]),
      p4 = plogis(.data[[wk_cols[4]]]),
      pL = pmin(1, pmax(0, 1 - (p1 + p2 + p3 + p4)))
    ) %>%
    dplyr::transmute(
      outcome, time, id, lat, lon, year,
      cv_week1 = p1, cv_week2 = p2, cv_week3 = p3, cv_week4 = p4, cv_later = pL
    )
}

make_raw_preds_from_wide_logit <- function(wide_df,
                                           base_col_prefix,
                                           holdout_years,
                                           earlier_col,
                                           earlier_is_logit = TRUE,
                                           add_cv_earlier = TRUE,
                                           renormalize_6 = TRUE) {
  
  wk_cols <- paste0(base_col_prefix, "_week", 1:4)
  missing <- setdiff(wk_cols, names(wide_df))
  if (length(missing) > 0) {
    warning("Skipping raw-logit for ", base_col_prefix, "; missing cols: ",
            paste(missing, collapse = ", "))
    return(NULL)
  }
  
  if (isTRUE(add_cv_earlier) && (is.null(earlier_col) || !earlier_col %in% names(wide_df))) {
    stop("unc_clim_raw requires earlier_col and it must exist in wide_df.")
  }
  
  out <- wide_df %>%
    dplyr::filter(.data$year %in% holdout_years) %>%
    dplyr::mutate(
      p1 = plogis(.data[[wk_cols[1]]]),
      p2 = plogis(.data[[wk_cols[2]]]),
      p3 = plogis(.data[[wk_cols[3]]]),
      p4 = plogis(.data[[wk_cols[4]]]),
      pE = dplyr::case_when(
        isTRUE(add_cv_earlier) & earlier_is_logit ~ plogis(.data[[earlier_col]]),
        isTRUE(add_cv_earlier)                   ~ as.numeric(.data[[earlier_col]]),
        TRUE                                     ~ 0
      ),
      pE = pmin(1, pmax(0, pE)),
      # later-only residual BEFORE renormalizing
      pL = 1 - (pE + p1 + p2 + p3 + p4),
      pL = pmin(1, pmax(0, pL))
    )
  
  if (isTRUE(renormalize_6) && isTRUE(add_cv_earlier)) {
    rs6 <- with(out, pE + p1 + p2 + p3 + p4 + pL)
    good <- is.finite(rs6) & rs6 > 0
    out <- out %>%
      dplyr::mutate(
        rs6 = ifelse(good, rs6, NA_real_),
        pE = pE / rs6,
        p1 = p1 / rs6,
        p2 = p2 / rs6,
        p3 = p3 / rs6,
        p4 = p4 / rs6,
        pL = pL / rs6
      )
  }
  
  out %>%
    dplyr::transmute(
      outcome, time, id, lat, lon, year,
      cv_week1 = p1,
      cv_week2 = p2,
      cv_week3 = p3,
      cv_week4 = p4,
      cv_later = pL,
      cv_earlier = if (isTRUE(add_cv_earlier)) pE else NULL
    )
}

# Clean 5-bin probability matrix:
# - replaces non-finite with 0
# - truncates negatives to 0
# - ensures each row sums to 1 (uniform if row sum is 0).
clean_probs5 <- function(P) {
  P <- as.matrix(P)
  P[!is.finite(P)] <- 0
  P <- pmax(P, 0)
  rs <- rowSums(P)
  zero <- rs <= 0
  if (any(zero)) {
    P[zero, ] <- 1 / ncol(P)
    rs <- rowSums(P)
  }
  P / rs
}

# Pooled ranked probability score (RPS) for 5 bins.
# P: n x 5 prob matrix; y_onehot: n x 5 one-hot outcomes; returns scalar.
pooled_rps5 <- function(P, y_onehot) {
  CP <- t(apply(P, 1, cumsum))
  CY <- t(apply(y_onehot, 1, cumsum))
  mean(rowSums((CP - CY)^2))
}

# Fair Brier score for 5-category CV predictions (finite-sample correction).
# cv_preds must have outcome and cv_week1..cv_week4..cv_later.
# allowed_cells used via restrict_to_allowed(); returns scalar.
compute_fair_brier5 <- function(cv_preds, allowed_cells) {
  bins <- c("week1", "week2", "week3", "week4", "later")
  cols <- paste0("cv_", bins)
  
  restrict_to_allowed(cv_preds, allowed_cells) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      fair_brier_row = {
        p <- dplyr::c_across(dplyr::all_of(cols))
        y <- as.integer(bins == outcome)
        bs_row <- sum((p - y)^2)
        corr_row <- sum(p * (1 - p)) / (30 - 1)
        bs_row - corr_row
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(fair_brier = mean(fair_brier_row, na.rm = TRUE)) %>%
    dplyr::pull(fair_brier)
}

# Fit one-vs-rest Platt calibration weights on training years.
# df must contain outcome_col, year_col, and columns named in prob_cols (raw probs in (0,1)).
# Returns list: weights_list (per-bin coefs) and weights_df (tidy table of intercept/slope).
fit_platt_weights_export <- function(df, prob_cols, training_years,
                                     outcome_col = "outcome", year_col = "year") {
  df_train <- df %>%
    dplyr::filter(.data[[year_col]] %in% training_years)
  
  fit_one_bin <- function(bin_name) {
    p_raw <- pmin(pmax(df_train[[bin_name]], 1e-6), 1 - 1e-6)
    logit_p <- qlogis(p_raw)
    y <- as.integer(df_train[[outcome_col]] == bin_name)
    fit <- glm(y ~ logit_p, family = binomial)
    
    list(bin = bin_name, coef = stats::coef(fit))
  }
  
  weights_list <- lapply(prob_cols, fit_one_bin)
  names(weights_list) <- prob_cols
  
  weights_df <- prob_cols %>%
    purrr::map_dfr(function(bin) {
      cf <- weights_list[[bin]]$coef
      tibble::tibble(
        bin       = bin,
        intercept = unname(cf[1]),
        slope     = unname(cf["logit_p"])
      )
    })
  
  list(weights_list = weights_list, weights_df = weights_df)
}

# Exclude cells for mapping (NE + SW islands rule-of-thumb).
# df must contain lat/lon; returns filtered df.
exclude_for_mapping <- function(df) {
  df %>% dplyr::filter(!(lon > 90 | (lon < 75 & lat < 12)))
}

# Build cell polygons for mapping.
# grid_centers: lat/lon centers; allowed_cells: lat/lon allowed subset.
# Returns list(polygons_df=..., allowed_polygons_df=...).
build_polygons_for_mapping <- function(grid_centers, allowed_cells, half_size = 1) {
  poly_data <- grid_centers %>%
    dplyr::distinct(id, lat, lon) %>%
    exclude_for_mapping()

  polygons_list <- lapply(seq_len(nrow(poly_data)), function(i) {
    cell <- poly_data[i, ]
    poly <- create_cell_polygon(cell$lat, cell$lon, half_size = half_size)
    poly$id <- cell$id
    poly$lat_center <- cell$lat
    poly$lon_center <- cell$lon
    poly
  })

  polygons_df <- do.call(rbind, polygons_list)

  allowed_polygons_df <- polygons_df %>%
    dplyr::inner_join(
      allowed_cells,
      by = "id"
    )

  list(polygons_df = polygons_df, allowed_polygons_df = allowed_polygons_df)
}

# Draw multiple ggplots in one row with a shared title using grid graphics.
# title: string; ...: ggplot objects; heights: c(title_height, plots_height).
arrange_ggplots_1row <- function(title, ..., heights = c(0.08, 0.92)) {
  plots <- list(...)
  n <- length(plots)
  
  grid::grid.newpage()
  lay <- grid::grid.layout(
    nrow = 2, ncol = n,
    heights = grid::unit(heights, c("null", "null"))
  )
  grid::pushViewport(grid::viewport(layout = lay))
  
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1:n))
  grid::grid.text(title, gp = grid::gpar(fontsize = 18, fontface = "bold"))
  grid::popViewport()
  
  for (i in seq_len(n)) {
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = i))
    grid::grid.draw(ggplot2::ggplotGrob(plots[[i]]))
    grid::popViewport()
  }
  
  grid::popViewport()
  invisible(NULL)
}

# Legend title based on metric column name.
legend_title_for_metric <- function(metric_col) {
  if (grepl("skill|rpss", metric_col, ignore.case = TRUE)) "Skill" else "Difference"
}

# Build formulas_list from spec settings
# Reads formula text from spec$models$formulas. Optionally generates windowed
# variants using spec$models$window_variants.
# @return Named list of formulas.
build_formulas_from_spec <- function(spec, cutoff_mode) {
  # ---- Base formulas: entirely from spec ----
  formula_cfg <- spec$models$formulas %||% list()
  
  if (length(formula_cfg) == 0) {
    stop("Spec must define models$formulas with named entries containing $text.")
  }
  
  base_texts <- purrr::imap_chr(formula_cfg, function(v, nm) {
    if (isTRUE(v$enabled %||% TRUE)) {
      txt <- v$text %||% NA_character_
      if (!is.character(txt) || is.na(txt) || !nzchar(txt)) {
        stop("Formula '", nm, "' is enabled but has no non-empty 'text' in spec.")
      }
      return(txt)
    } else {
      return(NA_character_)
    }
  })
  
  base_texts <- base_texts[!is.na(base_texts)]
  
  # Expand shortcuts and convert to formulas
  formulas_list <- base_texts %>%
    purrr::map(expand_formula)
  
  # ---- Optional window variants (spec-driven template) ----
  wcfg_all <- spec$models$window_variants %||% NULL
  
  # allow either a single mapping or a list of mappings
  if (!is.null(wcfg_all)) {
    wcfg_list <- if (is.list(wcfg_all) && !is.null(names(wcfg_all)) && !is.null(wcfg_all$enabled)) {
      list(wcfg_all)
    } else if (is.list(wcfg_all) && is.null(wcfg_all$enabled)) {
      wcfg_all
    } else {
      list(wcfg_all)
    }
    
    for (wcfg in wcfg_list) {
      if (!isTRUE(wcfg$enabled %||% FALSE)) next
      
      only_if <- wcfg$only_if_cutoff_mode %||% NULL
      if (!is.null(only_if) && !identical(cutoff_mode, only_if)) next
      
      base_name <- wcfg$base_name %||% NULL
      if (is.null(base_name) || !base_name %in% names(base_texts)) {
        stop("window_variants$base_name must match a name in models$formulas.")
      }
      
      start_years <- as.integer(unlist(wcfg$start_years %||% integer(0)))
      end_year <- as.integer(wcfg$end_year %||% NA_integer_)
      if (length(start_years) == 0 || is.na(end_year)) {
        stop("window_variants must define start_years and end_year.")
      }
      
      include_baseline_1900 <- isTRUE(wcfg$include_baseline_1900 %||% TRUE)
      
      from <- wcfg$replace$from %||% NULL
      to   <- wcfg$replace$to %||% NULL
      if (is.null(from) || is.null(to)) {
        stop("window_variants$replace must define 'from' and 'to'.")
      }
      
      base_txt <- base_texts[[base_name]]
      
      for (sy in start_years) {
        if (sy == 1900) {
          if (include_baseline_1900) {
            # already present as base_name
          }
          next
        }
        window <- make_window_suffix(sy, end_year)
        windowed_txt <- gsub(from, gsub("\\{window\\}", window, to), base_txt, fixed = TRUE)
        nm <- paste0(base_name, window)
        formulas_list[[nm]] <- expand_formula(windowed_txt)
      }
    }
  }
  
  formulas_list
}

# Input RDS filename by cutoff mode and optional resolution.
input_rds_from_cutoff <- function(cutoff_mode, resolution = "") {
  prefix <- if (nzchar(resolution)) paste0(resolution, "_") else ""
  dplyr::case_when(
    identical(cutoff_mode, "clim_mok_date")  ~ paste0("cv_data_", prefix, "clim_mok_date_new_pipeline.rds"),
    identical(cutoff_mode, "no_mok_filter")  ~ paste0("cv_data_", prefix, "no_mok_filter_new_pipeline.rds"),
    TRUE                                     ~ paste0("cv_data_", prefix, "new_pipeline.rds")
  )
}

# Summarize pooled metrics and compute skills vs a baseline model.
# all_cells: bound cell metrics (expects rows with id == "ALL"); baseline_model used for reference.
# Returns tibble with added columns: brier_skill, rps_skill, and `AUC diff` (percentage points).
summarize_models_pooled <- function(all_cells, baseline_model = "unc_clim_raw") {
  model_avg <- all_cells %>% dplyr::filter(id == "ALL")

  clim_row <- model_avg %>% dplyr::filter(model == baseline_model) %>% dplyr::slice(1)
  
  brier_clim <- clim_row$brier[1]
  rps_clim   <- clim_row$rps[1]
  auc_clim   <- clim_row$auc[1]
  
  model_avg %>%
    dplyr::mutate(
      brier_skill = 1 - (brier / brier_clim),
      rps_skill   = 1 - (rps   / rps_clim),
      `AUC diff`  = (auc - auc_clim) * 100,
      cr_one_diff = (cr_one - clim_row$cr_one[1]) * 100,
      cr_two_diff = (cr_two - clim_row$cr_two[1]) * 100
    )
}

# Per-cell comparison of a final model vs climatology baseline.
# all_cells must contain lat/lon/model and brier/rps/auc columns.
# Returns tibble with per-cell auc_diff, brier_skill, rps_skill, and cv_method.
summarize_maps_compare <- function(all_cells, method, clim_model = "clim_raw", final_model = "blended_model") {
  clim_df <- all_cells %>%
    dplyr::filter(model == clim_model) %>%
    dplyr::rename(brier_clim = brier, rps_clim = rps, auc_clim = auc) %>%
    dplyr::select(id, brier_clim, rps_clim, auc_clim)

  final_df <- all_cells %>%
    dplyr::filter(model == final_model) %>%
    dplyr::rename(brier_final = brier, rps_final = rps, auc_final = auc) %>%
    dplyr::select(id, brier_final, rps_final, auc_final)

  dplyr::full_join(clim_df, final_df, by = "id") %>%
    dplyr::mutate(
      auc_diff    = auc_final - auc_clim,
      brier_skill = 1 - (brier_final / brier_clim),
      rps_skill   = 1 - (rps_final   / rps_clim),
      cv_method   = method
    )
}

# fits the blended model in a cross-validated way 
compute_cv_global <- function(formula,
                              data_train,
                              holdout_years,
                              multinomial = TRUE,
                              true_holdout_years = c(),
                              data_pred = NULL) {
  if (is.null(data_pred)){
    data_pred = data_train
  }
  cv_results <- mclapply(
    mc.cores = getOption("mc.cores", 25L),
    holdout_years,
    function(test_year) {
      
      train <- dplyr::filter(data_train, year != test_year & !(year %in% true_holdout_years))
      test  <- dplyr::filter(data_pred,  year == test_year)
      
      # shapley testing hook (keep as-is)
      if ("train_q_outcome" %in% names(train) && all(train$train_q_outcome, na.rm = TRUE)) {
        train <- train %>% mutate(outcome = q_outcome)
      }
      
      if (nrow(train) == 0 || nrow(test) == 0) return(NULL)
      
      fit <- multinom(formula, data = train, maxit = 1000)
      preds <- predict(fit, newdata = test, type = "probs")

      if (multinomial) {
        colnames(preds) <- paste0("cv_", colnames(preds))
        cbind(test, preds)
      } else {
        cbind(test, cv_week1 = preds)
      }
    }
  )
  
  do.call(rbind, cv_results)
}

compute_cv_local <-
  function(formula, data, holdout_years, multinomial = TRUE) {
    unique_cells <- data %>% dplyr::distinct(id, lat, lon)
    cv_results <-
      mclapply(
        1:nrow(unique_cells),
        function(i) {
          cell <- unique_cells[i, ]
          cell_data <- dplyr::filter(data, id == cell$id)
          cell_holdout_years <- intersect(holdout_years, unique(cell_data$year))
          cell_cv <- lapply(cell_holdout_years, function(test_year) {
            train <- dplyr::filter(cell_data, year != test_year)
            test  <- dplyr::filter(cell_data, year == test_year)
            if(nrow(train) == 0 || nrow(test) == 0) return(NULL)
            if(length(unique(train$outcome) == 1)) return(cbind(test, cv_week1 = NA))
            fit <- multinom(formula, data = train, maxit = 1000)
            preds <- predict(fit, newdata = test, type = "probs")
            if(multinomial){
              if(nrow(test) == 1){
                preds <- as.data.frame(t(preds))
              }
              colnames(preds) <- paste0("cv_", colnames(preds))
              return(cbind(test,preds))
            }else{
              cbind(test, cv_week1 = preds)}
          })
          do.call(rbind, cell_cv)
        })
    do.call(rbind, cv_results)
  }

compute_cv_neighbors <- function(formula, data, holdout_years, multinomial = TRUE) {
  unique_cells <- data %>% dplyr::distinct(id, lat, lon)
  cv_results <- mclapply(1:nrow(unique_cells), function(i) {
    target_cell <- unique_cells[i, ]
    cell_data <- dplyr::filter(data, id == target_cell$id)
    cell_holdout_years <- intersect(holdout_years, unique(cell_data$year))
    cell_cv <- lapply(cell_holdout_years, function(test_year) {
      train <- dplyr::filter(data,
                             year != test_year,
                             abs(lat - target_cell$lat) <= 2,
                             abs(lon - target_cell$lon) <= 2)
      test  <- dplyr::filter(data, year == test_year, id == target_cell$id)
      if(nrow(train) == 0 || nrow(test) == 0) return(NULL)
      if(length(unique(train$outcome)) == 1){
        return(cbind(test, cv_week1 = NA))
      }
      fit <- multinom(formula, data = train, maxit=1000)
      preds <- predict(fit, newdata = test, type = "probs")
      if(multinomial){
        if(nrow(test) == 1){
          preds <- as.data.frame(t(preds))
        }

        colnames(preds) <- paste0("cv_", colnames(preds))
        return(cbind(test,preds))
      }else{
        cbind(test, cv_week1 = preds)}
    })
    do.call(rbind, cell_cv)
  })
  do.call(rbind, cv_results)
}

compute_cv_clusters <- function(formula, data, holdout_years, cluster_var, multinomial = TRUE) {
  unique_cells <- data %>% dplyr::distinct(id, lat, lon)
  cv_results <- mclapply(1:nrow(unique_cells), function(i) {
    target_cell <- unique_cells[i, ]
    cell_data <- dplyr::filter(data, id == target_cell$id)
    cluster = cell_data[[cluster_var]][1]
    cell_holdout_years <- intersect(holdout_years, unique(cell_data$year))
    cell_cv <- lapply(cell_holdout_years, function(test_year) {
      train <- dplyr::filter(data,
                             year != test_year,
                             .data[[cluster_var]] == cluster)
      test  <- dplyr::filter(data, year == test_year, id == target_cell$id)
      if(nrow(train) == 0 || nrow(test) == 0) return(NULL)
      if(length(unique(train$outcome)) == 1){
        return(cbind(test, cv_week1 = NA))
      }
      fit <- multinom(formula, data = train, maxit=1000)
      preds <- predict(fit, newdata = test, type = "probs")
      if(multinomial){
        if(nrow(test) == 1){
          preds <- as.data.frame(t(preds))
        }

        colnames(preds) <- paste0("cv_", colnames(preds))
        return(cbind(test,preds))
      }else{
        cbind(test, cv_week1 = preds)}
    })
    do.call(rbind, cell_cv)
  })
  do.call(rbind, cv_results)
}

create_cell_polygon <- function(lat, lon, half_size = 1) {
  data.frame(
    lon = c(lon - half_size, lon + half_size, lon + half_size, lon - half_size, lon - half_size),
    lat = c(lat - half_size, lat - half_size, lat + half_size, lat + half_size, lat - half_size)
  )
}

platt_cv_multibin <- function(df,
                              prob_cols,
                              holdout_years,
                              true_holdout_years = integer(),
                              outcome_col = "outcome",
                              year_col    = "year",
                              cv_prefix   = "cv_",
                              allowed_cells = NULL) {
  df <- data.table::as.data.table(df)

  bins <- prob_cols

  out_list <- lapply(holdout_years, function(test_year) {

    train <- df[df[[year_col]] != test_year &
                  !(df[[year_col]] %in% true_holdout_years), , drop = FALSE]
    test  <- df[df[[year_col]] == test_year, , drop = FALSE]
    if (nrow(test) == 0) {
      return(NULL)
    }
    if (!is.null(allowed_cells)) {
      train <- restrict_to_allowed(train, allowed_cells)
    }

    cal_mat <- matrix(
      NA_real_,
      nrow = nrow(test),
      ncol = length(bins),
      dimnames = list(NULL, bins)
    )

    for (b in bins) {
      p_tr <- as.numeric(train[[b]])
      y_tr <- as.integer(as.character(train[[outcome_col]]) == b)

      ok <- is.finite(p_tr) & !is.na(y_tr)

      # Skip Platt scaling if fewer than 10 usable training samples or only one
      # class present; fall back to raw probabilities for this bin.
      if (sum(ok) < 10 || length(unique(y_tr[ok])) < 2) {
        cal_mat[, b] <- as.numeric(test[[b]])
        next
      }

      train_df <- data.frame(
        y = y_tr[ok],
        p = p_tr[ok]
      )

      fit <- stats::glm(y ~ p, data = train_df, family = stats::binomial())

      test_p <- as.numeric(test[[b]])
      pred <- stats::predict(
        fit,
        newdata = data.frame(p = test_p),
        type = "response"
      )

      cal_mat[, b] <- as.numeric(pred)
    }

    row_sums <- rowSums(cal_mat, na.rm = TRUE)
    cal_mat <- sweep(cal_mat, 1, row_sums, FUN = "/")
    cal_mat[!is.finite(cal_mat)] <- NA_real_

    meta_cols <- setdiff(names(test), bins)
    res <- cbind(
      test[, ..meta_cols],
      as.data.frame(cal_mat)
    )
    res
  })

  out <- dplyr::bind_rows(out_list)

  for (b in bins) {
    names(out)[names(out) == b] <- paste0(cv_prefix, b)
  }

  out
}

sanitize_filename <- function(x) {
  x %>%
    gsub("[[:space:]]+", "_", .) %>%
    gsub("[^A-Za-z0-9_\\-\\.]", "", .) %>%
    substr(1, 180)
}

calibration_from_wide_probs_quantiles <- function(cv_preds,
                                                  outcome_col = "outcome",
                                                  prob_prefix = "cv_",
                                                  bins = c("earlier","week1","week2","week3","week4","later"),
                                                  nbins = 10) {

  if (!outcome_col %in% names(cv_preds)) {
    stop("Outcome column '", outcome_col, "' not found in cv_preds.")
  }

  prob_cols <- paste0(prob_prefix, bins)
  keep <- prob_cols %in% names(cv_preds)
  bins <- bins[keep]
  prob_cols <- prob_cols[keep]

  if (length(prob_cols) == 0) {
    stop("No probability columns found among: ", paste(paste0(prob_prefix, bins), collapse = ", "))
  }

  long <- cv_preds %>%
    dplyr::select(dplyr::all_of(outcome_col), dplyr::all_of(prob_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(prob_cols),
      names_to = "forecast_bin",
      values_to = "p"
    ) %>%
    dplyr::mutate(
      forecast_bin = sub(paste0("^", prob_prefix), "", forecast_bin),
      outcome_chr  = as.character(.data[[outcome_col]]),
      y = as.integer(outcome_chr == forecast_bin)
    ) %>%
    dplyr::filter(is.finite(p), !is.na(y))

  cal <- long %>%
    dplyr::group_by(forecast_bin) %>%
    dplyr::mutate(qbin = dplyr::ntile(p, nbins)) %>%
    dplyr::group_by(forecast_bin, qbin) %>%
    dplyr::summarise(
      n         = dplyr::n(),
      pred_mean = mean(p, na.rm = TRUE),
      obs_frac  = mean(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(forecast_bin, pred_mean)

  cal
}

restrict_to_allowed <- function(df, allowed_cells) {
  allowed_ids <- allowed_cells %>%
    dplyr::distinct(id)

  df %>%
    dplyr::semi_join(allowed_ids, by = "id")
}

save_reliability_and_metrics <- function(cv_preds, model_name, method,
                                         path_box,
                                         cv_metrics_list, yearly_metrics,
                                         nbins = 10,
                                         prob_prefix = "cv_",
                                         allowed_cells = NULL,
                                         run_tag = "") {
  cv_preds <- data.table::as.data.table(cv_preds)
  tag <- as.character(ifelse(is.null(run_tag), "", run_tag))
  if (nzchar(tag) && !startsWith(tag, "_")) tag <- paste0("_", tag)

  calib_dir <- file.path(path_box, "calibration plots")
  if (!dir.exists(calib_dir)) dir.create(calib_dir, recursive = TRUE)

  bybin_list <- NULL

  tryCatch({
    prob_cols_hist <- paste0(prob_prefix, c("week1","week2","week3","week4","later"))
    prob_cols_hist <- intersect(prob_cols_hist, names(cv_preds))

    preds_hist <- as.numeric(unlist(cv_preds[, ..prob_cols_hist], use.names = FALSE))
    preds_hist <- preds_hist[is.finite(preds_hist) & !is.na(preds_hist)]
    preds_hist <- preds_hist[preds_hist >= 0 & preds_hist <= 1]

    cal_df <- calibration_from_wide_probs_pooled_quantiles(
      cv_preds     = cv_preds,
      outcome_col  = "outcome",
      prob_prefix  = prob_prefix,
      bins         = c("earlier","week1","week2","week3","week4","later"),
      nbins        = nbins
    )

    chart_df <- build_reliability_chart_data(
      cv_preds     = cv_preds,
      cal_df       = cal_df,
      prob_prefix  = prob_prefix,
      bins         = c("earlier","week1","week2","week3","week4","later"),
      hist_bins    = nrow(cal_df)
    )

    stub <- sanitize_filename(sprintf("reliability_quantiles_%s_%s_%dbins%s", method, model_name, nbins, tag))
    chart_rds <- file.path(calib_dir, paste0(stub, "_chartdata.rds"))
    save_reliability_chart_rds(chart_df, chart_rds)

    p_cal <- plot_reliability_from_chart_rds(
      chart_rds,
      title_text = sprintf("Reliability (quantile bins) | model=%s | cv=%s", model_name, method),
      show_diag = TRUE
    )

    png_file <- file.path(
      calib_dir,
      sprintf("reliability_quantiles_%s_%s%s.png", method, model_name, tag)
    )
    ggplot2::ggsave(png_file, p_cal, width = 8, height = 6, dpi = 300)

  }, error = function(e) {
    warning("Calibration plotting failed for model=", model_name, " method=", method, " : ", conditionMessage(e))
    bybin_list <<- NULL
  })

  cell_metrics <- compute_cell_metrics_fast(cv_preds, allowed_cells = allowed_cells)
  cell_metrics$model <- model_name
  cv_metrics_list[[method]][[model_name]] <- cell_metrics

  yearly_metrics[[method]][[model_name]] <- compute_yearly_from_preds(cv_preds, model_name, allowed_cells)

  list(
    cv_metrics_list = cv_metrics_list,
    yearly_metrics  = yearly_metrics
  )
}

fast_auc <- function(y01, score) {
  ok <- is.finite(score) & !is.na(y01)
  y  <- as.integer(y01[ok])
  s  <- score[ok]
  n1 <-  as.numeric(sum(y == 1L))
  n0 <-  as.numeric(sum(y == 0L))
  r <- rank(s, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

compute_cell_metrics_fast <- function(df, multinomial = TRUE, allowed_cells = NULL,
                                      workers = 6L, parallel_auc = TRUE,
                                      print_outputs = TRUE,
                                      only_all = FALSE) {
  stopifnot(is.data.frame(df))
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Please install matrixStats.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install data.table.")
  }
  df <- data.table::as.data.table(df)

  outcome_levels <- c("week1","week2","week3","week4","later")
  desired_rps <- c("earlier", outcome_levels)

  n <- nrow(df)
  if (n == 0L) return(tibble::tibble())

  prob_cols_rps <- paste0("cv_", desired_rps)
  prob_mat <- matrix(0, nrow = n, ncol = length(prob_cols_rps))
  colnames(prob_mat) <- desired_rps

  present <- intersect(prob_cols_rps, names(df))
  if (length(present)) {
    tmp <- as.matrix(df[, ..present])
    storage.mode(tmp) <- "double"
    tmp[!is.finite(tmp)] <- 0
    prob_mat[, sub("^cv_", "", present)] <- tmp
  }

  out_chr <- df$outcome
  if (is.factor(out_chr)) out_chr <- as.character(out_chr)
  out_norm <- tolower(gsub("\\s+", "", out_chr))

  obs_mat <- matrix(0L, nrow = n, ncol = ncol(prob_mat))
  colnames(obs_mat) <- colnames(prob_mat)
  j <- match(out_norm, colnames(prob_mat))
  ok <- which(!is.na(j))
  obs_mat[cbind(ok, j[ok])] <- 1L

  cdf_p <- matrixStats::rowCumsums(prob_mat)
  cdf_o <- matrixStats::rowCumsums(obs_mat)
  cdiff <- cdf_p - cdf_o
  rps_row <- matrixStats::rowSums2(cdiff[, -ncol(cdiff), drop = FALSE]^2)
  df$rps_row <- rps_row

  prob5_cols <- paste0("cv_", outcome_levels)
  prob5 <- matrix(0, nrow = n, ncol = length(outcome_levels))
  colnames(prob5) <- outcome_levels

  present5 <- intersect(prob5_cols, names(df))
  if (length(present5)) {
    tmp5 <- as.matrix(df[, ..present5])
    storage.mode(tmp5) <- "double"
    tmp5[!is.finite(tmp5)] <- 0
    prob5[, sub("^cv_", "", present5)] <- tmp5
  }

  obs5 <- matrix(0L, nrow = n, ncol = length(outcome_levels))
  colnames(obs5) <- outcome_levels
  j5 <- match(out_norm, outcome_levels)
  ok5 <- which(!is.na(j5))
  obs5[cbind(ok5, j5[ok5])] <- 1L

  err5 <- (prob5 - obs5)^2
  brier_row <- matrixStats::rowSums2(err5)
  df$brier_row <- brier_row

  # cr_one: does the highest-probability bin match the true outcome?
  cr_one_row <- as.integer(max.col(prob5, ties.method = "first") == j5)
  cr_one_row[is.na(j5)] <- NA_integer_
  df$cr_one_row <- cr_one_row

  # cr_two: does the best contiguous 2-bin pair contain the true outcome?
  n_bins <- ncol(prob5)
  pair_sums <- prob5[, -n_bins, drop = FALSE] + prob5[, -1L, drop = FALSE]  # n × 4
  best_pair_start <- max.col(pair_sums, ties.method = "first")              # which pair
  cr_two_row <- as.integer(j5 == best_pair_start | j5 == best_pair_start + 1L)
  cr_two_row[is.na(j5)] <- NA_integer_
  df$cr_two_row <- cr_two_row

  idx_all <- seq_len(n)
  if (!is.null(allowed_cells)) {
    allowed_ids <- allowed_cells %>%
      dplyr::distinct(id)

    keep <- df %>%
      dplyr::transmute(
        id,
        row_id = dplyr::row_number()
      ) %>%
      dplyr::semi_join(allowed_ids, by = "id")

    idx_all <- keep$row_id
  }

  err5_all  <- err5[idx_all, , drop = FALSE]
  obs5_all  <- obs5[idx_all, , drop = FALSE]
  prob5_all <- prob5[idx_all, , drop = FALSE]
  rps_all   <- rps_row[idx_all]

  brier_by <- colMeans(err5_all, na.rm = TRUE)
  names(brier_by) <- paste0("brier_", names(brier_by))

  auc_by <- vapply(seq_along(outcome_levels), function(jj) {
    fast_auc(obs5_all[, jj], prob5_all[, jj])
  }, numeric(1))
  names(auc_by) <- paste0("auc_", outcome_levels)

  auc_all <- fast_auc(as.integer(c(obs5_all)), as.numeric(c(prob5_all)))

  preds_all <- as.numeric(c(prob5_all))
  obs_all   <- as.integer(c(obs5_all))
  pos <- preds_all[obs_all == 1L]
  neg <- preds_all[obs_all == 0L]
  if(min(pos, neg) >= 5){  
    pietra <- as.numeric(suppressWarnings(stats::ks.test(pos, neg)$statistic))
  }else{
    pietra <- NA
  }

  full_metrics <- data.table::data.table(
    id = "ALL", lat = NA_character_, lon = NA_character_,
    brier  = sum(err5_all, na.rm = TRUE) / nrow(err5_all),
    rps    = mean(rps_all, na.rm = TRUE),
    auc    = auc_all,
    cr_one = mean(cr_one_row[idx_all], na.rm = TRUE),
    cr_two = mean(cr_two_row[idx_all], na.rm = TRUE),
    n      = length(idx_all),
    pietra = pietra
  )

  for (nm in names(brier_by)) full_metrics[[nm]] <- brier_by[[nm]]
  for (nm in names(auc_by))   full_metrics[[nm]] <- auc_by[[nm]]

  if (print_outputs) print(full_metrics)

  if (isTRUE(only_all)) {
    return(tibble::as_tibble(full_metrics))
  }

  dt <- data.table::as.data.table(df)

  cell_base <- dt[, .(
    lat   = lat[1L],
    lon   = lon[1L],
    brier  = mean(brier_row, na.rm = TRUE),
    rps    = mean(rps_row, na.rm = TRUE),
    cr_one = mean(cr_one_row, na.rm = TRUE),
    cr_two = mean(cr_two_row, na.rm = TRUE),
    n      = .N
  ), by = .(id)]

  split_idx <- split(seq_len(n), dt$id)

  auc_fun <- function(idx) {
    obs_vec <- as.integer(c(obs5[idx, , drop = FALSE]))
    pred_vec <- as.numeric(c(prob5[idx, , drop = FALSE]))
    fast_auc(obs_vec, pred_vec)
  }

  if (parallel_auc && length(split_idx) > 1L) {
    if (!requireNamespace("future.apply", quietly = TRUE) ||
        !requireNamespace("future", quietly = TRUE)) {
      stop("Please install future and future.apply for parallel_auc=TRUE.")
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    auc_by_cell <- future.apply::future_lapply(split_idx, auc_fun)
  } else {
    auc_by_cell <- lapply(split_idx, auc_fun)
  }

  auc_dt <- data.table::data.table(
    id = names(auc_by_cell),
    auc = as.numeric(auc_by_cell)
  )

  cell_metrics <- merge(cell_base, auc_dt, by = "id", all.x = TRUE)

  out <- data.table::rbindlist(list(cell_metrics, full_metrics), fill = TRUE)
  tibble::as_tibble(out)
}

calibration_from_wide_probs_pooled_quantiles <- function(cv_preds,
                                                         outcome_col = "outcome",
                                                         prob_prefix = "cv_",
                                                         bins = c("earlier","week1","week2","week3","week4","later"),
                                                         nbins = 10) {
  stopifnot(is.data.frame(cv_preds))
  if (!outcome_col %in% names(cv_preds)) {
    stop("Outcome column '", outcome_col, "' not found in cv_preds.")
  }

  prob_cols <- paste0(prob_prefix, bins)
  keep <- prob_cols %in% names(cv_preds)
  bins <- bins[keep]
  prob_cols <- prob_cols[keep]

  if (length(prob_cols) == 0) {
    stop("No probability columns found among: ", paste(prob_cols, collapse = ", "))
  }

  long <- cv_preds %>%
    dplyr::select(dplyr::all_of(outcome_col), dplyr::all_of(prob_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(prob_cols),
      names_to = "forecast_bin",
      values_to = "p"
    ) %>%
    dplyr::mutate(
      forecast_bin = sub(paste0("^", prob_prefix), "", forecast_bin),
      outcome_chr  = as.character(.data[[outcome_col]]),
      y = as.integer(outcome_chr == forecast_bin)
    ) %>%
    dplyr::filter(is.finite(p), !is.na(y))

  cal <- long %>%
    dplyr::mutate(qbin = dplyr::ntile(p, nbins)) %>%
    dplyr::group_by(qbin) %>%
    dplyr::summarise(
      n         = dplyr::n(),
      pred_mean = mean(p, na.rm = TRUE),
      obs_frac  = mean(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(pred_mean)

  cal
}

plot_calibration_pooled_with_hist <- function(cal_df,
                                              preds,
                                              title_text = "Reliability (pooled)",
                                              show_diag = TRUE,
                                              bins = nrow(cal_df),
                                              hist_alpha = 0.25) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2: install.packages('ggplot2')")
  }

  ggplot2::ggplot() +
    ggplot2::geom_histogram(
      data = data.frame(pred = preds),
      ggplot2::aes(x = pred, y = ggplot2::after_stat(count / max(count))),
      bins = bins,
      boundary = 0,
      closed = "left",
      alpha = hist_alpha
    )+
    { if (isTRUE(show_diag)) ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") } +
    ggplot2::geom_line(
      data = cal_df,
      ggplot2::aes(x = pred_mean, y = obs_frac),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      data = cal_df,
      ggplot2::aes(x = pred_mean, y = obs_frac),
      size = 2.2
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = title_text,
      x = "Forecast probability",
      y = "Observed frequency"
    ) +
    ggplot2::theme_minimal()
}

build_reliability_chart_data <- function(cv_preds,
                                         cal_df,
                                         prob_prefix = "cv_",
                                         bins = c("earlier","week1","week2","week3","week4","later"),
                                         hist_bins = nrow(cal_df)) {
  stopifnot(is.data.frame(cv_preds), is.data.frame(cal_df))
  cv_preds <- data.table::as.data.table(cv_preds)

  prob_cols <- paste0(prob_prefix, bins)
  prob_cols <- intersect(prob_cols, names(cv_preds))
  if (length(prob_cols) == 0) stop("No prob columns found for histogram.")

  preds <- as.numeric(unlist(cv_preds[, ..prob_cols], use.names = FALSE))
  preds <- preds[is.finite(preds) & !is.na(preds)]
  preds <- preds[preds >= 0 & preds <= 1]

  edges <- seq(0, 1, length.out = hist_bins + 1)
  b <- cut(preds, breaks = edges, include.lowest = TRUE, right = TRUE)
  h <- as.data.frame(table(b), stringsAsFactors = FALSE)
  h$count <- as.integer(h$Freq)
  h$Freq <- NULL

  rng <- gsub("\\[|\\]|\\(|\\)", "", h$b)
  parts <- strsplit(rng, ",", fixed = TRUE)
  left  <- as.numeric(vapply(parts, `[`, "", 1))
  right <- as.numeric(vapply(parts, `[`, "", 2))
  mid   <- (left + right) / 2

  hist_df <- dplyr::tibble(
    bin_left  = left,
    bin_right = right,
    bin_mid   = mid,
    count     = h$count
  ) %>%
    dplyr::arrange(.data$bin_left)

  denom <- max(hist_df$count, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) denom <- 1

  hist_df <- hist_df %>%
    dplyr::mutate(count_scaled_01 = .data$count / denom)

  cal_out <- cal_df %>%
    dplyr::mutate(
      series = "calibration",
      x = .data$pred_mean,
      y = .data$obs_frac
    ) %>%
    dplyr::select(series, dplyr::everything(), x, y)

  hist_out <- hist_df %>%
    dplyr::mutate(
      series = "hist",
      x = .data$bin_mid,
      y = .data$count_scaled_01
    )

  dplyr::bind_rows(cal_out, hist_out)
}

save_reliability_chart_rds <- function(chart_df, file) {
  stopifnot(is.data.frame(chart_df))
  saveRDS(chart_df, file)
  invisible(file)
}

plot_reliability_from_chart_rds <- function(rds_file,
                                            title_text = "Reliability (pooled)",
                                            show_diag = TRUE,
                                            hist_alpha = 0.25) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2: install.packages('ggplot2')")
  }

  df <- readRDS(rds_file)

  cal_df  <- df[df$series == "calibration", , drop = FALSE]
  hist_df <- df[df$series == "hist", , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = hist_df,
      ggplot2::aes(x = .data$bin_mid, y = .data$y),
      alpha = hist_alpha,
      width = (hist_df$bin_right[1] - hist_df$bin_left[1])
    ) +
    { if (isTRUE(show_diag)) ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") } +
    ggplot2::geom_line(
      data = cal_df,
      ggplot2::aes(x = .data$pred_mean, y = .data$obs_frac),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      data = cal_df,
      ggplot2::aes(x = .data$pred_mean, y = .data$obs_frac),
      size = 2.2
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = title_text,
      x = "Forecast probability",
      y = "Observed frequency"
    ) +
    ggplot2::theme_minimal()
}

compute_yearly_from_preds <- function(pred_df, model_name, allowed_cells = NULL) {
  pred_df %>%
    dplyr::group_by(year) %>%
    dplyr::group_modify(
      ~ compute_cell_metrics_fast(.x,
                                  allowed_cells = allowed_cells,
                                  only_all = TRUE,
                                  print_outputs = FALSE) %>%
        dplyr::select(brier, rps, auc, cr_one, cr_two) %>%
        dplyr::slice(1)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(model = model_name)
}
