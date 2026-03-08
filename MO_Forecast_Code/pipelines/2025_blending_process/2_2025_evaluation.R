# ==============================================================================
# File: 2_2025_evaluation.R
# ==============================================================================
# Purpose
#   Out-of-sample evaluation of 2025 monsoon onset forecasts. Computes Brier,
#   AUC, and RPS for the blended model, conditional and unconditional
#   climatologies (from the SENT forecast file), and Platt-calibrated NGCM
#   forecasts. Scores are computed against three ground-truth variants
#   (mr_mok_year, mr_mok_clim, mr_may1).
#
# Inputs
#   - Monsoon_Data/evaluation_2025/blended_forecasts_2025_sent.csv
#   - Monsoon_Data/evaluation_2025/ngcm_forecasts_2025_sent.csv
#   - Monsoon_Data/evaluation_2025/ngcm_forecasts_2025_sent_no_mok_filter.csv
#   - Monsoon_Data/evaluation_2025/mr_onset_with_mok_year_2025.csv
#   - Monsoon_Data/evaluation_2025/mr_onset_with_mok_clim_2025.csv
#   - Monsoon_Data/evaluation_2025/mr_onset_with_may1_2025.csv
#   - Monsoon_Data/results/2025_model_evaluation/platt_weights_*_calibrated_df_*_2000_2024.rds
#
# Outputs
#   - Monsoon_Data/results/2025_model_evaluation/evaluation/metrics_min_brier_rps_auc_models.csv
#   - Monsoon_Data/results/2025_model_evaluation/evaluation/model_metrics_sent_vs_<on_lbl>.rds
#   - Monsoon_Data/results/2025_model_evaluation/yearly_metrics_2025{_suffix}.rds
#
# Dependencies
#   - R/2025_blending_process/evaluation_2025_utils.R
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(pROC)
  library(purrr)
  library(conflicted)
})
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("na_if", "dplyr")
conflict_prefer("map", "purrr")

source("R/2025_blending_process/evaluation_2025_utils.R")

# ---- PATHS ----
eval_dir <- "Monsoon_Data/evaluation_2025"
out_dir  <- "Monsoon_Data/results/2025_model_evaluation/evaluation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

platt_dir <- "Monsoon_Data/results/2025_model_evaluation"

forecast_files <- list(
  sent = file.path(eval_dir, "blended_forecasts_2025_sent.csv"),
  ngcm = file.path(eval_dir, "ngcm_forecasts_2025_sent.csv"),
  ngcm_no_mok_filter = file.path(eval_dir, "ngcm_forecasts_2025_sent_no_mok_filter.csv")
)

onset_files <- list(
  mr_mok_year = file.path(eval_dir, "mr_onset_with_mok_year_2025.csv"),
  mr_mok_clim = file.path(eval_dir, "mr_onset_with_mok_clim_2025.csv"),
  mr_may1     = file.path(eval_dir, "mr_onset_with_may1_2025.csv")
)

# ==============================================================================
# MAIN
# ==============================================================================

# Iterate over ground-truth variants (SENT is the baseline forecast source)
combos <- tibble(on_lbl = names(onset_files))

metrics_all <- purrr::pmap_dfr(
  combos,
  function(on_lbl) {

    lbl <- paste0("sent_vs_", on_lbl)

    # --- Read SENT baseline forecast file ---
    if (!file.exists(forecast_files[["sent"]])) stop("Forecast file missing: ", forecast_files[["sent"]])

    sent_fc <- readr::read_csv(forecast_files[["sent"]], show_col_types = FALSE) %>%
      mutate(time = parse_any_date(time),
             id = paste0(lat, "_", lon)) %>%
      filter(!is.na(week1))

    # Require that SENT contains both climatology families as separate columns
    clim_cols <- c(paste0("clim_week", 1:4), "clim_later")
    clim_con_cols <- c(paste0("clim_con_week", 1:4), "clim_con_later")

    missing_clim <- setdiff(clim_cols, names(sent_fc))
    missing_clim_con <- setdiff(clim_con_cols, names(sent_fc))

    if (length(missing_clim) > 0) {
      stop("SENT file is missing required climatology columns: ", paste(missing_clim, collapse = ", "))
    }
    if (length(missing_clim_con) > 0) {
      stop("SENT file is missing required conditional climatology columns: ", paste(missing_clim_con, collapse = ", "))
    }

    # --- Read onsets ---
    on <- read_onsets(onset_files[[on_lbl]])

    # Write pre-onset SENT file
    sent_pre_onset <- sent_fc %>%
      inner_join(on, by = "id") %>%
      filter(!is.na(time), !is.na(true_onset)) %>%
      filter(time < true_onset) %>%
      select(-true_onset)
    write_csv(
      sent_pre_onset,
      file.path(out_dir, paste0("blended_forecasts_2025_sent_preonset_", on_lbl, ".csv"))
    )

    # --- Join & outcome (SENT baseline scoring frame) ---
    df <- sent_fc %>%
      left_join(on, by = "id") %>%
      filter(!is.na(true_onset), !is.na(time)) %>%
      mutate(
        diff_days = as.integer(true_onset - time),
        outcome   = ifelse(diff_days >= 0, diff_to_bin(diff_days), NA_character_)
      ) %>%
      filter(!is.na(outcome), diff_days > 0)

    # --- Multiclass Brier + macro AUC for all baseline models ---
    # Uses the same matrix-based computation as 1_blend_evaluation.R
    # (rowSums of squared errors across 5 bins, averaged over samples).
    outcome_chr_base <- df$outcome
    Y_base <- vapply(INTERVAL_BINS_5, function(b) as.integer(outcome_chr_base == b), numeric(length(outcome_chr_base)))
    Y_base <- matrix(Y_base, nrow = length(outcome_chr_base), ncol = length(INTERVAL_BINS_5), byrow = FALSE)

    # Build probability matrices for each model family
    P_forecast <- df %>% select(all_of(INTERVAL_BINS_5)) %>% as.matrix()
    P_clim     <- df %>% select(all_of(clim_cols)) %>% as.matrix()
    P_clim_con <- df %>% select(all_of(clim_con_cols)) %>% as.matrix()

    # Multiclass Brier: mean(rowSums((P - Y)^2))
    brier_forecast <- mean(rowSums((P_forecast - Y_base)^2))
    brier_clim     <- mean(rowSums((P_clim - Y_base)^2))
    brier_clim_con <- mean(rowSums((P_clim_con - Y_base)^2))

    # Pooled AUC: stack all bins into one long binary dataset
    # (matches compute_cell_metrics_fast in blend_evaluation_utils.R)
    pooled_auc <- function(P, Y) {
      y_long <- as.integer(c(Y))
      p_long <- as.numeric(c(P))
      if (length(unique(y_long)) < 2) return(NA_real_)
      as.numeric(pROC::auc(pROC::roc(y_long, p_long, quiet = TRUE)))
    }

    auc_forecast <- pooled_auc(P_forecast, Y_base)
    auc_clim     <- pooled_auc(P_clim, Y_base)
    auc_clim_con <- pooled_auc(P_clim_con, Y_base)

    # --- Skills vs unc_clim_raw ---
    brier_skill_forecast <- skill_score(brier_forecast, brier_clim)
    brier_skill_clim_con <- skill_score(brier_clim_con, brier_clim)

    # --- RPS (6 bins) + skills vs unc_clim_raw ---
    rps_forecast <- rps_frame(make_probs_for_rps(df, "forecast"), df$outcome) %>% mean(na.rm = TRUE)
    rps_clim     <- rps_frame(make_probs_for_rps(df, "clim"),     df$outcome) %>% mean(na.rm = TRUE)
    rps_clim_con <- rps_frame(make_probs_for_rps(df, "clim_con"), df$outcome) %>% mean(na.rm = TRUE)

    rps_skill_forecast <- skill_score(rps_forecast, rps_clim)
    rps_skill_clim_con <- skill_score(rps_clim_con, rps_clim)

    # --- Baseline 3-row output (SENT only) ---
    out <- tibble(
      dataset = lbl,
      model   = c("blended_model", "clim_raw", "unc_clim_raw"),
      brier_skill = c(brier_skill_forecast, brier_skill_clim_con, 0),
      rps_skill   = c(rps_skill_forecast,   rps_skill_clim_con,   0),
      brier       = c(brier_forecast, brier_clim_con, brier_clim),
      rps         = c(rps_forecast, rps_clim_con, rps_clim),
      auc         = c(auc_forecast, auc_clim_con, auc_clim)
    )

    # --- NGCM Platt-calibrated "calibrated" row ---
    ngcm_path <- ngcm_path_from_onset(on_lbl, forecast_files)
    if (!file.exists(ngcm_path)) {
      message("NGCM forecast file missing for ", lbl, ": ", ngcm_path, " (skipping ngcm_calibrated)")
    } else {

      ngcm_fc <- readr::read_csv(ngcm_path, show_col_types = FALSE) %>%
        mutate(time = parse_any_date(time),
               id = paste0(lat, "_", lon)) %>%
        filter(!is.na(week1))

      # Ensure NGCM frame has clim columns (join from SENT if missing)
      missing_clim_ngcm <- setdiff(clim_cols, names(ngcm_fc))
      if (length(missing_clim_ngcm) > 0) {
        ngcm_fc <- ngcm_fc %>%
          left_join(
            sent_fc %>% select(any_of(c("id", "time", clim_cols))),
            by = c("id", "time")
          )
      }

      df_ngcm <- ngcm_fc %>%
        left_join(on, by = "id") %>%
        filter(!is.na(true_onset), !is.na(time)) %>%
        mutate(
          diff_days = as.integer(true_onset - time),
          outcome   = ifelse(diff_days >= 0, diff_to_bin(diff_days), NA_character_)
        ) %>%
        filter(!is.na(outcome), diff_days > 0)

      platt_cfg <- platt_config_from_onset(on_lbl)
      platt_obj <- read_platt_weights(platt_dir, platt_cfg$label, platt_cfg$cutoff_tag)

      if (!platt_obj$ok) {
        message("Platt weights not found for ", lbl, ": ", platt_obj$path, " (skipping ngcm_calibrated)")
      } else if (nrow(df_ngcm) == 0) {
        message("No NGCM rows to score for ", lbl, " (skipping ngcm_calibrated)")
      } else {

        # Raw NGCM probs -> calibrated -> renormalized to sum 1 across 5 bins
        P_raw <- df_ngcm %>% select(all_of(INTERVAL_BINS_5)) %>% as.matrix()
        P_cal <- apply_platt_5(P_raw, platt_obj$df)

        # Brier + pooled AUC on 5 bins
        outcome_chr <- df_ngcm$outcome
        Y <- vapply(INTERVAL_BINS_5, function(b) as.integer(outcome_chr == b), numeric(length(outcome_chr)))
        Y <- matrix(Y, nrow = length(outcome_chr), ncol = length(INTERVAL_BINS_5), byrow = FALSE)

        brier_ngcm <- mean(rowSums((P_cal - Y)^2))
        auc_ngcm   <- pooled_auc(P_cal, Y)

        # 6-bin RPS
        probs_ngcm_rps <- tibble(
          week1 = P_cal[, 1],
          week2 = P_cal[, 2],
          week3 = P_cal[, 3],
          week4 = P_cal[, 4],
          later = P_cal[, 5]
        ) %>%
          mutate(earlier = pmax(0, 1 - (week1 + week2 + week3 + week4 + later))) %>%
          select(all_of(RPS_BINS))

        rps_ngcm <- rps_frame(probs_ngcm_rps, df_ngcm$outcome) %>% mean(na.rm = TRUE)

        # Skills relative to SENT unc_clim_raw (same reference as baseline)
        brier_ref <- out %>% filter(model == "unc_clim_raw") %>% pull(brier) %>% as.double()
        rps_ref   <- out %>% filter(model == "unc_clim_raw") %>% pull(rps)   %>% as.double()

        ngcm_row <- tibble(
          dataset = lbl,
          model   = "ngcm_calibrated",
          brier_skill = skill_score(brier_ngcm, brier_ref),
          rps_skill   = skill_score(rps_ngcm,   rps_ref),
          brier       = brier_ngcm,
          rps         = rps_ngcm,
          auc         = auc_ngcm
        )

        out <- bind_rows(out, ngcm_row)
      }
    }

    # Save per-dataset RDS
    saveRDS(out, file.path(out_dir, paste0("model_metrics_", lbl, ".rds")))

    # Save yearly metrics for 3_produce_figures.R (matches yearly_metrics_global format)
    # Remap NGCM model name to match 1_blend_evaluation.R convention for bind_rows
    ngcm_yearly_nm <- ngcm_yearly_name_from_onset(on_lbl)
    yearly_suffix  <- yearly_suffix_from_onset(on_lbl)
    yearly_out <- out %>%
      mutate(model = dplyr::if_else(model == "ngcm_calibrated", ngcm_yearly_nm, model)) %>%
      transmute(year = 2025L, model, brier, brier_skill, rps, rpss = rps_skill, auc)
    saveRDS(yearly_out, file.path(platt_dir, paste0("yearly_metrics_2025", yearly_suffix, ".rds")))

    out
  }
)

# ---- SAVE FINAL OUTPUTS ----
write_csv(metrics_all, file.path(out_dir, "metrics_min_brier_rps_auc_models.csv"))

# Convenience split
write_csv(metrics_all %>% filter(str_detect(dataset, "^sent_vs_")), file.path(out_dir, "metrics_min_SENT_plus_NGCMcalibrated.csv"))

