# ==============================================================================
# File: evaluation_2025_utils.R
# ==============================================================================
# Purpose
#   Helper functions for 2_2025_evaluation.R. Provides date parsing, onset
#   reading, Brier/AUC/RPS scoring, Platt calibration lookup and application,
#   and onset-label-to-path/suffix mappings.
#
# Dependencies
#   - dplyr, readr, lubridate, pROC, tibble
# ==============================================================================

INTERVAL_BINS_5 <- c("week1", "week2", "week3", "week4", "later")
RPS_BINS <- c("earlier", "week1", "week2", "week3", "week4", "later")

skill_score <- function(score_model, score_ref) 1 - (score_model / score_ref)

# Map onset label to yearly output file suffix
yearly_suffix_from_onset <- function(on_lbl) {
  switch(on_lbl,
         mr_mok_year = "",
         mr_mok_clim = "_clim_mok_date",
         mr_may1     = "_no_mok_filter",
         stop("Unknown onset label: ", on_lbl)
  )
}

parse_any_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  as.Date(lubridate::parse_date_time(
    x,
    orders = c("ymd", "ymd HMS", "mdy", "mdy HMS", "dmy", "dmy HMS")
  ))
}

diff_to_bin <- function(diff_days) {
  dplyr::case_when(
    diff_days < -28 ~ NA_character_,
    diff_days <= -21 ~ "weekm3",
    diff_days <= -14 ~ "weekm2",
    diff_days <=  -7 ~ "weekm1",
    diff_days <=   0 ~ "weekm0",
    diff_days <=   7 ~ "week1",
    diff_days <=  14 ~ "week2",
    diff_days <=  21 ~ "week3",
    diff_days <=  28 ~ "week4",
    TRUE             ~ "later"
  )
}

read_onsets <- function(path) {
  if (!file.exists(path)) stop("Onset file missing: ", path)
  readr::read_csv(path, show_col_types = FALSE) %>%
    dplyr::mutate(
      true_onset = dplyr::na_if(as.character(true_onset), "NaT"),
      true_onset = as.Date(lubridate::parse_date_time(
        true_onset,
        orders = c("dmy", "dmy HMS", "ymd", "ymd HMS", "mdy", "mdy HMS")
      ))
    ) %>%
    dplyr::mutate(id = paste0(lat, "_", lon)) %>%
    dplyr::filter(!is.na(true_onset)) %>%
    dplyr::distinct(id, true_onset)
}

# Choose which NGCM forecast file to use based on ground truth label
ngcm_path_from_onset <- function(on_lbl, forecast_files) {
  if (identical(on_lbl, "mr_may1")) forecast_files[["ngcm_no_mok_filter"]] else forecast_files[["ngcm"]]
}

# Brier + AUC for long-form (bin vs not-bin)
brier_auc <- function(df, prob_col = "prob") {
  d <- df %>% dplyr::filter(!is.na(.data[[prob_col]]), !is.na(outcome), !is.na(bin))
  y <- (d$bin == d$outcome)
  bs <- mean((d[[prob_col]] - y)^2)
  
  auc_val <- NA_real_
  if (length(unique(y)) == 2) {
    roc_obj <- pROC::roc(response = y, predictor = d[[prob_col]], quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
  }
  
  tibble::tibble(brier = bs, auc = auc_val)
}

# RPS for a single row (6 bins), with internal renormalization
rps_row <- function(prob_vec, outcome_bin) {
  if (is.na(outcome_bin)) return(NA_real_)
  if (any(is.na(prob_vec))) return(NA_real_)
  
  s <- sum(prob_vec)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  prob_vec <- prob_vec / s
  
  k <- match(outcome_bin, RPS_BINS)
  if (is.na(k)) return(NA_real_)
  
  ok <- as.numeric(seq_along(RPS_BINS) >= k)
  ck <- cumsum(prob_vec)
  sum((ck - ok)^2)
}

rps_frame <- function(df_probs, outcome_vec) {
  as.numeric(mapply(
    function(i, ob) rps_row(as.numeric(df_probs[i, RPS_BINS]), ob),
    seq_len(nrow(df_probs)), outcome_vec
  ))
}

# Build 6-bin probs for RPS from a df that contains forecast/clim/clim_con cols
make_probs_for_rps <- function(df, family = c("forecast", "clim", "clim_con")) {
  family <- match.arg(family)
  
  df %>%
    dplyr::transmute(
      week1 = dplyr::case_when(
        family == "forecast" ~ week1,
        family == "clim"     ~ clim_week1,
        family == "clim_con" ~ clim_con_week1
      ),
      week2 = dplyr::case_when(
        family == "forecast" ~ week2,
        family == "clim"     ~ clim_week2,
        family == "clim_con" ~ clim_con_week2
      ),
      week3 = dplyr::case_when(
        family == "forecast" ~ week3,
        family == "clim"     ~ clim_week3,
        family == "clim_con" ~ clim_con_week3
      ),
      week4 = dplyr::case_when(
        family == "forecast" ~ week4,
        family == "clim"     ~ clim_week4,
        family == "clim_con" ~ clim_con_week4
      ),
      later = dplyr::case_when(
        family == "forecast" ~ later,
        family == "clim"     ~ clim_later,
        family == "clim_con" ~ clim_con_later
      )
    ) %>%
    dplyr::mutate(
      earlier = pmax(0, 1 - (week1 + week2 + week3 + week4 + later))
    ) %>%
    dplyr::select(dplyr::all_of(RPS_BINS))
}

# ==============================================================================
# Platt calibration helpers
# ==============================================================================

# Platt weight configuration per ground-truth variant.
# Returns list(label, cutoff_tag) matching the file exported by 1_blend_evaluation.R
# for the corresponding cutoff mode and NGCM variant.
# label = base_label used in Platt weight filename (forecast_variant convention)
# cutoff_tag = cutoff mode suffix for the output_tag
platt_config_from_onset <- function(on_lbl) {
  switch(on_lbl,
         mr_mok_year = list(label = "ngcm_clim_mok_date",  cutoff_tag = ""),
         mr_mok_clim = list(label = "ngcm_clim_mok_date",  cutoff_tag = "_clim_mok_date"),
         mr_may1     = list(label = "ngcm",                 cutoff_tag = "_no_mok_filter"),
         stop("Unknown onset label: ", on_lbl)
  )
}

# Model name for NGCM in yearly metrics, matching the convention in 1_blend_evaluation.R
# (forecast_calibrated_variant). Used to align 2025 yearly metrics with 2000-2024 series.
ngcm_yearly_name_from_onset <- function(on_lbl) {
  switch(on_lbl,
         mr_mok_year = "ngcm_calibrated_clim_mok_date",
         mr_mok_clim = "ngcm_calibrated_clim_mok_date",
         mr_may1     = "ngcm_calibrated",
         stop("Unknown onset label: ", on_lbl)
  )
}

# Apply saved Platt weights to a 5-bin probability matrix and renormalize to sum 1.
# weights_df must have columns: bin, intercept, slope
apply_platt_5 <- function(prob_mat, weights_df) {
  bins <- INTERVAL_BINS_5
  cal <- prob_mat
  for (i in seq_along(bins)) {
    b <- bins[i]
    row <- weights_df[weights_df$bin == b, , drop = FALSE]
    if (nrow(row) == 0) next
    p_raw <- pmin(pmax(prob_mat[, i], 1e-6), 1 - 1e-6)
    logit_p <- qlogis(p_raw)
    cal[, i] <- plogis(row$intercept + row$slope * logit_p)
  }
  rs <- rowSums(cal)
  rs[rs == 0] <- 1
  cal / rs
}

# Read Platt weights from file exported by 1_blend_evaluation.R.
# platt_label and cutoff_tag together determine the filename.
read_platt_weights <- function(platt_dir, platt_label, cutoff_tag = "") {
  output_tag_cv <- paste0(cutoff_tag, "_2000_2024")
  
  platt_path <- file.path(
    platt_dir,
    paste0("platt_weights_", platt_label, "_calibrated_df", output_tag_cv, ".rds")
  )
  
  if (!file.exists(platt_path)) {
    return(list(ok = FALSE, path = platt_path, df = NULL))
  }
  
  list(ok = TRUE, path = platt_path, df = readRDS(platt_path))
}