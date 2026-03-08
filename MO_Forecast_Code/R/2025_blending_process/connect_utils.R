# ==============================================================================
# File: connect_utils.R
# ==============================================================================
# Purpose
#   Helper functions for 0_connect_prepare_data_to_2025_pipeline.R. Provides
#   day-to-week aggregation, logit-winsorization, rolling-window rain summaries,
#   and the main make_cv_rds_from_daylevel() converter.
#
# Function index
#   winsor_weekp(p, lo, hi)
#     Winsorize probability to [lo, hi]. Default [0.0001, 0.9999].
#
#   logit_winsor(p, lo, hi)
#     Winsorize then apply qlogis (logit transform).
#
#   sum_week_probs_from_dayprefix(df, day_prefix, out_prefix, day_max, ...)
#     Sum daily columns <day_prefix>1..<day_prefix>N into weekly bins.
#     df: data.table; day_prefix: string. Returns: tibble with week columns.
#
#   sum_week_probs(df, prefix, day_max, days_per_week, n_weeks)
#     Sum <prefix>_p_onset_day_1..<day_max> into <prefix>_p_onset_week1..N.
#
#   sum_week_probs_with_day0(df, prefix, ...)
#     Like sum_week_probs but also extracts day_0 column for "earlier" bin.
#     Returns list(day0, week).
#
#   make_clim_logits_from_prefix(raw, input_prefix, output_tag, ...)
#     Aggregate climatology day probs to weeks and apply logit_winsor.
#
#   roll_sums_mat(mat, k)
#     Rolling k-day row sums from a wide matrix. Returns [nrow x (ncol-k+1)].
#
#   week_max_over_starts(roll_mat, week_start_days)
#     Per-row max of rolling sums at specified start days.
#
#   week_min_over_starts(roll_mat, week_start_days)
#     Per-row min of rolling sums at specified start days.
#
#   make_cv_rds_from_daylevel(spec)
#     Main converter: reads daily combined RDS, builds weekly bins, onset
#     outcomes, climatology logits, rain-based predictors, and writes wide
#     RDS for the 2025 blending pipeline.
#     spec: parsed YAML list with mode, input_rds, output_rds, forecast_models, etc.
#
# Dependencies
#   - dplyr, stringr, lubridate, tidyr, purrr, tibble, data.table
#   - R/_shared/misc.R (for %||% operator; must be sourced by the calling script)
# ==============================================================================
library(data.table)
winsor_weekp <- function(p, lo = .0001, hi = 0.9999) {
  pmin(pmax(p, lo), hi)
}

logit_winsor <- function(p, lo = .0001, hi = 0.9999) {
  qlogis(winsor_weekp(p, lo, hi))
}

# Sum probs by week from columns like <day_prefix>1 .. <day_prefix>N
# where day_prefix already includes the "day_" part, e.g. "ngcm_p_onset_mok_day_"
sum_week_probs_from_dayprefix <- function(df, day_prefix, out_prefix,
                                          day_max = 28, days_per_week = 7, n_weeks = 4) {
  cols <- paste0(day_prefix, 1:day_max)
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop("Missing required columns for ", out_prefix, ": ", paste(miss, collapse = ", "))
  }
  mat <- as.matrix(df[, ..cols])

  out <- vector("list", n_weeks)
  for (w in seq_len(n_weeks)) {
    lo <- (w - 1) * days_per_week + 1
    hi <- w * days_per_week
    out[[w]] <- rowSums(mat[, lo:hi, drop = FALSE], na.rm = FALSE)
  }
  names(out) <- paste0(out_prefix, "_week", 1:n_weeks)
  dplyr::as_tibble(out)
}

make_clim_logits_from_prefix <- function(raw, input_prefix, output_tag,
                                         day_max, days_per_week, n_weeks) {
  wk <- sum_week_probs(
    raw,
    prefix = input_prefix,
    day_max = day_max,
    days_per_week = days_per_week,
    n_weeks = n_weeks
  )

  wk %>%
    dplyr::transmute(
      "{paste0('prob_clim_mr_', output_tag, '_week1')}" :=
        logit_winsor(.data[[paste0(input_prefix, "_p_onset_week1")]]),
      "{paste0('prob_clim_mr_', output_tag, '_week2')}" :=
        logit_winsor(.data[[paste0(input_prefix, "_p_onset_week2")]]),
      "{paste0('prob_clim_mr_', output_tag, '_week3')}" :=
        logit_winsor(.data[[paste0(input_prefix, "_p_onset_week3")]]),
      "{paste0('prob_clim_mr_', output_tag, '_week4')}" :=
        logit_winsor(.data[[paste0(input_prefix, "_p_onset_week4")]])
    )
}

# Sum probs by week from p_onset_day_1..p_onset_day_N columns
sum_week_probs <- function(df, prefix, day_max = 28, days_per_week = 7, n_weeks = 4) {
  cols <- paste0(prefix, "_p_onset_day_", 1:day_max)
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop("Missing required columns for ", prefix, ": ", paste(miss, collapse = ", "))
  }
  mat <- as.matrix(df[, ..cols])

  out <- vector("list", n_weeks)
  for (w in seq_len(n_weeks)) {
    lo <- (w - 1) * days_per_week + 1
    hi <- w * days_per_week
    out[[w]] <- rowSums(mat[, lo:hi, drop = FALSE], na.rm = FALSE)
  }
  names(out) <- paste0(prefix, "_p_onset_week", 1:n_weeks)
  dplyr::as_tibble(out)
}

# Sum probs by week for an "unc clim" series that also has day_0
sum_week_probs_with_day0 <- function(df, prefix, day_max = 28, days_per_week = 7, n_weeks = 4) {
  col0 <- paste0(prefix, "_p_onset_day_0")
  if (!(col0 %in% names(df))) stop("Missing required column: ", col0)

  week_tbl <- sum_week_probs(df, prefix, day_max = day_max, days_per_week = days_per_week, n_weeks = n_weeks)
  list(
    day0 = df[[col0]],
    week = week_tbl
  )
}

# Rolling k-day sums from day_1..day_N columns, returned as a matrix [nrow x (N-k+1)]
roll_sums_mat <- function(mat, k) {
  n <- ncol(mat)
  if (n < k) {
    return(matrix(NA_real_, nrow = nrow(mat), ncol = 0))
  }
  starts <- 1:(n - k + 1)
  out <- vapply(starts, function(s) rowSums(mat[, s:(s + k - 1), drop = FALSE], na.rm = FALSE),
                numeric(nrow(mat)))
  t(t(out))  # ensure [nrow x nstarts]
}

week_max_over_starts <- function(roll_mat, week_start_days) {
  if (ncol(roll_mat) == 0) return(rep(NA_real_, nrow(roll_mat)))
  ok <- week_start_days[week_start_days >= 1 & week_start_days <= ncol(roll_mat)]
  if (length(ok) == 0) return(rep(NA_real_, nrow(roll_mat)))
  apply(roll_mat[, ok, drop = FALSE], 1, max, na.rm = TRUE)
}

week_min_over_starts <- function(roll_mat, week_start_days) {
  if (ncol(roll_mat) == 0) return(rep(NA_real_, nrow(roll_mat)))
  ok <- week_start_days[week_start_days >= 1 & week_start_days <= ncol(roll_mat)]
  if (length(ok) == 0) return(rep(NA_real_, nrow(roll_mat)))
  apply(roll_mat[, ok, drop = FALSE], 1, min, na.rm = TRUE)
}

# Main conversion: daily combined data -> weekly RDS
make_cv_rds_from_daylevel <- function(spec) {
  mode <- spec$mode
  input_rds <- spec$input_rds
  output_rds <- spec$output_rds
  day_max <- spec$day_max %||% 28
  days_per_week <- spec$days_per_week %||% 7
  n_weeks <- spec$n_weeks %||% 4

  raw <- readRDS(input_rds)

  # Parse lat/lon from id like "10_78"
  if (!("id" %in% names(raw))) stop("Expected an 'id' column like 'lat_lon' (e.g., 10_78).")
  id_parts <- stringr::str_split(raw$id, "_", simplify = TRUE)
  if (ncol(id_parts) < 2) stop("Couldn't parse lat/lon from id. Expected format 'lat_lon'.")

  raw <- raw %>%
    dplyr::mutate(
      time = as.Date(time),
      lat = as.numeric(id_parts[, 1]),
      lon = as.numeric(id_parts[, 2]),
      year = lubridate::year(time)
    )

  # Onset threshold (used for diff_* computations and also needed downstream)
  # Use first available forecast model's onset_thresh column
  first_model <- spec$forecast_models[[1]]$name
  thresh_col <- paste0(first_model, "_onset_thresh")
  if (!(thresh_col %in% names(raw))) stop("Missing ", thresh_col, " (onset threshold).")
  raw <- raw %>% dplyr::mutate(onset_threshold = .data[[thresh_col]])

  # Outcome: bin true onset date relative to forecast init date into week1..week4/later
  if (!("true_onset_date" %in% names(raw))) stop("Missing true_onset_date.")
  raw <- raw %>%
    dplyr::mutate(
      true_onset_date = as.Date(true_onset_date),
      lead_day = as.integer(true_onset_date - time),
      outcome = dplyr::case_when(
        lead_day <= 0 ~ NA,
        lead_day <= 7 ~ "week1",
        lead_day <= 14 ~ "week2",
        lead_day <= 21 ~ "week3",
        lead_day <= 28 ~ "week4",
        !is.na(lead_day) ~ "later",
        TRUE ~ NA
      )
    )

  # Climatology: base prefix week probabilities
  base_prefix <- spec$climatology$base_prefix %||% "clim"
  unc_prefix <- spec$climatology$unconditional_prefix %||% "clim_unc"

  clim_week_probs <- sum_week_probs(raw, base_prefix, day_max = day_max, days_per_week = days_per_week, n_weeks = n_weeks)

  # Climatology window variants (logit-transformed)
  window_tags <- unlist(spec$climatology$window_tags %||% list())
  clim_variant_logits <- if (length(window_tags) > 0) {
    clim_variant_logits_list <- purrr::set_names(window_tags) %>%
      purrr::map(function(tag) {
        pref <- paste0(base_prefix, "_", tag)
        make_clim_logits_from_prefix(
          raw,
          input_prefix = pref,
          output_tag   = tag,
          day_max = day_max,
          days_per_week = days_per_week,
          n_weeks = n_weeks
        )
      })

    dplyr::bind_cols(clim_variant_logits_list)
  } else {
    tibble::tibble(.row = seq_len(nrow(raw))) %>% dplyr::select()
  }

  # Forecast model week probabilities (loop over spec$forecast_models)
  model_week_cols_list <- list()
  for (fm in spec$forecast_models) {
    model_name <- fm$name

    # Base probabilities
    model_week_cols_list[[model_name]] <- sum_week_probs(
      raw, model_name, day_max = day_max, days_per_week = days_per_week, n_weeks = n_weeks
    )

    # Variant probabilities (e.g., clim_mok_date)
    for (variant in fm$variants) {
      variant_key <- paste0(model_name, "_", variant)
      model_week_cols_list[[variant_key]] <- sum_week_probs_from_dayprefix(
        raw,
        day_prefix = paste0(model_name, "_p_onset_", variant, "_day_"),
        out_prefix = paste0(model_name, "_p_onset_", variant),
        day_max = day_max,
        days_per_week = days_per_week,
        n_weeks = n_weeks
      )
    }
  }
  model_week_cols <- dplyr::bind_cols(model_week_cols_list)

  # Unconditional climatology: has day_0 which maps to "earlier"
  unc <- sum_week_probs_with_day0(raw, unc_prefix, day_max = day_max, days_per_week = days_per_week, n_weeks = n_weeks)
  unc_day0 <- unc$day0
  unc_week_probs <- unc$week

  # Logit after winsorization to [0.0001, 0.9999]
  clim_logits <- clim_week_probs %>%
    dplyr::transmute(
      prob_clim_mr_week1 = logit_winsor(.data[[paste0(base_prefix, "_p_onset_week1")]]),
      prob_clim_mr_week2 = logit_winsor(.data[[paste0(base_prefix, "_p_onset_week2")]]),
      prob_clim_mr_week3 = logit_winsor(.data[[paste0(base_prefix, "_p_onset_week3")]]),
      prob_clim_mr_week4 = logit_winsor(.data[[paste0(base_prefix, "_p_onset_week4")]]),

      prob_clim_mr_unc_earlier = logit_winsor(unc_day0),
      prob_clim_mr_unc_week1   = logit_winsor(unc_week_probs[[paste0(unc_prefix, "_p_onset_week1")]]),
      prob_clim_mr_unc_week2   = logit_winsor(unc_week_probs[[paste0(unc_prefix, "_p_onset_week2")]]),
      prob_clim_mr_unc_week3   = logit_winsor(unc_week_probs[[paste0(unc_prefix, "_p_onset_week3")]]),
      prob_clim_mr_unc_week4   = logit_winsor(unc_week_probs[[paste0(unc_prefix, "_p_onset_week4")]])
    )

  # Rain-based predictors (loop over spec$forecast_models)
  rain_predictors_list <- list()
  week_start_days <- lapply(seq_len(n_weeks), function(w) ((w - 1) * days_per_week + 1):(w * days_per_week))

  for (fm in spec$forecast_models) {
    model_name <- fm$name
    rain_preds <- fm$rain_predictors %||% character(0)
    if (length(rain_preds) == 0) next

    # Build rain matrix
    need_rain <- paste0(model_name, "_rain_mean_day_", 1:(day_max + 10))
    if (!all(need_rain %in% names(raw))) {
      stop("Missing ", model_name, " rain columns: ", paste(setdiff(need_rain, names(raw)), collapse = ", "))
    }
    rain_mat <- as.matrix(raw[, ..need_rain])

    # Pre-compute rolling sums needed by this model's predictor set
    needs_5 <- any(rain_preds %in% c("diff_5day", "max_5day"))
    needs_10 <- "min_10day" %in% rain_preds
    if (needs_5) {
      roll5 <- roll_sums_mat(rain_mat, 5)
      max5_by_week <- lapply(week_start_days, function(sd) week_max_over_starts(roll5, sd))
    }
    if (needs_10) {
      roll10 <- roll_sums_mat(rain_mat, 10)
      min10_by_week <- lapply(week_start_days, function(sd) week_min_over_starts(roll10, sd))
    }

    for (pred_type in rain_preds) {
      if (pred_type == "diff_5day") {
        for (w in seq_len(n_weeks)) {
          col_name <- paste0("diff_", model_name, "_week", w)
          rain_predictors_list[[col_name]] <- max5_by_week[[w]] - raw$onset_threshold
        }
      } else if (pred_type == "min_10day") {
        for (w in seq_len(n_weeks)) {
          col_name <- paste0("min_", model_name, "_10day_week", w)
          rain_predictors_list[[col_name]] <- min10_by_week[[w]]
        }
      } else if (pred_type == "max_5day") {
        for (w in seq_len(n_weeks)) {
          col_name <- paste0("max_", model_name, "_5day_week", w)
          rain_predictors_list[[col_name]] <- max5_by_week[[w]]
        }
      } else {
        stop("Unknown rain predictor type: ", pred_type)
      }
    }
  }
  rain_predictors <- tibble::as_tibble(rain_predictors_list)

  wide_df <- raw %>%
    dplyr::transmute(
      id, time, year, lat, lon, onset_threshold, outcome
    ) %>%
    dplyr::bind_cols(
      clim_logits,
      clim_variant_logits,
      model_week_cols,
      rain_predictors
    )

  wide_df$outcome <- as.character(wide_df$outcome)
  out_dir <- dirname(output_rds)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  saveRDS(wide_df, output_rds)
  invisible(wide_df)
}
