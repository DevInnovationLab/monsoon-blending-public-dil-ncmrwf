# ==============================================================================
# File: climatology_utils.R
#
# Purpose:
#   Utilities for fitting a per-grid-cell climatology model of monsoon onset day
#   from historical (ground-truth) onset dates, and producing issue-date
#   probability forecasts over lead days. Can fit multiple climatology models in one script
#
# High-level workflow:
#   1) Locate inputs/outputs from a YAML spec and a particular choice of climatology:
#        - get_paths_clim()
#        - get_climatology_options_from_run()
#   2) Read and filter ground-truth onset data for KDE training:
#        - read_gt_onset()
#        - filter_gt_training()
#   3) Build the set of issue dates (times) to forecast on:
#        - season_dates_for_year()
#        - build_issue_grid()
#   4) Determine the forecast horizon per issue date (optional piecewise horizons):
#        - resolve_forecast_window_by_time()
#   5) Fit KDEs per cell and compute issue-date forecasts:
#        - fit_kde()
#        - fit_kdes_by_cell()
#        - compute_d0()
#        - predict_from_kde()
#        - compute_forecasts_for_cell()
#        - compute_all_forecasts()
#   conditional (run option)
#     If TRUE (default), forecast probabilities are conditional on onset not having
#     occurred by the issue date. If FALSE, probabilities
#     are unconditional day-mass values from the KDE.
#
# Function reference:
#   get_paths_clim(spec)
#     Derives the input ground-truth onset file path and the output directory/stem
#     for climatology products, using spec$paths overrides when provided and
#     otherwise defaulting to Stage-2 outputs under spec$output.
#
#   get_climatology_options_from_run(co)
#     Extracts climatology configuration for a single run entry from
#     spec$climatologies, returning a standardized options list with fields:
#       train_year_min, conditional, train_year_max, test_year_min, test_year_max,
#       season_start_md, issue_end_md, onset_col, forecast_window, horizons.
#     Mirrors get_climatology_options(), but operates on one run instead of the
#     full spec to support multiple climatology runs per script execution.
#
#
#   read_gt_onset_from_tbl(gt_tbl, onset_col = "mr_onset_day")
#     Reads and standardizes ground-truth onset data from an already-loaded wide
#     ground-truth table (typically *_wide.csv read once by the caller), returning
#     a tibble with columns:
#       id, year, onset_day
#     Errors if id, year, or onset_col are missing. 
#
#   filter_gt_training(gt, y_min, y_max)
#     Filters the ground-truth table to the training year range [y_min, y_max]
#     and drops rows with missing onset_day.
#
#   season_dates_for_year(year, start_md, end_md)
#     For a given year, returns a Date sequence from YYYY-start_md through
#     YYYY-end_md (inclusive), where start_md/end_md are "MM-DD". Validates that
#     the dates parse and that end >= start.
#
#   build_issue_grid(test_year_min, test_year_max, season_start_md, issue_end_md)
#     Constructs the full grid of issue dates ("time") across test/output years,
#     returning a tibble with columns:
#       year, time
#     containing all dates from season_start_md through issue_end_md for each year.
#
#   resolve_forecast_window_by_time(time, forecast_window, horizons)
#     Returns the forecast horizon H(time) as an integer. If horizons is NULL,
#     returns forecast_window (constant horizon). If horizons is provided as a
#     list of {start_md, end_md, forecast_window}, selects the rule whose within-
#     year date range contains time and returns its forecast_window; returns NA if
#     no rule applies.
#
#   fit_kde(x)
#     Fits a 1D KDE to onset_day samples x using stats::density(bw="SJ") and
#     returns the density object. Returns NULL if there are fewer than 10 samples
#     or density fitting fails.
#
#   compute_d0(time, season_start_md)
#     Converts an issue date into an integer offset d0 = (time - season_start),
#     where season_start is YYYY-season_start_md. d0 = 0 on season start.
#
#   predict_from_kde(dens, d0, forecast_window)
#     Given a fitted KDE density object dens (for a single cell), produces a
#     length-H vector of issue-date probabilities over lead days k=1..H:
#
#       p_k = P(onset on day d0+k | onset > d0),
#
#     where probabilities are approximated using differences of the KDE CDF and
#     conditioning on survival past d0. Returns NA vector if dens is NULL or if
#     the conditional denominator is numerically degenerate.
#
#   fit_kdes_by_cell(gt_train)
#     For each (lat, lon) in gt_train, fits a KDE of onset_day and returns a
#     named list of density objects keyed by "lat lon".
#
#   compute_forecasts_for_cell(lat, lon, issue_grid, kdes,
#                              season_start_md, forecast_window, horizons)
#     For a single cell, computes predicted probabilities for every issue date
#     in issue_grid and every lead day 1..H(time), returning a long tibble with:
#       time, day, predicted_prob, model, lat, lon
#
#   compute_all_forecasts(gt_train, issue_grid, season_start_md, forecast_window, horizons)
#     Fits KDEs for all cells from gt_train and computes forecasts for all cells
#     over issue_grid, returning a list:
#       forecasts: long tibble of all cell/date/lead predictions
#       kdes: named list of KDE objects by cell key
# ==============================================================================

get_paths_clim <- function(spec) {
  gt_path <- file.path(spec$output$out_dir, paste0(spec$id, "_wide.rds"))
  
  out_dir  <- spec$paths$climatology_out_dir %||%
    file.path(dirname(spec$output$out_dir), "Climatology")
  
  out_stem <- spec$paths$climatology_out_stem %||%
    "climatology_issue"
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  list(gt_path = gt_path, out_dir = out_dir, out_stem = out_stem)
}

get_climatology_options_from_run <- function(co) {
  list(
    train_year_min  = as.integer(co$train_year_min),
    train_year_max  = as.integer(co$train_year_max),
    test_year_min   = as.integer(co$test_year_min),
    test_year_max   = as.integer(co$test_year_max),
    season_start_md = as.character(co$season_start_md),
    issue_end_md    = as.character(co$issue_end_md),
    onset_col       = as.character(co$onset_col %||% "mr_onset_day"),
    forecast_window = if (!is.null(co$forecast_window)) as.integer(co$forecast_window) else NULL,
    horizons        = co$horizons %||% NULL,
    conditional     = if (!is.null(co$conditional)) isTRUE(co$conditional) else TRUE,
    cv_by_year      = if (!is.null(co$cv_by_year)) isTRUE(co$cv_by_year) else TRUE
  )
}



# ------------------------------------------------------------------------------
# Ground-truth IO (ID-based)
# ------------------------------------------------------------------------------



read_gt_onset_from_tbl <- function(gt_tbl,
                                   onset_col = "mr_onset_day",
                                   na_sentinel = NA) {
  if (!onset_col %in% names(gt_tbl))
    stop("Missing onset column '", onset_col, "' in ground-truth table.", call. = FALSE)
  
  if (!("id" %in% names(gt_tbl)))
    stop("Missing required column 'id' in ground-truth table.", call. = FALSE)
  
  if (!("year" %in% names(gt_tbl)))
    stop("Missing required column 'year' in ground-truth table.", call. = FALSE)
  
  gt_tbl %>%
    transmute(
      id        = as.character(id),
      year      = as.integer(year),
      onset_day = as.integer(.data[[onset_col]])
    ) %>%
    mutate(
      onset_day = if (!is.null(na_sentinel)) tidyr::replace_na(onset_day, as.integer(na_sentinel)) else onset_day
    )
}



filter_gt_training <- function(gt, y_min, y_max) {
  gt %>%
    filter(!is.na(onset_day)) %>%
    filter(year >= y_min, year <= y_max)
}

max_forecast_window <- function(forecast_window, horizons) {
  if (!is.null(horizons)) {
    hs <- purrr::map_int(horizons, ~ as.integer(.x$forecast_window))
    return(max(hs, na.rm = TRUE))
  }
  as.integer(forecast_window)
}

# ----------------------------
# Issue-date grid
# ----------------------------

season_dates_for_year <- function(year, start_md, end_md) {
  start <- as.Date(paste0(year, "-", start_md))
  end   <- as.Date(paste0(year, "-", end_md))
  
  if (is.na(start) || is.na(end)) {
    stop("Bad season dates; expected 'MM-DD'. Got: ", start_md, " / ", end_md, call. = FALSE)
  }
  if (end < start) {
    stop("issue_end_md is before season_start_md for year ", year, ": ", start_md, " to ", end_md, call. = FALSE)
  }
  
  seq.Date(start, end, by = "day")
}

build_issue_grid <- function(test_year_min, test_year_max, season_start_md, issue_end_md) {
  yrs <- seq.int(test_year_min, test_year_max)
  
  purrr::map_dfr(yrs, function(y) {
    tibble(
      year = as.integer(y),
      time = season_dates_for_year(y, season_start_md, issue_end_md)
    )
  }) %>%
    distinct(year, time) %>%
    arrange(year, time)
}
# ----------------------------
# Horizon: H(time)
# ----------------------------

resolve_forecast_window_by_time <- function(time, forecast_window, horizons) {
  if (is.null(horizons)) return(as.integer(forecast_window))
  
  yr <- as.integer(format(time, "%Y"))
  
  rules <- purrr::map_dfr(horizons, function(h) {
    tibble(
      start = as.Date(paste0(yr, "-", as.character(h$start_md))),
      end   = as.Date(paste0(yr, "-", as.character(h$end_md))),
      H     = as.integer(h$forecast_window)
    )
  })
  
  hit <- rules %>% filter(time >= start, time <= end)
  if (nrow(hit) == 0L) return(NA_integer_)
  hit$H[[1]]
}

# ----------------------------
# KDE + forecast math
# ----------------------------

fit_kde <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 10) return(NULL)
  
  dens <- try(stats::density(x, bw = "SJ"), silent = TRUE)
  if (inherits(dens, "try-error")) return(NULL)
  
  dens
}

compute_d0 <- function(time, season_start_md) {
  yr <- as.integer(format(time, "%Y"))
  season_start <- as.Date(paste0(yr, "-", season_start_md))
  as.integer(time - season_start)
}

predict_from_kde <- function(dens, d0, forecast_window,
                             conditional = TRUE,
                             include_day0 = FALSE,
                             min_prob = 0.0000005,
                             eps = 1e-12) {
  
  # helper: enforce per-entry lower bound while matching a target sum (<= 1)
  enforce_floor_with_target_sum <- function(p_raw, target_sum, lb) {
    H <- length(p_raw)
    
    if (!is.finite(target_sum) || target_sum < 0) target_sum <- 0
    if (target_sum > 1) target_sum <- 1
    
    p_raw <- pmax(p_raw, 0)
    p_raw[!is.finite(p_raw)] <- 0
    
    # If we cannot satisfy the floor given the target sum, revert to uniform at target_sum/H
    if (target_sum <= eps) return(rep(0, H))
    if (target_sum < H * lb) return(rep(target_sum / H, H))
    
    p <- pmax(p_raw, lb)
    
    # Reduce only the mass above the floor to hit target_sum exactly
    excess <- sum(p) - target_sum
    if (excess <= eps) return(p)
    
    reducible <- p - lb
    reducible_sum <- sum(reducible)
    
    if (reducible_sum <= eps) {
      return(rep(target_sum / H, H))
    }
    
    p <- p - excess * (reducible / reducible_sum)
    p <- pmax(p, lb)  # numerical safety
    p
  }
  
  if (is.null(dens)) {
    out_len <- forecast_window + if (include_day0 && !conditional) 1L else 0L
    return(rep(NA_real_, out_len))
  }
  
  cdf <- stats::approxfun(dens$x, cumsum(dens$y) / sum(dens$y), rule = 2)
  
  days <- seq_len(forecast_window)
  
  # Unconditional day masses for day_1..day_H (relative to issue date)
  num <- cdf(d0 + days) - cdf(d0 + days - 1)
  num <- pmax(num, 0)
  num[!is.finite(num)] <- 0
  
  if (!conditional) {
    # total mass within horizon (do NOT force to 1)
    target_sum <- sum(num)
    p_adj <- enforce_floor_with_target_sum(num, target_sum, lb = min_prob)
    
    if (include_day0) {
      day0 <- cdf(d0) # P(onset <= d0)
      return(c(day0, p_adj))
    }
    return(p_adj)
  }
  
  # Conditional normalization: P(onset on d0+k | onset > d0)
  base_prob <- cdf(d0)
  denom <- 1 - base_prob
  
  # If conditioning is degenerate, avoid NA: fall back to a "reasonable" distribution
  # Here: keep the in-horizon mass as if denom were 1 (i.e., use num) but smooth it.
  if (!is.finite(denom) || denom <= eps) {
    target_sum <- sum(num) # still <= 1
    return(enforce_floor_with_target_sum(num, target_sum, lb = min_prob))
  }
  
  p_raw <- num / denom
  p_raw <- pmax(p_raw, 0)
  p_raw[!is.finite(p_raw)] <- 0
  
  # Keep the *in-horizon* conditional mass (don’t force to 1)
  target_sum <- sum(p_raw)
  if (!is.finite(target_sum) || target_sum <= eps) {
    # fallback: spread a small amount of mass (or zero if truly none)
    return(rep(0, forecast_window))
  }
  
  enforce_floor_with_target_sum(p_raw, target_sum, lb = min_prob)
}

# ------------------------------------------------------------------------------
# KDEs by id (replaces lat/lon cell keys)
# ------------------------------------------------------------------------------

fit_kdes_by_cell <- function(gt_train) {
  cells <- gt_train %>% distinct(id)
  
  kdes_df <- cells %>%
    rowwise() %>%
    mutate(
      cell_key = as.character(id),
      kde = {
        cur_id <- as.character(id)  # capture the rowwise id into a plain variable
        list(
          fit_kde(
            gt_train %>%
              filter(id == cur_id) %>%
              pull(onset_day)
          )
        )
      }
    ) %>%
    ungroup()
  
  out <- kdes_df$kde
  names(out) <- kdes_df$cell_key
  out
}

# ------------------------------------------------------------------------------
# Forecasts for a single id (replaces lat/lon signature)
# ------------------------------------------------------------------------------

compute_forecasts_for_cell <- function(id, issue_grid, kdes,
                                       season_start_md, forecast_window, horizons,
                                       conditional = TRUE,
                                       cv_by_year = FALSE,
                                       gt_train = NULL) {
  
  cell_key <- as.character(id)
  
  max_H <- max_forecast_window(forecast_window, horizons)
  
  # Only unconditional runs get day_0
  include_day0 <- isFALSE(conditional)
  max_cols <- max_H + if (include_day0) 1L else 0L
  
  # Choose density source:
  # - non-CV: single prefit KDE from `kdes`
  # - CV:     per-year KDEs trained excluding that year
  dens_static <- NULL
  dens_by_year <- NULL
  
  if (!cv_by_year) {
    dens_static <- kdes[[cell_key]]
  } else {
    if (is.null(gt_train)) stop("cv_by_year=TRUE requires gt_train.", call. = FALSE)
    
    years_needed <- sort(unique(as.integer(format(issue_grid$time, "%Y"))))
    
    gt_id <- gt_train %>%
      filter(id == cell_key)
    
    dens_by_year <- purrr::map(
      years_needed,
      function(y) {
        x <- gt_id %>%
          filter(year != y) %>%
          pull(onset_day)
        fit_kde(x)
      }
    )
    names(dens_by_year) <- as.character(years_needed)
  }
  
  issue_grid %>%
    mutate(
      id  = cell_key,
      d0  = compute_d0(time, season_start_md),
      H   = purrr::map_int(time, ~ resolve_forecast_window_by_time(.x, forecast_window, horizons)),
      issue_year = as.integer(format(time, "%Y")),
      probs = purrr::pmap(list(d0, H, issue_year), function(d0_i, H_i, y_i) {
        out <- rep(NA_real_, max_cols)
        
        if (is.na(H_i) || H_i <= 0) return(out)
        
        dens_use <- if (!cv_by_year) {
          dens_static
        } else {
          dens_by_year[[as.character(y_i)]]
        }
        
        p <- predict_from_kde(
          dens_use, d0_i, H_i,
          conditional  = conditional,
          include_day0 = include_day0
        )
        
        # Fill from left; p is length H (conditional) or H+1 (unconditional+day0)
        out[seq_len(min(length(p), max_cols))] <- p[seq_len(min(length(p), max_cols))]
        out
      })
    ) %>%
    transmute(
      time,
      year  = as.integer(format(time, "%Y")),
      id,
      model = dplyr::case_when(
        conditional & !cv_by_year  ~ "clim_kde",
        !conditional & !cv_by_year ~ "clim_kde_unc",
        conditional & cv_by_year   ~ "clim_kde_cv",
        TRUE                       ~ "clim_kde_unc_cv"
      ),
      probs = probs
    ) %>%
    mutate(
      probs = purrr::map(probs, ~ {
        nm <- if (include_day0) {
          c("predicted_prob_day_0", paste0("predicted_prob_day_", seq_len(max_H)))
        } else {
          paste0("predicted_prob_day_", seq_len(max_H))
        }
        rlang::set_names(.x, nm)
      })
    ) %>%
    tidyr::unnest_wider(probs)
}

compute_all_forecasts <- function(gt_train, issue_grid, season_start_md,
                                  forecast_window, horizons,
                                  conditional = TRUE,
                                  cv_by_year = TRUE) {
  
  cells <- gt_train %>% distinct(id)
  
  kdes <- if (!cv_by_year) fit_kdes_by_cell(gt_train) else NULL
  
  res <- purrr::map_dfr(
    cells$id,
    ~ compute_forecasts_for_cell(
      .x, issue_grid, kdes,
      season_start_md, forecast_window, horizons,
      conditional = conditional,
      cv_by_year  = cv_by_year,
      gt_train    = gt_train
    )
  )
  
  list(forecasts = res, kdes = kdes)
}
