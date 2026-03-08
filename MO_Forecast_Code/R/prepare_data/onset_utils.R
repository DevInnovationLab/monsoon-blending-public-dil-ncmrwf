# ==============================================================================
# File: onset_utils.R
#
# Purpose:
#   Low-level utilities for computing monsoon onset indices from
#   rainfall time series.
#
# Overview:
#   This file provides fast, side-effect–free functions used by both
#   forecast and ground-truth pipelines to identify the first
#   "onset day" satisfying rainfall accumulation and follow-up rules.
#
# Functions:
#
#   read_mok_dates(spec)
#     - Reads and prepares Monsoon Onset Key (MOK) dates when
#       spec$screen$mode == "MOK".
#     - Returns a data.table with columns:
#         * year
#         * mok_date
#     - Returns NULL when MOK screening is disabled.
#
#   read_thresholds(spec)
#     - Reads grid-cell–specific onset thresholds from
#       spec$thresholds$file.
#     - Standardizes column names to:
#         * lat, lon, onset_thresh
#     - Returns NULL if no threshold file is specified.
#
#   .roll_sum_na_rm_left(x, k)
#     - Internal helper.
#     - Computes left-aligned rolling sums of length k with NA removed
#       (NA treated as zero).
#
#   .roll_sum_na_propagate_left(x, k)
#     - Internal helper.
#     - Computes left-aligned rolling sums of length k where any NA
#       in the window propagates to the sum.
#
#   find_onset(series, window, thresh,
#              reject_if_short_followup = FALSE,
#              start_day = 0)
#     - Core onset-detection routine.
#     - Returns the index of the earliest window satisfying:
#         * point rainfall > 1
#         * rolling sum over `window` exceeds `thresh`
#         * optional follow-up constraints (10-day / 30-day rules)
#         * optional start-day restriction (e.g., June or MOK)
#     - Returns NA_real_ if no valid onset is found.
#
# Notes:
#   - All functions are pure (no I/O except CSV reads in read_* helpers)
#   - Designed for O(n) performance on long daily series
#   - All configuration is passed explicitly via `spec`
#
# Used By:
#   - 1_process_raw_nc_files.R
#
# ==============================================================================


read_mok_dates <- function(spec) {
  if (is.null(spec$mok$file)) stop("screen$mode == 'MOK' requires spec$mok$$mok_file.")
  if (!file.exists(spec$mok$file)) stop("MOK file not found: ", spec$mok$file)
  
  mok <- data.table::fread(spec$mok$file)
  ycol <- spec$mok$year_col
  dcol <- spec$mok$day_col
  if (!(ycol %in% names(mok)) || !(dcol %in% names(mok))) {
    stop("MOK file must contain columns ", ycol, " and ", dcol)
  }
  
  # Mok date is (year-<mok_base_date>) + mok_day
  base_md <- spec$mok$base_date
  mok[, mok_date := as.Date(paste0(get(ycol), "-", base_md)) + as.integer(get(dcol))]
  mok[, .(year = as.integer(get(ycol)), mok_date)]
}
read_thresholds <- function(spec) {
  # spec$thresholds can be missing entirely
  if (is.null(spec$thresholds$file)) return(NULL)
  
  f <- spec$thresholds$file
  if (!file.exists(f)) stop("Threshold file not found: ", f, call. = FALSE)
  
  ext <- tolower(tools::file_ext(f))
  
  # Helper: normalize to data.table(lat, lon, onset_thresh)
  normalize_thr <- function(lat, lon, onset_thresh) {
    lat <- as.vector(lat)
    lon <- as.vector(lon)
    
    # onset_thresh might be:
    # - vector length = length(lat)*length(lon), OR
    # - matrix [lat x lon] or [lon x lat]
    ot <- onset_thresh
    
    if (is.matrix(ot) || length(dim(ot)) == 2L) {
      d <- dim(ot)
      # try interpret as [lat x lon]
      if (identical(d, c(length(lat), length(lon)))) {
        grid <- CJ(lat = lat, lon = lon, unique = TRUE)
        grid[, onset_thresh := as.vector(ot)]
        return(unique(grid))
      }
      # try interpret as [lon x lat]
      if (identical(d, c(length(lon), length(lat)))) {
        grid <- CJ(lat = lat, lon = lon, unique = TRUE)
        grid[, onset_thresh := as.vector(t(ot))]
        return(unique(grid))
      }
      # fall back: just vectorize
      ot <- as.vector(ot)
    } else {
      ot <- as.vector(ot)
    }
    
    if (length(ot) == length(lat) * length(lon)) {
      grid <- CJ(lat = lat, lon = lon, unique = TRUE)
      grid[, onset_thresh := ot]
      return(unique(grid))
    }
    
    stop(
      "Threshold array shape doesn't match lat/lon lengths. ",
      "length(lat)=", length(lat), ", length(lon)=", length(lon),
      ", length(onset_thresh)=", length(ot),
      call. = FALSE
    )
  }
  
  if (ext == "mat") {
    # ---- .mat thresholds (2deg style) ----
    if (!requireNamespace("R.matlab", quietly = TRUE)) {
      stop("Need package R.matlab to read .mat thresholds. Install with install.packages('R.matlab').", call. = FALSE)
    }
    
    md <- R.matlab::readMat(f)
    
    # Allow override of object names via spec, else use your historic names
    lat_name <- spec$thresholds$mat_lat_name %||% "lat"
    lon_name <- spec$thresholds$mat_lon_name %||% "lon"
    thr_name <- spec$thresholds$mat_thresh_name %||% "onset.day.thres"
    
    if (is.null(md[[lat_name]]) || is.null(md[[lon_name]]) || is.null(md[[thr_name]])) {
      stop(
        "MAT file missing one of: ", lat_name, ", ", lon_name, ", ", thr_name,
        ". Available objects: ", paste(names(md), collapse = ", "),
        call. = FALSE
      )
    }
    
    mat_lat <- as.vector(md[[lat_name]])
    mat_lon <- as.vector(md[[lon_name]])
    onset_thresh <- md[[thr_name]]
    
    return(normalize_thr(mat_lat, mat_lon, onset_thresh))
  }
  
  if (ext %in% c("nc", "nc4", "netcdf")) {
    # ---- NetCDF thresholds (quarter-degree style) ----
    nc <- ncdf4::nc_open(f)
    on.exit(ncdf4::nc_close(nc), add = TRUE)
    
    # Allow override via spec; else use your historic/default names
    lat_name <- spec$thresholds$lat_var %||% spec$thresholds$lat_col %||% "lat"
    lon_name <- spec$thresholds$lon_var %||% spec$thresholds$lon_col %||% "lon"
    thr_name <- spec$thresholds$thresh_var %||% spec$thresholds$thresh_col %||% "MWmean"
    
    # Use robust coord getter if the vars aren’t present under those names
    mat_lat <- tryCatch(ncdf4::ncvar_get(nc, lat_name), error = function(e) NULL)
    if (is.null(mat_lat)) mat_lat <- get_nc_coord(nc, "lat")
    
    mat_lon <- tryCatch(ncdf4::ncvar_get(nc, lon_name), error = function(e) NULL)
    if (is.null(mat_lon)) mat_lon <- get_nc_coord(nc, "lon")
    
    if (!(thr_name %in% names(nc$var))) {
      stop("Threshold NetCDF missing variable '", thr_name, "'. Vars: ", paste(names(nc$var), collapse = ", "), call. = FALSE)
    }
    onset_thresh <- ncdf4::ncvar_get(nc, thr_name)
    
    return(normalize_thr(mat_lat, mat_lon, onset_thresh))
  }
  
  # ---- CSV/TSV thresholds (your current style) ----
  th <- data.table::fread(f)
  
  lc <- spec$thresholds$lat_col %||% "lat"
  oc <- spec$thresholds$lon_col %||% "lon"
  tc <- spec$thresholds$thresh_col %||% "onset_thresh"
  
  if (!(lc %in% names(th)) || !(oc %in% names(th)) || !(tc %in% names(th))) {
    stop("Threshold file must contain columns ", lc, ", ", oc, ", ", tc, call. = FALSE)
  }
  
  data.table::setnames(th, c(lc, oc, tc), c("lat", "lon", "onset_thresh"))
  unique(th[, .(lat, lon, onset_thresh)])
}


# ----------------------------
# Fast rolling helpers (base R, O(n))
# ----------------------------

.roll_sum_na_rm_left <- function(x, k) {
  # left-aligned rolling sum with na.rm=TRUE (length n-k+1)
  n <- length(x)
  if (k <= 0L || n < k) return(numeric(0))
  x0 <- x
  x0[is.na(x0)] <- 0
  cs <- c(0, cumsum(x0))
  cs[(k + 1):(n + 1)] - cs[1:(n - k + 1)]
}

.roll_sum_na_propagate_left <- function(x, k) {
  # left-aligned rolling sum where any NA in window => NA (length n-k+1)
  n <- length(x)
  if (k <= 0L || n < k) return(numeric(0))
  na <- is.na(x)
  x0 <- x
  x0[na] <- 0

  cs  <- c(0, cumsum(x0))
  cna <- c(0, cumsum(na))

  s <- cs[(k + 1):(n + 1)] - cs[1:(n - k + 1)]
  na_ct <- cna[(k + 1):(n + 1)] - cna[1:(n - k + 1)]
  s[na_ct > 0] <- NA_real_
  s
}

# ----------------------------
# Unified onset finder
# ----------------------------

find_onset <- function(series,
                       win = 5,
                       thresh,
                       reject_if_short_followup = FALSE,
                       start_day = 0) {
  n <- length(series)
  if (n < win || is.na(thresh)) return(NA_real_)

  max_candidate <- n - win + 1L
  if (max_candidate < 1L) return(NA_real_)

  # start_day: reject any candidate indices < start_day
  # (negative start_day => no restriction)
  min_i <- as.integer(max(1L, ceiling(start_day)))
  if (min_i > max_candidate) return(NA_real_)

  idx <- min_i:max_candidate

  # base window sums (na.rm=TRUE, like your original find_onset)
  wsum_all <- .roll_sum_na_rm_left(series, win)   # length max_candidate
  wsum <- wsum_all[idx]

  base_ok <- (series[idx] > 1) & (wsum > thresh)
  base_ok[is.na(base_ok)] <- FALSE
  cand <- idx[base_ok]
  if (!length(cand)) return(NA_real_)

  # Follow-up rule:
  # If follow length >= 10, require no 10-day sum < 5 within follow segment,
  # where 10-day sums propagate NA (like zoo::rollapply(sum) default), and we ignore NA.
  if (n < 10L) return(cand[1])

  sum10 <- .roll_sum_na_propagate_left(series, 10L)   # length n-9
  bad10 <- !is.na(sum10) & (sum10 < 5)
  pre_bad <- c(0L, cumsum(as.integer(bad10)))
  last10start <- n - 10L + 1L

  if (reject_if_short_followup) {
    # require full 30 days follow-up
    cand <- cand[cand + 30L <= n]
    if (!length(cand)) return(NA_real_)

    end_start <- cand + 21L  # last 10-window start within [i, i+30]
    has_bad <- (pre_bad[end_start + 1L] - pre_bad[cand]) > 0L
    ok <- !has_bad
    if (any(ok)) cand[which(ok)[1]] else NA_real_
  } else {
    # allow truncated follow-up; still apply the 10-day check when follow length >= 10
    follow_end <- pmin(n, cand + 30L)
    follow_len <- follow_end - cand + 1L

    # only candidates with follow_len >= 10 get checked; others pass
    needs_check <- follow_len >= 10L
    ok <- rep(TRUE, length(cand))

    if (any(needs_check)) {
      c2 <- cand[needs_check]
      fe <- follow_end[needs_check]
      end_start <- pmin(last10start, fe - 10L + 1L)  # last 10-window start in follow segment

      has_bad <- (pre_bad[end_start + 1L] - pre_bad[c2]) > 0L
      ok[needs_check] <- !has_bad
    }

    if (any(ok)) cand[which(ok)[1]] else NA_real_
  }
}
find_onset_precomp <- function(series, win, thresh,
                               wsum_all, pre_bad, last10start,
                               reject_if_short_followup = FALSE,
                               start_day = 0) {
  n <- length(series)
  if (n < win || is.na(thresh)) return(NA_real_)
  
  max_candidate <- n - win + 1L
  if (max_candidate < 1L) return(NA_real_)
  
  min_i <- as.integer(max(1L, ceiling(start_day)))
  if (min_i > max_candidate) return(NA_real_)
  
  idx <- min_i:max_candidate
  wsum <- wsum_all[idx]
  
  base_ok <- (series[idx] > 1) & (wsum > thresh)
  base_ok[is.na(base_ok)] <- FALSE
  cand <- idx[base_ok]
  if (!length(cand)) return(NA_real_)
  
  if (n < 10L) return(cand[1])
  
  if (reject_if_short_followup) {
    cand <- cand[cand + 30L <= n]
    if (!length(cand)) return(NA_real_)
    
    end_start <- cand + 21L
    has_bad <- (pre_bad[end_start + 1L] - pre_bad[cand]) > 0L
    ok <- !has_bad
    if (any(ok)) cand[which(ok)[1]] else NA_real_
  } else {
    follow_end <- pmin(n, cand + 30L)
    follow_len <- follow_end - cand + 1L
    
    needs_check <- follow_len >= 10L
    ok <- rep(TRUE, length(cand))
    
    if (any(needs_check)) {
      c2 <- cand[needs_check]
      fe <- follow_end[needs_check]
      end_start <- pmin(last10start, fe - 10L + 1L)
      has_bad <- (pre_bad[end_start + 1L] - pre_bad[c2]) > 0L
      ok[needs_check] <- !has_bad
    }
    
    if (any(ok)) cand[which(ok)[1]] else NA_real_
  }
}


