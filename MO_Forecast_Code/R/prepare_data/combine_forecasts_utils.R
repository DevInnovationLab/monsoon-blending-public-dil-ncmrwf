# ==============================================================================
# File: combine_forecasts_utils.R
#
# Purpose:
#   Utilities for combining multiple monsoon-onset forecast products into a
#   single, consistent *WIDE* schema suitable for merging and downstream output.
#
# Overview (current wide-only pipeline):
#   The driver script (e.g., 3_combine_datasets.R) needs to:
#     1) Read ground-truth onset (wide by year: lat/lon/year).
#     2) Read climatology issue forecasts (WIDE by issue date; day horizons as columns).
#     3) Read one or more rainfall-forecast products (stage-2 WIDE CSVs),
#        potentially from multiple sources and year subsets, then:
#          - extract/rename "constant" variables (one value per lat/lon/time/year)
#          - extract/rename "daily" variables stored as <var>_day_<k> columns
#          - optionally normalize/create a "plus" bin column: <prefix>_day_<max_day>plus
#     4) Return per-forecast-family outputs as:
#          - daily (WIDE) keyed by (lat, lon, time, year)
#          - constants (WIDE) keyed by (lat, lon, time, year)
#
# Conventions:
#   - Primary join keys for WIDE issue-date tables:
#       * (lat, lon, time, year)
#   - Ground-truth tables join on:
#       * (lat, lon, year)
#   - Daily horizon columns follow the pattern:
#       <prefix>_day_<k>
#     and an optional plus-bin column:
#       <prefix>_day_<max_day>plus
#     where <prefix> is typically either:
#       - out_prefix (for climatology), or
#       - "<forecast_name>_<out>" (for forecast families)
#
# Dependencies:
#   - data.table (fread, rbindlist, merge, setkeyv, setnames)
#   - A null-coalescing operator %||% is assumed to exist in the caller
#     (or define: `%||% <- function(x,y) if (!is.null(x)) x else y`).
#
# Functions:
#
#   parse_years(x)
#     Parse a "years" spec string into an integer vector.
#     Accepts comma-separated years and ranges like "2001,2003:2007".
#     Returns unique years; empty/NULL/NA returns integer(0) (meaning no filter).
#
#   detect_plus_bin(day_cols)
#     Given a character vector of day-column names, find columns that end with "plus".
#     Errors if more than one is found. Returns character(0) if none exist.
#
#   normalize_plus_bin(dt, prefix)
#     Ensure a WIDE plus-bin column exists and is consistent with other day columns
#     for a given prefix.
#       - Finds day columns matching `^<prefix>_day_`
#       - If no plus column exists, infers max day from numeric suffixes and creates
#         `<prefix>_day_<max_day>plus` initialized to NA_real_
#       - Sets the plus-bin value row-wise to `pmax(0, 1 - rowSums(days, na.rm=TRUE))`
#         using all non-plus day columns.
#     Returns the modified data.table.
#
#   read_ground_truth_wide(path)
#     Read stage-2 ground-truth wide CSV (requires lat, lon, year).
#     If present, maps:
#       mr_onset_day  -> true_onset_day (numeric)
#       mr_onset_date -> true_onset_date (Date)
#     If missing, fills with NA. Returns:
#       (lat, lon, year, true_onset_day, true_onset_date)
#
#   read_and_format_climatology_wide(path, out_prefix = "clim_p_onset")
#     Read climatology WIDE CSV (requires lat, lon, time, year).
#     Must NOT contain a `day` column.
#     Expects probability columns named `predicted_prob_day_<k>`, coerces them numeric,
#     renames to `<out_prefix>_day_<k>`, keeps only keys + renamed day columns, and
#     calls normalize_plus_bin(dt, out_prefix).
#     Returns a keyed WIDE data.table.
#
#   read_one_forecast_source(source, daily_cols, const_cols)
#     Read a single forecast WIDE CSV at `source$file` and optionally filter by years
#     from `source$years` (via parse_years()).
#     Requires key columns (lat, lon, time, year), constant columns (const_cols), and
#     for each daily var in daily_cols, columns matching `<var>_day_<k>`.
#     Keeps only keys + constants + all matched day columns and coerces day columns
#     numeric. Returns a WIDE data.table.
#
#   format_forecast_family(forecast_name, conf, plus_label_template)
#     Read+bind all WIDE sources for one forecast family and produce:
#       - $daily: WIDE table keyed by (lat, lon, time, year) with columns
#           `<forecast_name>_<out>_day_<k>` for each daily spec in conf$daily,
#           and (if add_plus=TRUE) a normalized `<...>_day_<max_day>plus` column.
#       - $constants: WIDE table keyed by (lat, lon, time, year) with columns
#           `<forecast_name>_<out>` for each constant spec in conf$constants.
#     Notes:
#       - Requires conf$max_day > 0 and conf$daily non-empty.
#       - Errors if any bound source produces a `day` column (wide-only enforcement).
#       - If duplicate (lat, lon, time, year) rows occur across sources, keeps the
#         first occurrence (with a warning).
#
# ==============================================================================
#
# NOTE (change):
#   This version uses a single column `id` as the grid-cell key (instead of lat/lon).
#   Primary join keys for WIDE issue-date tables:
#     * (id, time, year)
#   Ground-truth tables join on:
#     * (id, year)
# ==============================================================================

parse_years <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(integer(0))
  x <- gsub("\\s+", "", as.character(x))
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  out <- integer(0)
  for (p in parts) {
    if (p == "") next
    if (grepl(":", p, fixed = TRUE)) {
      rr <- unlist(strsplit(p, ":", fixed = TRUE))
      if (length(rr) != 2) stop("Bad years range: ", p, call. = FALSE)
      a <- as.integer(rr[1]); b <- as.integer(rr[2])
      if (is.na(a) || is.na(b)) stop("Bad years range: ", p, call. = FALSE)
      out <- c(out, seq.int(a, b))
    } else {
      y <- as.integer(p)
      if (is.na(y)) stop("Bad year: ", p, call. = FALSE)
      out <- c(out, y)
    }
  }
  unique(out)
}

# ----------------------------
# Plus-bin normalization
# ----------------------------
detect_plus_bin <- function(day_cols) {
  plus <- grep("plus$", day_cols, value = TRUE)
  if (length(plus) > 1L)
    stop("Multiple plus-bin columns found: ",
         paste(plus, collapse = ", "),
         call. = FALSE)
  plus
}

normalize_plus_bin <- function(dt, prefix) {
  day_pat <- paste0("^", prefix, "_day_")
  day_cols <- grep(day_pat, names(dt), value = TRUE)
  
  if (!length(day_cols))
    return(dt)
  
  plus_col <- detect_plus_bin(day_cols)
  
  if (length(plus_col) == 0L) {
    nums <- suppressWarnings(as.integer(sub(".*_day_", "", day_cols)))
    nums <- nums[!is.na(nums)]
    if (!length(nums))
      stop("Cannot infer day numbers for prefix ", prefix, call. = FALSE)
    
    max_day <- max(nums)
    plus_col <- paste0(prefix, "_day_", max_day + 1, "plus")
    dt[, (plus_col) := NA_real_]
    day_cols <- c(day_cols, plus_col)
  }
  
  base_cols <- setdiff(day_cols, plus_col)
  
  dt[, (plus_col) := {
    s <- rowSums(.SD, na.rm = TRUE)
    pmax(0, 1 - s)
  }, .SDcols = base_cols]
  
  dt[]
}

# ------------------------------------------------------------------------------
# Ground truth (id-based)
# ------------------------------------------------------------------------------

read_ground_truth_wide <- function(path) {
  dt <- data.table::as.data.table(readRDS(path))
  
  req <- c("id", "year")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0L)
    stop("Ground truth missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  
  true_onset_day  <- if ("mr_onset_day"  %in% names(dt)) as.numeric(dt[["mr_onset_day"]])  else NA_real_
  true_onset_date <- if ("mr_onset_date" %in% names(dt)) as.Date(dt[["mr_onset_date"]])   else as.Date(NA)
  
  data.table(
    id             = as.character(dt[["id"]]),
    year           = as.integer(dt[["year"]]),
    true_onset_day = true_onset_day,
    true_onset_date = true_onset_date
  )
}

# ------------------------------------------------------------------------------
# Climatology WIDE (id-based)
# ------------------------------------------------------------------------------

read_and_format_climatology_wide <- function(
    path,
    out_prefix = "clim_p_onset"
) {
  if (!file.exists(path))
    stop("Climatology file not found: ", path, call. = FALSE)

  dt <- data.table::as.data.table(readRDS(path))

  # ---- required keys
  key_req <- c("id", "time", "year")
  miss_keys <- setdiff(key_req, names(dt))
  if (length(miss_keys) > 0L) {
    stop(
      "Climatology missing key columns: ",
      paste(miss_keys, collapse = ", "),
      call. = FALSE
    )
  }
  
  # ---- disallow long climatology
  if ("day" %in% names(dt)) {
    stop(
      "Climatology must be WIDE (no `day` column). ",
      "Found `day` in: ", path,
      call. = FALSE
    )
  }
  
  # ---- find wide probability columns
  prob_cols <- grep("^predicted_prob_day_\\d+$", names(dt), value = TRUE)
  k <- as.integer(sub("^predicted_prob_day_", "", prob_cols))
  prob_cols <- prob_cols[order(k)]
  if (!length(prob_cols)) {
    stop(
      "No predicted_prob_day_<k> columns found in climatology: ",
      path,
      call. = FALSE
    )
  }
  
  # ---- coerce key types
  dt[, `:=`(
    id   = as.character(id),
    time = as.Date(time),
    year = as.integer(year)
  )]
  
  # ---- coerce probability columns
  for (cc in prob_cols) {
    dt[, (cc) := suppressWarnings(as.numeric(get(cc)))]
  }
  
  # ---- rename to combine-friendly names
  new_names <- sub(
    "^predicted_prob_day_",
    paste0(out_prefix, "_day_"),
    prob_cols
  )
  data.table::setnames(dt, prob_cols, new_names)
  
  # ---- keep only what Stage 4 needs
  keep <- c(key_req, new_names)
  dt <- dt[, ..keep]
  
  dt <- normalize_plus_bin(dt, out_prefix)
  
  data.table::setkeyv(dt, c("id", "time", "year"))
  dt[]
}

# ------------------------------------------------------------------------------
# Forecast sources (id-based)
# ------------------------------------------------------------------------------

read_one_forecast_source <- function(source, daily_cols, const_cols) {
  file_path <- as.character(source$file)
  if (!file.exists(file_path)) stop("Forecast file not found: ", file_path, call. = FALSE)

  yrs <- parse_years(source$years %||% "")

  dt <- data.table::as.data.table(readRDS(file_path))

  key_req <- c("id", "time", "year")
  miss_keys <- setdiff(key_req, names(dt))
  if (length(miss_keys) > 0L) {
    stop("Forecast file missing key columns: ", paste(miss_keys, collapse = ", "),
         "\nFile: ", file_path, call. = FALSE)
  }
  
  dt[, id := as.character(id)]
  dt[, time := as.Date(time)]
  dt[, year := as.integer(year)]
  if (length(yrs) > 0L) dt <- dt[year %in% yrs]
  
  # WIDE format: must have constants and <daily>_day_<n> columns
  miss_const <- setdiff(const_cols, names(dt))
  if (length(miss_const) > 0L) {
    stop("Forecast wide file missing constant columns: ", paste(miss_const, collapse = ", "),
         "\nFile: ", file_path, call. = FALSE)
  }
  
  # validate wide day columns exist for each daily var
  for (dv in daily_cols) {
    pat <- paste0("^", dv, "_day_\\d+$")
    dcols <- grep(pat, names(dt), value = TRUE)
    if (!length(dcols)) {
      stop("Forecast wide file has no day columns for daily var '", dv,
           "'. Expected columns like ", dv, "_day_<n>.\nFile: ", file_path, call. = FALSE)
    }
  }
  
  # keep only keys + constants + all wide day columns
  all_day_cols <- unlist(lapply(
    daily_cols,
    function(dv) grep(paste0("^", dv, "_day_\\d+$"), names(dt), value = TRUE)
  ))
  keep <- unique(c(key_req, const_cols, all_day_cols))
  dt <- dt[, ..keep]
  
  # coerce day cols numeric
  for (cc in all_day_cols) dt[, (cc) := suppressWarnings(as.numeric(get(cc)))]
  
  dt[]
}

# ------------------------------------------------------------------------------
# Forecast family formatter (id-based)
# ------------------------------------------------------------------------------

format_forecast_family <- function(forecast_name, conf, plus_label_template) {
  
  max_day <- as.integer(conf$max_day)
  if (is.na(max_day) || max_day <= 0)
    stop("Forecast '", forecast_name, "' must define a positive max_day.", call. = FALSE)
  
  if (is.null(conf$daily) || length(conf$daily) == 0L)
    stop("Forecast '", forecast_name, "' must provide a non-empty daily list.", call. = FALSE)
  
  daily <- lapply(conf$daily, function(x) {
    list(
      col      = as.character(x$col),
      out      = as.character(x$out),
      add_plus = isTRUE(x$add_plus)
    )
  })
  
  constants <- lapply(conf$constants %||% list(), function(x) {
    list(col = as.character(x$col), out = as.character(x$out))
  })
  
  cols_daily <- unique(vapply(daily, `[[`, character(1), "col"))
  cols_const <- unique(vapply(constants, `[[`, character(1), "col"))
  
  # ----------------------------
  # Read WIDE sources only
  # ----------------------------
  parts <- lapply(
    conf$sources,
    read_one_forecast_source,
    daily_cols = cols_daily,
    const_cols = cols_const
  )
  
  dt <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  
  if ("day" %in% names(dt))
    stop("Forecast '", forecast_name, "' produced a `day` column; forecasts must be wide-only.",
         call. = FALSE)
  
  key_cols <- c("id", "time", "year")
  data.table::setkeyv(dt, key_cols)
  
  # ----------------------------
  # Overlap guard (id/time/year)
  # ----------------------------
  if (dt[, anyDuplicated(paste(id, time, year))] > 0) {
    warning("Forecast '", forecast_name,
            "' has overlapping rows across sources; keeping first occurrence.")
    dt <- unique(dt, by = key_cols)
  }
  
  # ----------------------------
  # Constants (wide)
  # ----------------------------
  if (length(constants)) {
    const_dt <- unique(dt[, c(key_cols, vapply(constants, `[[`, "", "col")), with = FALSE])
    for (x in constants) {
      data.table::setnames(
        const_dt,
        x$col,
        paste0(forecast_name, "_", x$out)
      )
    }
  } else {
    const_dt <- unique(dt[, ..key_cols])
  }
  
  data.table::setkeyv(const_dt, key_cols)
  
  # ----------------------------
  # Daily outputs (wide)
  # ----------------------------
  daily_wide <- unique(dt[, ..key_cols])
  data.table::setkeyv(daily_wide, key_cols)
  
  for (d in daily) {
    src <- d$col
    out <- d$out
    
    pat <- paste0("^", src, "_day_")
    dcols <- grep(pat, names(dt), value = TRUE)
    if (!length(dcols))
      stop("Forecast '", forecast_name, "' missing day columns for ", src, call. = FALSE)
    
    # horizon check (only plain numeric suffixes)
    bad <- grep(paste0("_day_(\\d+)$"), dcols, value = TRUE)
    if (length(bad)) {
      k <- as.integer(sub(".*_day_", "", bad))
      if (any(k > max_day))
        stop("Forecast '", forecast_name, "' exceeds max_day for ", src, call. = FALSE)
    }
    
    tmp <- dt[, c(key_cols, dcols), with = FALSE]
    
    new_names <- sub(
      paste0("^", src, "_day_"),
      paste0(forecast_name, "_", out, "_day_"),
      dcols
    )
    data.table::setnames(tmp, dcols, new_names)
    
    daily_wide <- merge(daily_wide, tmp, by = key_cols, all = TRUE)
    
    # plus-bin normalization (if requested)
    if (isTRUE(d$add_plus)) {
      daily_wide <- normalize_plus_bin(
        daily_wide,
        paste0(forecast_name, "_", out)
      )
    }
  }
  
  list(
    daily     = daily_wide,
    constants = const_dt
  )
}
