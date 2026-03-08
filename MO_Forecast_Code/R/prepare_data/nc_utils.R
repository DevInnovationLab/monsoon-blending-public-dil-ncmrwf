# ==============================================================================
# File: nc_utils.R
# ==============================================================================
# Purpose
#   Shared utilities for the *single-pass* NetCDF -> onset-output pipeline used by
#   1_process_raw_nc_files.R.
#
#   This file supports BOTH pipeline modes:
#     - type: "rainfall_forecast"
#         * reads each NetCDF forecast file into a *wide-by-lead-day* table
#           (memory-safe vs a huge long table)
#         * optionally transforms rainfall from a source grid to a target grid
#           using a linear-combination weights file (by cell id)
#         * runs onset detection per ensemble member (if 'number' exists), then
#           aggregates to ensemble summaries
#
#     - type: "ground_truth_rainfall"
#         * reads each NetCDF observed rainfall file into a *long* table
#           (time x cell)
#         * optionally transforms rainfall from a source grid to a target grid
#           using the same weights-file interface
#         * runs onset detection per (id, year)
#
# Key conventions
#   - Spatial key used throughout is `id`, defined as:
#       id = paste0(lat, "_", lon)
#     (created by add_id_from_latlon() when needed).
#
#   - Thresholds are joined by `id`:
#       * If the threshold file has lat/lon only, we derive id from them.
#       * If a cell-transform is enabled, thresholds must correspond to the
#         *target* grid ids (since onset is computed after transformation).
#
#   - Years are obtained from filenames (not from NetCDF contents):
#       list_nc_files_with_year() extracts a 4-digit year from each filename and
#       run_single_pipeline() overwrites/sets dt$year to this value.
#
#   - Latitude/longitude coordinate names can vary across NetCDF sources.
#     get_nc_coord() quietly tries common variants for lat/lon as variables, then
#     as dimensions, avoiding noisy ncdf4 console output during failures.
#
# Dependencies / Imports
#   - data.table, dplyr, stringr, purrr, tidyr, ncdf4, lubridate
#   - Also assumes helper functions exist elsewhere:
#       * %||% (from R/_shared/misc.R)
#       * load_spec(), validate_spec() (from R/_shared/read_spec.R)
#       * find_onset(), read_mok_dates(), read_thresholds() (from onset_utils.R)
#
# ==============================================================================
# Function index (what each function does)
# ==============================================================================
#
# ---- Spec + naming helpers ---------------------------------------------------
# validate_spec_single(spec)
#   Validates the YAML spec for required sections and fields, and applies basic
#   type-specific checks:
#     - rainfall_forecast requires options: min_day, max_day, window
#     - ground_truth_rainfall requires options: window, cutoff_month_day
#
# get_value_var(spec)
#   Returns the NetCDF variable name to read (spec$input$value_col) and errors
#   if missing.
#
# rename_dimensions(df, rename_map)
#   Case-insensitive rename of dimension/coordinate column names according to
#   spec$dimensions$rename (mapping old -> new). Returns the renamed data frame.
#
# ---- Input file discovery ----------------------------------------------------
# list_nc_files_with_year(spec)
#   Lists NetCDF files in input.nc_folder filtered by input.file_regex, extracts
#   a 4-digit year from each filename, and returns a data.table:
#     nc_path, year
#   Applies optional year bounds options.min_year/max_year.
#
# ---- NetCDF coordinate + time utilities -------------------------------------
# get_nc_coord(nc, coord_type = c("lat","lon"))
#   Robust coordinate getter that tries common variable name variants quietly
#   (LATITUDE/latitude/lat/LAT and LONGITUDE/longitude/lon/LON/long/LONG).
#   If not found as variables, tries as dimensions (nc$dim[[name]]$vals).
#   Returns numeric vector; errors if nothing found.
#
# nc_time_to_posixct(time_num, time_units, tz = "UTC")
#   Converts numeric NetCDF time values to POSIXct using the "units" attribute
#   of the time coordinate (e.g., "days since 1900-01-01").
#
# ---- ID + thresholds ---------------------------------------------------------
# add_id_from_latlon(dt, lat_col="lat", lon_col="lon", id_col="id")
#   Adds an id column as "<lat>_<lon>" if not already present.
#
# prep_thresholds_id(thr_dt)
#   Normalizes a threshold table to be keyed by (id, onset_thresh):
#     - lowercases column names
#     - if id missing but lat/lon present, derives id
#     - ensures onset_thresh is numeric
#     - returns unique (id, onset_thresh) with data.table key set on id
#
# attach_thresholds_id(dt, thr_dt)
#   Ensures dt has id, then left-joins onset_thresh from thr_dt by id.
#   If thr_dt is NULL, adds onset_thresh := NA_real_.
#
# ---- Optional: grid/cell transform (linear combination by id) ----------------
# read_cell_transform(spec)
#   Reads a long-form weights file when options.cell_transform_enabled is TRUE.
#   Expects columns: target_id, source_id, weight.
#   Returns a data.table or NULL if disabled.
#
# transform_forecast_wide(dt, weights_dt)
#   Applies the weights transform to forecast *wide-by-day* data.
#   Input must contain day columns named "*_day_<n>".
#   Produces rainfall on target_id cells (replacing id := target_id), preserving
#   other key columns, then casts back to wide format.
#
# transform_groundtruth_long(dt, weights_dt, value_col)
#   Applies the weights transform to ground-truth *long* data for the given value
#   column. Produces rainfall on target_id cells (replacing id := target_id).
#
# ---- NetCDF readers (stage-1 replacement) -----------------------------------
# nc_read_forecast_wide(nc_path, var_name, dim_rename_map, spec, day_dim="day",
#                       prefix=NULL, add_year=TRUE)
#   Reads a single forecast NetCDF into a wide-by-lead-day data.table:
#     - identifies coordinates/dim values (with robust lat/lon handling)
#     - applies dimension renaming (spec$dimensions$rename)
#     - reshapes the variable array so day is the column dimension
#     - filters days by options.min_day/options.max_day
#     - names day columns as "<prefix>_day_<day>"
#     - converts time to Date + adds year (if time exists and add_year=TRUE)
#
# nc_read_groundtruth_long(nc_path, var_name, dim_rename_map, add_year=TRUE)
#   Reads a single ground-truth NetCDF into a long data.table:
#     - expands the full dimension grid
#     - flattens the variable array into one value column
#     - applies dimension renaming and lowercasing
#     - drops NA values
#     - converts time if numeric and adds year (if add_year=TRUE)
#
# ---- Stage-2 onset helpers (id-keyed) ---------------------------------------
# order_day_cols(wide_dt, key_cols)
#   Given a wide forecast table where lead-day columns are numeric names
#   ("1","2",...), returns them ordered along with their integer values.
#
# .first_pos_ge_date(cal_dates, cutoff_date)
#   Returns first index where cal_dates >= cutoff_date, or 9999 if none.
#
# calc_member_onsets(row_dt, day_cols, day_ints, win, t0, yr, th, mk)
#   Computes onset indices for one ensemble member under three rules:
#     - raw: no start restriction
#     - june: start at June 2 of that year
#     - mok: start at mok_date (if provided)
#   Uses find_onset() from onset_utils.R and returns indices + the series.
#
# process_rainfall_forecast_id(dt, spec, mok_dt=NULL, thr_dt=NULL)
#   Forecast pipeline (id-keyed):
#     - expects wide-by-day columns "<wide_prefix>_day_<n>"
#     - renames day columns to numeric ("1","2",...)
#     - attaches thresholds (by id) and MOK (by year)
#     - computes per-member onset indices, then aggregates across members to
#       produce one row per (id, time, year, onset_thresh, mok_date) with:
#         * ensemble mean rainfall by day
#         * ensemble sd by day
#         * fraction raining by day (>1 threshold)
#         * onset probability mass function by day (raw/june/mok)
#         * rainfall quantiles by day (0, .2, .4, .6, .8, 1)
#
# process_ground_truth_rainfall_id(dt, spec, mok_dt=NULL, thr_dt=NULL, value_col)
#   Ground-truth pipeline (id-keyed):
#     - expects long daily series with columns: id, time, year, <value_col>
#     - attaches thresholds (by id) and MOK (by year)
#     - uses cutoff_month_day to restrict earliest candidate date each year
#     - computes mr_onset_idx, mr_onset_date, mr_onset_day, and an annotated
#       long table with mr_onset_flag
#
#
# ---- Pipeline entrypoint -----------------------------------------------------
# run_single_pipeline(spec_id)
#   Main driver used by single_pipeline_nc_to_onset.R:
#     - loads + validates YAML spec
#     - sets future plan (multisession or sequential)
#     - discovers files + years
#     - reads aux tables (MOK, thresholds) and optional weights
#     - loops by year and processes one year at a time
#     - appends outputs across years into final *_wide.csv (and *_long.csv for
#       ground truth)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(ncdf4)
  library(lubridate)
})

# ---------------------------------------------------------------------------
# Spec helpers
# ---------------------------------------------------------------------------

validate_spec_single <- function(spec) {
  spec <- validate_spec(
    spec,
    required_top = c("input", "dimensions", "output", "options", "type"),
    required_paths = c("input.nc_folder", "input.file_regex", "input.value_col",
                       "output.out_dir"),
    checks = list(function(s) {
      if (!dir.exists(s$input$nc_folder)) {
        stop("input$nc_folder does not exist: ", s$input$nc_folder, call. = FALSE)
      }
      
      if (s$type == "rainfall_forecast") {
        for (nm in c("min_day", "max_day", "window")) {
          if (is.null(s$options[[nm]])) stop("Missing options.", nm, call.=FALSE)
        }
      }
      
      if (s$type == "ground_truth_rainfall") {
        if (is.null(s$options$window)) stop("Missing options.window", call.=FALSE)
        if (is.null(s$options$cutoff_month_day)) stop("Missing options.cutoff_month_day", call.=FALSE)
      }
    }),
    label = "spec"
  )
  spec
}

get_value_var <- function(spec) {
  v <- spec$input$value_col
  if (!is.null(v) && nzchar(v)) return(as.character(v))
  stop("Missing spec$input$value_col in YAML.", call. = FALSE)
}

rename_dimensions <- function(df, rename_map) {
  if (length(rename_map) == 0) return(df)
  
  old <- names(rename_map)
  new <- unname(unlist(rename_map))
  
  df_names <- names(df)
  idx <- match(tolower(old), tolower(df_names))
  keep <- which(!is.na(idx))
  if (!length(keep)) return(df)
  
  from <- df_names[idx[keep]]
  to   <- new[keep]
  mapping <- stats::setNames(from, to) # new = old
  
  df %>% dplyr::rename(!!!mapping)
}

# ---------------------------------------------------------------------------
# Years from filenames
# ---------------------------------------------------------------------------

list_nc_files_with_year <- function(spec) {
  files <- list.files(spec$input$nc_folder, full.names = TRUE)
  files <- files[str_detect(basename(files), spec$input$file_regex)]
  if (!length(files)) stop("No .nc files matched file_regex in nc_folder.", call. = FALSE)
  
  yrs_chr <- stringr::str_extract(basename(files), "(19|20)\\d{2}")
  yrs <- suppressWarnings(as.integer(yrs_chr))
  
  if (anyNA(yrs)) {
    bad <- files[is.na(yrs)]
    stop(
      "Could not extract a 4-digit year from some filenames. First few:\n  - ",
      paste(head(bad, 5), collapse = "\n  - "),
      call. = FALSE
    )
  }
  
  dt <- data.table::data.table(nc_path = files, year = yrs)
  
  min_year <- spec$options$min_year %||% NA_integer_
  max_year <- spec$options$max_year %||% NA_integer_
  if (!is.na(min_year)) dt <- dt[year >= min_year]
  if (!is.na(max_year)) dt <- dt[year <= max_year]
  if (!nrow(dt)) stop("After min_year/max_year filtering, no files remain.", call. = FALSE)
  
  dt[]
}

# ---------------------------------------------------------------------------
# NetCDF coord getter (quiet + supports coords stored as dims)
# ---------------------------------------------------------------------------

get_nc_coord <- function(nc, coord_type = "lat") {
  variants <- if (coord_type == "lat") {
    c("LATITUDE", "latitude", "lat", "LAT")
  } else {
    c("LONGITUDE", "longitude", "lon", "LON", "long", "LONG")
  }
  
  quiet_ncvar_get <- function(nc, v) {
    tryCatch({
      zz <- tempfile()
      sink(zz)
      on.exit(sink(), add = TRUE)
      suppressWarnings(suppressMessages(ncvar_get(nc, v)))
    }, error = function(e) NULL)
  }
  
  # 1) try as variable (quiet)
  for (v in variants) {
    out <- quiet_ncvar_get(nc, v)
    if (!is.null(out)) return(as.vector(out))
  }
  
  # 2) try as dimension (common for rectilinear coords)
  for (v in variants) {
    if (v %in% names(nc$dim)) {
      d <- nc$dim[[v]]
      if (!is.null(d$vals) && length(d$vals) > 0) return(as.vector(d$vals))
    }
  }
  
  stop(paste("Could not find", coord_type, "variable/dimension in NetCDF file"), call. = FALSE)
}


# ---------------------------------------------------------------------------
# NetCDF coord getter (quietly checks different possible "time" variables (e.g. capitalization, etc. ))
# ---------------------------------------------------------------------------
get_nc_time <- function(nc) {
  variants <- c("TIME", "time", "Time", "t", "T")
  
  quiet_ncvar_get <- function(nc, v) {
    tryCatch({
      zz <- tempfile()
      sink(zz)
      on.exit(sink(), add = TRUE)
      suppressWarnings(suppressMessages(ncdf4::ncvar_get(nc, v)))
    }, error = function(e) NULL)
  }
  
  quiet_ncatt_get_units <- function(nc, v) {
    tryCatch({
      u <- ncdf4::ncatt_get(nc, v, "units")$value
      if (!is.null(u) && nzchar(u)) u else NULL
    }, error = function(e) NULL)
  }
  
  # 1) try as variable (quiet)
  for (v in variants) {
    x <- quiet_ncvar_get(nc, v)
    if (!is.null(x)) {
      u <- quiet_ncatt_get_units(nc, v)
      return(list(values = as.vector(x), units = u, name = v, source = "var"))
    }
  }
  
  # 2) try case-insensitive match among variables
  for (v in variants) {
    hit <- names(nc$var)[tolower(names(nc$var)) == tolower(v)][1]
    if (!is.na(hit) && nzchar(hit)) {
      x <- quiet_ncvar_get(nc, hit)
      if (!is.null(x)) {
        u <- quiet_ncatt_get_units(nc, hit)
        return(list(values = as.vector(x), units = u, name = hit, source = "var_ci"))
      }
    }
  }
  
  # 3) try as dimension values
  for (v in variants) {
    if (v %in% names(nc$dim)) {
      d <- nc$dim[[v]]
      if (!is.null(d$vals) && length(d$vals) > 0) {
        # units might still live on a var (same name, different case)
        u <- quiet_ncatt_get_units(nc, v)
        if (is.null(u)) {
          hit <- names(nc$var)[tolower(names(nc$var)) == tolower(v)][1]
          if (!is.na(hit) && nzchar(hit)) u <- quiet_ncatt_get_units(nc, hit)
        }
        return(list(values = as.vector(d$vals), units = u, name = v, source = "dim"))
      }
    }
  }
  
  stop("Could not find time variable/dimension in NetCDF file", call. = FALSE)
}

coerce_time_col <- function(dt, nc, time_col = "time") {
  if (!(time_col %in% names(dt))) return(dt)
  
  if (is.numeric(dt[[time_col]])) {
    tinfo <- get_nc_time(nc)
    if (is.null(tinfo$units) || !nzchar(tinfo$units)) {
      stop("Time column is numeric but NetCDF time 'units' attribute is missing; cannot convert.", call. = FALSE)
    }
    dt[, (time_col) := as.Date(nc_time_to_posixct(get(time_col), tinfo$units))]
  } else {
    dt[, (time_col) := as.Date(get(time_col))]
  }
  
  dt
}

# ---------------------------------------------------------------------------
# Time conversion
# ---------------------------------------------------------------------------

nc_time_to_posixct <- function(time_num, time_units, tz = "UTC") {
  if (is.null(time_units) || !nzchar(time_units)) {
    stop("NetCDF time variable is numeric but has no 'units' attribute.", call. = FALSE)
  }
  
  m <- regexec("^\\s*(seconds|second|minutes|minute|hours|hour|days|day)\\s+since\\s+(.+)\\s*$",
               time_units, ignore.case = TRUE)
  parts <- regmatches(time_units, m)[[1]]
  if (length(parts) < 3) stop("Unrecognized NetCDF time units format: ", time_units, call. = FALSE)
  
  unit <- tolower(parts[2])
  origin_str <- parts[3]
  
  origin <- as.POSIXct(origin_str, tz = tz, format = "%Y-%m-%d %H:%M:%S")
  if (is.na(origin)) origin <- as.POSIXct(origin_str, tz = tz, format = "%Y-%m-%d %H:%M")
  if (is.na(origin)) origin <- as.POSIXct(origin_str, tz = tz, format = "%Y-%m-%d")
  if (is.na(origin)) stop("Could not parse NetCDF time origin: ", origin_str, call. = FALSE)
  
  mult <- switch(unit,
                 "second" = 1, "seconds" = 1,
                 "minute" = 60, "minutes" = 60,
                 "hour"   = 3600, "hours"   = 3600,
                 "day"    = 86400, "days"    = 86400,
                 stop("Unsupported time unit: ", unit, call. = FALSE)
  )
  
  origin + time_num * mult
}

# ---------------------------------------------------------------------------
# id helpers
# ---------------------------------------------------------------------------

add_id_from_latlon <- function(dt, lat_col = "lat", lon_col = "lon", id_col = "id") {
  dt <- copy(dt)
  if (id_col %in% names(dt)) return(dt)
  
  if (!all(c(lat_col, lon_col) %in% names(dt))) {
    stop("Cannot create id: missing columns ", lat_col, " and/or ", lon_col, call. = FALSE)
  }
  dt[, (id_col) := paste0(
    format(get(lat_col), scientific = FALSE, trim = TRUE),
    "_",
    format(get(lon_col), scientific = FALSE, trim = TRUE)
  )]
  dt
}

# ---------------------------------------------------------------------------
# Thresholds keyed by id
# ---------------------------------------------------------------------------

prep_thresholds_id <- function(thr_dt) {
  if (is.null(thr_dt)) return(NULL)
  thr_dt <- copy(thr_dt)
  
  if (!("onset_thresh" %in% names(thr_dt))) {
    stop("Thresholds table must contain onset_thresh.", call. = FALSE)
  }
  
  nms <- names(thr_dt)
  setnames(thr_dt, nms, tolower(nms))
  
  if (!("id" %in% names(thr_dt))) {
    if (!all(c("lat", "lon") %in% names(thr_dt))) {
      stop("Thresholds table must have either id or (lat, lon).", call. = FALSE)
    }
    thr_dt <- add_id_from_latlon(thr_dt, lat_col = "lat", lon_col = "lon", id_col = "id")
  }
  
  thr_dt[, onset_thresh := as.numeric(onset_thresh)]
  thr_dt <- unique(thr_dt[, .(id, onset_thresh)])
  setkey(thr_dt, id)
  thr_dt
}

attach_thresholds_id <- function(dt, thr_dt) {
  dt <- copy(dt)
  dt <- add_id_from_latlon(dt)
  
  if (is.null(thr_dt)) {
    dt[, onset_thresh := NA_real_]
    return(dt)
  }
  
  merge(dt, thr_dt, by = "id", all.x = TRUE)
}

# ---------------------------------------------------------------------------
# Cell transform (optional)
# ---------------------------------------------------------------------------

read_cell_transform <- function(spec) {
  enabled <- isTRUE(spec$options$cell_transform_enabled %||% FALSE)
  f <- spec$options$cell_transform_file %||% NULL
  if (!enabled) return(NULL)
  
  if (is.null(f) || !nzchar(f)) stop("cell_transform_enabled=TRUE but options.cell_transform_file is empty.", call. = FALSE)
  if (!file.exists(f)) stop("cell_transform_file not found: ", f, call. = FALSE)
  
  w <- data.table::fread(f)
  req <- c("target_id", "source_id", "weight")
  if (!all(req %in% names(w))) stop("Transform file must have: target_id, source_id, weight", call. = FALSE)
  w[, weight := as.numeric(weight)]
  w[]
}

transform_forecast_wide <- function(dt, weights_dt) {
  if (is.null(weights_dt)) return(dt)
  
  dt <- data.table::copy(dt)
  dt <- add_id_from_latlon(dt)  # ensures id exists if only lat/lon present
  
  day_cols <- grep("_day_\\d+$", names(dt), value = TRUE)
  if (!length(day_cols)) stop("No forecast day columns found to transform (expected *_day_<n>).", call. = FALSE)
  
  # Keep only non-spatial metadata as keys
  drop_cols <- intersect(names(dt), c("lat", "lon", "latitude", "longitude"))
  meta_cols <- setdiff(setdiff(names(dt), c(day_cols, "id")), drop_cols)
  
  long <- data.table::melt(
    dt,
    id.vars = c("id", meta_cols),
    measure.vars = day_cols,
    variable.name = "day_col",
    value.name = "rain"
  )
  
  data.table::setkey(weights_dt, source_id)
  long <- merge(long, weights_dt, by.x = "id", by.y = "source_id", allow.cartesian = TRUE)
  
  # Aggregate across ALL source cells contributing to the same target cell
  trans <- long[, .(rain = sum(weight * rain)),
                by = c(meta_cols, "target_id", "day_col")]
  
  # rename target_id -> id
  trans[, id := target_id]
  trans[, target_id := NULL]
  
  # cast back to wide
  out <- data.table::dcast(
    trans,
    formula = as.formula(paste(paste(meta_cols, collapse = " + "), "+ id ~ day_col")),
    value.var = "rain"
  )
  
  out
}

transform_groundtruth_long <- function(dt, weights_dt, value_col) {
  if (is.null(weights_dt)) return(dt)
  
  dt <- data.table::copy(dt)
  dt <- add_id_from_latlon(dt)
  if (!(value_col %in% names(dt))) stop("Value col not found: ", value_col, call. = FALSE)
  
  drop_cols <- intersect(names(dt), c("lat", "lon", "latitude", "longitude"))
  meta_cols <- setdiff(setdiff(names(dt), c("id", value_col)), drop_cols)
  
  data.table::setkey(weights_dt, source_id)
  x <- merge(dt, weights_dt, by.x = "id", by.y = "source_id", allow.cartesian = TRUE)
  
  x[, (value_col) := as.numeric(get(value_col))]
  
  trans <- x[, .(tmp = sum(weight * get(value_col))),
             by = c(meta_cols, "target_id")]
  
  trans[, id := target_id]
  trans[, target_id := NULL]
  data.table::setnames(trans, "tmp", value_col)
  
  trans[]
}

# ---------------------------------------------------------------------------
# NetCDF readers (stage-1 replacement)
# ---------------------------------------------------------------------------

nc_read_forecast_wide <- function(nc_path, var_name, dim_rename_map, spec,
                                  day_dim = "day", prefix = NULL, add_year = TRUE) {
  prefix <- prefix %||% var_name
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  if (!(var_name %in% names(nc$var))) return(NULL)
  v <- nc$var[[var_name]]
  dims <- v$dim
  if (!length(dims)) return(NULL)
  
dim_vals <- purrr::map(dims, function(d) {
  nm_lc <- tolower(d$name)

  if (nm_lc %in% c("lat", "latitude")) {
    if (!is.null(d$vals) && length(d$vals) > 0) return(as.vector(d$vals))
    return(as.vector(get_nc_coord(nc, "lat")))
  }

  if (nm_lc %in% c("lon", "longitude", "long")) {
    if (!is.null(d$vals) && length(d$vals) > 0) return(as.vector(d$vals))
    return(as.vector(get_nc_coord(nc, "lon")))
  }

  # Handle time coordinate robustly (variable or dimension), similar to lat/lon
  if (nm_lc %in% c("time")) {
    if (!is.null(d$vals) && length(d$vals) > 0) return(as.vector(d$vals))
    return(as.vector(get_nc_time(nc)))
  }

  x <- NULL
  if (!is.null(d$vals) && length(d$vals) > 0) x <- d$vals
  else if (d$name %in% names(nc$var)) x <- ncdf4::ncvar_get(nc, d$name)
  else x <- seq_len(d$len)

  as.vector(x)
})
  dim_names <- purrr::map_chr(dims, ~ .x$name)
  names(dim_vals) <- dim_names
  
  dummy <- tibble::as_tibble(setNames(lapply(dim_vals, function(x) x[1]), dim_names))
  dummy <- rename_dimensions(dummy, dim_rename_map)
  renamed_dim_names <- tolower(names(dummy))
  
  day_idx <- which(renamed_dim_names == tolower(day_dim))
  if (length(day_idx) != 1L) stop("Could not identify day dimension (after rename).", call. = FALSE)
  
  perm <- c(setdiff(seq_along(dims), day_idx), day_idx)
  arr <- ncdf4::ncvar_get(nc, var_name)
  arrp <- aperm(arr, perm)
  
  day_vals <- dim_vals[[ names(dim_vals)[day_idx] ]]
  day_vals_num <- suppressWarnings(as.integer(day_vals))
  
  min_day <- spec$options$min_day %||% NA_integer_
  max_day <- spec$options$max_day %||% NA_integer_
  keep <- rep(TRUE, length(day_vals))
  if (!is.na(min_day)) keep <- keep & !is.na(day_vals_num) & day_vals_num >= min_day
  if (!is.na(max_day)) keep <- keep & !is.na(day_vals_num) & day_vals_num <= max_day
  if (!any(keep)) stop("Day filtering removed all day columns for: ", basename(nc_path), call. = FALSE)
  
  d <- dim(arrp)
  idx <- rep(list(quote(expr = )), length(d))
  idx[[length(d)]] <- which(keep)
  arrp <- do.call(`[`, c(list(arrp), idx, list(drop = FALSE)))
  
  day_vals <- day_vals[keep]
  day_vals_num <- day_vals_num[keep]
  
  other_dim_vals <- dim_vals[ setdiff(seq_along(dim_vals), day_idx) ]
  grid <- do.call(expand.grid, c(other_dim_vals, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)))
  grid <- tibble::as_tibble(grid)
  grid <- rename_dimensions(grid, dim_rename_map)
  grid <- dplyr::rename_with(grid, tolower)
  
  mat <- matrix(as.vector(arrp), nrow = nrow(grid), ncol = length(day_vals), byrow = FALSE)
  day_lab <- ifelse(is.na(day_vals_num), as.character(day_vals), as.character(day_vals_num))
  colnames(mat) <- paste0(prefix, "_day_", day_lab)
  
  out <- data.table::as.data.table(grid)
  out <- cbind(out, data.table::as.data.table(mat))
  
  if (isTRUE(add_year) && ("time" %in% names(out))) {
    if (is.numeric(out$time)) {
      cand <- c(names(nc$var), names(nc$dim))
      hit <- cand[tolower(cand) == "time"][1]
      tu <- ncdf4::ncatt_get(nc, hit, "units")$value
      out[, time := as.Date(nc_time_to_posixct(time, tu))]
    } else {
      out[, time := as.Date(time)]
    }
    out[, year := lubridate::year(time)]
  }
  
  out
}

nc_read_groundtruth_long <- function(nc_path, var_name, dim_rename_map, add_year = TRUE) {
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  if (!(var_name %in% names(nc$var))) return(NULL)
  v <- nc$var[[var_name]]
  dims <- v$dim
  if (!length(dims)) return(NULL)
  
  dim_vals <- purrr::map(dims, function(d) {
    nm_lc <- tolower(d$name)
    
    if (nm_lc %in% c("lat", "latitude")) {
      if (!is.null(d$vals) && length(d$vals) > 0) return(as.vector(d$vals))
      return(as.vector(get_nc_coord(nc, "lat")))
    }
    
    if (nm_lc %in% c("lon", "longitude", "long")) {
      if (!is.null(d$vals) && length(d$vals) > 0) return(as.vector(d$vals))
      return(as.vector(get_nc_coord(nc, "lon")))
    }
    
    # Handle time coordinate robustly (variable or dimension), similar to lat/lon
    if (nm_lc %in% c("time")) {
      if (!is.null(d$vals) && length(d$vals) > 0) return(as.vector(d$vals))
      return(as.vector(get_nc_time(nc)))
    }
    
    x <- NULL
    if (!is.null(d$vals) && length(d$vals) > 0) x <- d$vals
    else if (d$name %in% names(nc$var)) x <- ncdf4::ncvar_get(nc, d$name)
    else x <- seq_len(d$len)
    
    as.vector(x)
  })
  dim_names <- purrr::map_chr(dims, ~ .x$name)
  names(dim_vals) <- dim_names
  
  grid <- do.call(expand.grid, c(dim_vals, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)))
  arr <- ncdf4::ncvar_get(nc, var_name)
  val <- as.vector(arr)
  if (nrow(grid) != length(val)) stop("Dimension grid mismatch for ", var_name, " in ", nc_path, call. = FALSE)
  
  df <- tibble::as_tibble(grid)
  df[[var_name]] <- val
  
  df <- rename_dimensions(df, dim_rename_map)
  df <- dplyr::rename_with(df, tolower)
  
  vcol <- tolower(var_name)
  if (!(vcol %in% names(df)) && (var_name %in% names(df))) names(df)[names(df) == var_name] <- vcol
  df <- dplyr::filter(df, !is.na(.data[[vcol]]))
  
  if (isTRUE(add_year) && ("time" %in% names(df))) {
    if (is.numeric(df$time)) {
      cand <- c(names(nc$var), names(nc$dim))
      hit <- cand[tolower(cand) == "time"][1]
      tu <- ncdf4::ncatt_get(nc, hit, "units")$value
      df <- dplyr::mutate(df, time = nc_time_to_posixct(.data$time, tu))
    }
    df <- dplyr::mutate(df, year = lubridate::year(.data$time))
  }
  
  data.table::as.data.table(df)
}

# ---------------------------------------------------------------------------
# Stage-2 helpers
# ---------------------------------------------------------------------------

order_day_cols <- function(wide_dt, key_cols) {
  cand <- setdiff(names(wide_dt), key_cols)
  day_ints <- suppressWarnings(as.integer(cand))
  ok <- !is.na(day_ints)
  ord <- order(day_ints[ok])
  list(day_cols = cand[ok][ord], day_ints = day_ints[ok][ord])
}

.first_pos_ge_date <- function(cal_dates, cutoff_date) {
  pos <- which(cal_dates >= cutoff_date)[1]
  if (is.na(pos)) 9999 else pos
}


# Vectorized row-wise onset computation
calc_onsets_rowwise <- function(dt, day_cols, day_ints, win) {
  # day_ints MUST be sorted increasing
  stopifnot(isTRUE(all(diff(day_ints) >= 0)))
  
  X <- as.matrix(dt[, ..day_cols])
  storage.mode(X) <- "double"
  
  t0 <- as.Date(dt$time)
  th <- as.numeric(dt$onset_thresh)
  
  # --- Vectorized start positions (FAST) ---
  # Want first j where day_ints[j] >= need
  # j = findInterval(need - 1, day_ints) + 1
  start_from_need <- function(need) {
    j <- findInterval(need - 1L, day_ints) + 1L
    # if need is bigger than max(day_ints), j becomes length(day_ints)+1 -> no valid start
    j[j > length(day_ints)] <- 9999L
    j
  }
  
  # Climatological MOK date cutoff (June 1)
  yr <- as.integer(dt$year)
  june2 <- as.Date(paste0(yr, "-06-02"))
  need_clim_mok_date <- as.integer(june2 - t0)
  start_clim_mok_date <- start_from_need(need_clim_mok_date)
  
  # MOK cutoff (if NA, treat as "no restriction" -> start at 1)
  mk <- as.Date(dt$mok_date)
  start_mok <- rep.int(1L, length(t0))
  ok_mk <- !is.na(mk)
  need_mok <- as.integer(mk[ok_mk] - t0[ok_mk])
  start_mok[ok_mk] <- start_from_need(need_mok)
  
  # --- Per-row onsets (still 3 calls/row, but less overhead around them) ---
  n <- nrow(X)
  onset_raw           <- integer(n)
  onset_clim_mok_date <- integer(n)
  onset_mok           <- integer(n)
  
  for (i in seq_len(n)) {
    s <- X[i, ]
    
    # precompute once
    wsum_all <- .roll_sum_na_rm_left(s, win)
    
    if (length(s) >= 10L) {
      sum10 <- .roll_sum_na_propagate_left(s, 10L)
      bad10 <- !is.na(sum10) & (sum10 < 5)
      pre_bad <- c(0L, cumsum(as.integer(bad10)))
      last10start <- length(s) - 10L + 1L
    } else {
      pre_bad <- 0L
      last10start <- 0L
    }
    
    onset_raw[i]           <- as.integer(find_onset_precomp(s, win, th[i], wsum_all, pre_bad, last10start,
                                                            start_day = 0))
    onset_clim_mok_date[i] <- as.integer(find_onset_precomp(s, win, th[i], wsum_all, pre_bad, last10start,
                                                            start_day = start_clim_mok_date[i]))
    onset_mok[i]           <- as.integer(find_onset_precomp(s, win, th[i], wsum_all, pre_bad, last10start,
                                                            start_day = start_mok[i]))
  }
  
  
  list(onset_raw = onset_raw, onset_clim_mok_date = onset_clim_mok_date, onset_mok = onset_mok)
}

process_rainfall_forecast_id <- function(dt, spec, mok_dt = NULL, thr_dt = NULL) {
  dt <- copy(dt)
  has_number <- "number" %in% names(dt)
  max_number <- spec$filter$max_number %||% NA_integer_
  
  
  dt <- add_id_from_latlon(dt)
  if ("time" %in% names(dt)) dt[, time := as.Date(time)]
  if ("year" %in% names(dt)) dt[, year := as.integer(year)]
  if (has_number) dt[, number := as.integer(number)]
  
  dt <- attach_thresholds_id(dt, thr_dt)
  
  if (!is.null(mok_dt) && "year" %in% names(dt)) {
    setkey(mok_dt, year)
    dt <- merge(dt, mok_dt, by = "year", all.x = TRUE)
  } else {
    dt[, mok_date := as.Date(NA)]
  }
  # Optional filtering by ensemble member number
  if (has_number && !is.na(max_number)) {
    before_n <- nrow(dt)
    dt <- dt[number <= max_number]
    after_n <- nrow(dt)
    
    if (after_n == 0L) {
      stop("After filtering number <= max_number (", max_number, 
           "), no rows remain.", call. = FALSE)
    }
    
    message("Filtered by number <= ", max_number, 
            " : ", before_n, " -> ", after_n, " rows")
  }
  
  
  wide_prefix <- spec$input$wide_prefix %||% tolower(spec$input$value_col)
  wide_pat <- paste0("^", wide_prefix, "_day_\\d+$")
  day_cols_pref <- grep(wide_pat, names(dt), value = TRUE)
  if (!length(day_cols_pref)) stop("No wide day columns matched: ", wide_pat, call. = FALSE)
  
  day_nums <- suppressWarnings(as.integer(sub(paste0("^", wide_prefix, "_day_"), "", day_cols_pref)))
  if (anyNA(day_nums)) stop("Could not parse day numbers from wide day columns.", call. = FALSE)
  setnames(dt, day_cols_pref, as.character(day_nums))
  
  key_base <- c("id", "time")
  if ("year" %in% names(dt)) key_base <- c(key_base, "year")
  key_base <- c(key_base, "onset_thresh", "mok_date")
  
  key_member <- key_base
  if (has_number) key_member <- c(key_member, "number")
  
  day_info <- order_day_cols(dt, key_member)
  day_cols <- day_info$day_cols
  day_ints <- day_info$day_ints
  
  win <- as.integer(spec$options$window)[1]
  min_day <- as.integer(spec$options$min_day)
  max_day <- as.integer(spec$options$max_day)
  
  keep <- day_ints >= min_day & day_ints <= (max_day + win - 1L)
  day_cols <- day_cols[keep]
  day_ints <- day_ints[keep]
  if (!length(day_cols)) stop("After min_day/max_day filtering, no day columns remain.", call. = FALSE)
  

  # --- Per-row onsets ---
  on <- calc_onsets_rowwise(dt, day_cols = day_cols, day_ints = day_ints, win = win)
  dt[, onset_raw           := on$onset_raw]
  dt[, onset_clim_mok_date := on$onset_clim_mok_date]
  dt[, onset_mok           := on$onset_mok]
  
  D <- length(day_ints)
  
  prob_from_idx <- function(idx, D) {
    denom <- length(idx)
    idx <- idx[!is.na(idx) & idx >= 1L & idx <= D]
    if (!length(idx)) return(rep(0, D))
    tabulate(idx, nbins = D) / denom
  }
  
  # --- stats (mean/sd/frac) in ONE pass ---
  stats_dt <- dt[, {
    m  <- lapply(.SD, mean)
    s  <- lapply(.SD, sd)
    fr <- lapply(.SD, function(x) mean(x > 1))
    c(m, s, fr)
  }, by = key_base, .SDcols = day_cols]
  
  K <- length(day_cols)
  
  # indices of the generated columns (after key_base)
  cols_start <- length(key_base) + 1L
  idx_mean <- cols_start:(cols_start + K - 1L)
  idx_sd   <- (cols_start + K):(cols_start + 2*K - 1L)
  idx_frac <- (cols_start + 2*K):(cols_start + 3*K - 1L)
  
  setnames(stats_dt, idx_mean, paste0("forecast_rain_day_",    day_ints))
  setnames(stats_dt, idx_sd,   paste0("forecast_rain_sd_day_", day_ints))
  setnames(stats_dt, idx_frac, paste0("frac_raining_day_",     day_ints))
  
  # --- probs (raw/clim_mok_date/mok) in ONE pass ---
  prob_dt <- dt[, {
    pr <- prob_from_idx(onset_raw,           D)
    pj <- prob_from_idx(onset_clim_mok_date, D)
    pm <- prob_from_idx(onset_mok,           D)
    as.list(c(pr, pj, pm))
  }, by = key_base]

  setnames(
    prob_dt,
    old = names(prob_dt)[-(1:length(key_base))],
    new = c(paste0("predicted_prob_day_",               day_ints),
            paste0("predicted_prob_clim_mok_date_day_", day_ints),
            paste0("predicted_prob_mok_day_",           day_ints))
  )
  
  setkeyv(stats_dt, key_base)
  setkeyv(prob_dt,  key_base)
  
  list(wide = stats_dt[prob_dt])
}

process_ground_truth_rainfall_id <- function(dt, spec, mok_dt = NULL, thr_dt = NULL, value_col) {
  dt <- copy(dt)
  dt <- add_id_from_latlon(dt)
  dt[, time := as.Date(time)]
  dt[, year := as.integer(year)]
  dt[, (value_col) := as.numeric(get(value_col))]
  
  dt <- attach_thresholds_id(dt, thr_dt)
  
  if (!is.null(mok_dt)) {
    setkey(mok_dt, year)
    dt <- merge(dt, mok_dt, by = "year", all.x = TRUE)
  } else {
    dt[, mok_date := as.Date(NA)]
  }
  
  cutoff_md <- spec$options$cutoff_month_day
  dt[, cutoff_date := as.Date(paste0(year, "-", cutoff_md))]
  dt[, start_date := cutoff_date]
  dt[!is.na(mok_date), start_date := pmax(cutoff_date, mok_date)]
  
  keys <- c("id", "year")
  setorderv(dt, c(keys, "time"))
  
  win <- as.integer(spec$options$window)[1]
  
  onset_dt <- dt[, {
    series <- get(value_col)
    th <- onset_thresh[1]
    dates <- time
    sd <- start_date[1]
    
    # start position (avoid Inf -> as.integer(Inf) -> NA crash)
    start_pos <- {
      p <- which(dates >= sd)[1]
      if (is.na(p)) 9999 else as.integer(p)
    }
    if (all(is.na(series)) || length(series) < win || is.na(th)) {
      list(mr_onset_idx = NA_real_, mr_onset_date = as.Date(NA))
    } else {
      # --- precompute once per group ---
      wsum_all <- .roll_sum_na_rm_left(series, win)
      
      if (length(series) >= 10L) {
        sum10 <- .roll_sum_na_propagate_left(series, 10L)
        bad10 <- !is.na(sum10) & (sum10 < 5)
        pre_bad <- c(0L, cumsum(as.integer(bad10)))
        last10start <- length(series) - 10L + 1L
      } else {
        # not used when n<10 (function returns early), but keep valid shapes
        pre_bad <- 0L
        last10start <- 0L
      }
      
      mr_idx <- find_onset_precomp(
        series = series,
        win = win,
        thresh = th,
        wsum_all = wsum_all,
        pre_bad = pre_bad,
        last10start = last10start,
        reject_if_short_followup = TRUE,
        start_day = start_pos
      )
      
      list(
        mr_onset_idx = as.numeric(mr_idx),
        mr_onset_date = if (!is.na(mr_idx)) as.Date(dates[mr_idx]) else as.Date(NA)
      )
    }
  }, by = keys]
  
  onset_dt[, cutoff_date := as.Date(paste0(year, "-", cutoff_md))]
  onset_dt[, mr_onset_day := as.numeric(mr_onset_date - cutoff_date)]
  
  annotated <- merge(dt, onset_dt, by = keys, all.x = TRUE)
  annotated[, mr_onset_flag := (!is.na(mr_onset_date) & time == mr_onset_date)]
  
  list(onset_dates = onset_dt, annotated = annotated)
}

run_single_pipeline <- function(spec_id) {
  spec <- load_spec(spec_id, type = "raw_data")
  spec <- validate_spec_single(spec)
  
  # parallel plan
  if (isTRUE(spec$input$parallel)) {
    workers <- spec$input$workers %||% 25
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::sequential)
  }
  
  files_dt <- list_nc_files_with_year(spec)
  years <- sort(unique(files_dt$year))
  
  dim_rename_map <- spec$dimensions$rename %||% list()
  var_name <- get_value_var(spec)
  
  # aux tables
  mok_dt <- NULL
  if (!is.null(spec$mok$file) && nzchar(spec$mok$file)) {
    mok_dt <- read_mok_dates(spec)  # Only if MOK filtering is desired
  }
  thr_dt <- read_thresholds(spec)       # From onset_utils.R
  thr_dt <- prep_thresholds_id(thr_dt)  # convert to id-keyed thresholds
  
  weights_dt <- read_cell_transform(spec)
  
  out_dir <- spec$output$out_dir
  out_base <- file.path(out_dir, spec_id)
  
  all_wide <- list()
  all_long <- list()
  
  map_fn <- if (isTRUE(spec$input$parallel)) furrr::future_map else purrr::map
  
  for (yr in years) {
    message("=== Processing year: ", yr, " ===")
    fpaths <- files_dt[year == yr, nc_path]
    
    if (spec$type == "rainfall_forecast") {
      wide_prefix <- spec$input$wide_prefix %||% tolower(spec$input$value_col)
      day_dim <- spec$input$wide_day_dim %||% "day"
      
      parts <- map_fn(fpaths, function(fp) {
        nc_read_forecast_wide(
          nc_path = fp,
          var_name = var_name,
          dim_rename_map = dim_rename_map,
          spec = spec,
          day_dim = day_dim,
          prefix = wide_prefix,
          add_year = TRUE
        )
      })
      
      dt <- data.table::rbindlist(purrr::compact(parts), fill = TRUE)
      if (!nrow(dt)) next
      
      # enforce year from filename (your requirement)
      dt[, year := as.integer(yr)]
      
      dt <- add_id_from_latlon(dt)
      
      # optional transform BEFORE onset detection
      dt <- transform_forecast_wide(dt, weights_dt)
      
      res <- process_rainfall_forecast_id(dt, spec, mok_dt = mok_dt, thr_dt = thr_dt)
      
      all_wide[[length(all_wide) + 1L]] <- res$wide
      
    } else if (spec$type == "ground_truth_rainfall") {
      
      parts <- map_fn(fpaths, function(fp) {
        nc_read_groundtruth_long(
          nc_path = fp,
          var_name = var_name,
          dim_rename_map = dim_rename_map,
          add_year = TRUE
        )
      })
      
      dt <- data.table::rbindlist(purrr::compact(parts), fill = TRUE)
      if (!nrow(dt)) next
      
      dt[, year := as.integer(yr)]
      dt <- add_id_from_latlon(dt)
      
      val_col <- tolower(var_name)
      if (!(val_col %in% names(dt)) && (var_name %in% names(dt))) {
        data.table::setnames(dt, var_name, val_col)
      }
      
      dt <- transform_groundtruth_long(dt, weights_dt, value_col = val_col)
      
      res <- process_ground_truth_rainfall_id(dt, spec, mok_dt = mok_dt, thr_dt = thr_dt, value_col = val_col)
      
      all_wide[[length(all_wide) + 1L]] <- res$onset_dates
      all_long[[length(all_long) + 1L]] <- res$annotated
      
    } else {
      stop("Unknown spec$type: ", spec$type, call. = FALSE)
    }
  }

  # Write accumulated results as RDS
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(all_wide) > 0L) {
    wide_path <- paste0(out_base, "_wide.rds")
    saveRDS(data.table::rbindlist(all_wide, fill = TRUE), wide_path)
    message("Wrote wide RDS: ", wide_path)
  }
  if (length(all_long) > 0L) {
    long_path <- paste0(out_base, "_long.rds")
    saveRDS(data.table::rbindlist(all_long, fill = TRUE), long_path)
    message("Wrote long RDS: ", long_path)
  }

  invisible(TRUE)
}

