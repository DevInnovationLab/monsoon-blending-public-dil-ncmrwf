# ==============================================================================
# File: regrid_utils.R
# ==============================================================================
# Purpose
#   Utility functions for building regridding weights between two rectilinear
#   grids. Used by 0_build_weights.R.
#
#   Two methods are supported:
#     - "conservative": area-weighted, correct for fine -> coarse regridding.
#       weight = overlap_area / target_area, normalized to sum to 1.
#       matches xesmf conservative regridding methodology.
#     - "equal_overlap": equal weights across all overlapping source cells,
#       appropriate for coarse -> fine disaggregation.
#
# Functions
#   make_grid_cells(lat_vec, lon_vec, lat_step, lon_step)
#     builds a data.table of cell bounds and ids from grid center coordinates.
#
#   calc_overlap_area(...)
#     vectorized area of overlap between two sets of rectangles on a sphere,
#     using cos(lat) scaling.
#
#   build_conservative_weights(source_nc_file, target_nc_file, ...)
#     builds area-weighted conservative regridding weights (fine -> coarse).
#
#   build_equal_overlap_weights(source_nc_file, target_nc_file, ...)
#     builds equal-overlap regridding weights (coarse -> fine).
#
#   build_weights_from_spec(spec)
#     dispatches to the appropriate function based on spec$method.
#     called by 0_build_weights.R.
#
# Dependencies
#   - ncdf4, data.table
#   - get_nc_coord() from R/prepare_data/nc_utils.R (must be sourced before this)
# ==============================================================================

# ---- helper: pick a representative nc file from a folder ----
get_representative_nc <- function(nc_folder) {
  # pick the first available nc file in the folder (alphabetical order)
  files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
  if (length(files) == 0L)
    stop("no .nc files found in: ", nc_folder, call. = FALSE)
  files[[1]]
}

# ---- helper: build cell bounds for rectilinear grid ----
make_grid_cells <- function(lat_vec, lon_vec, lat_step, lon_step) {
  lat <- as.numeric(lat_vec)
  lon <- as.numeric(lon_vec)
  
  lat_half <- lat_step / 2
  lon_half <- lon_step / 2
  
  # cartesian product of centers
  g <- CJ(lat = lat, lon = lon, unique = TRUE)
  
  g[, `:=`(
    lat_min = lat - lat_half,
    lat_max = lat + lat_half,
    lon_min = lon - lon_half,
    lon_max = lon + lon_half
  )]
  
  # id as latitude_longitude
  g[, id := paste0(
    format(round(lat, 2), scientific = FALSE, trim = TRUE, nsmall = 2),
    "_",
    format(round(lon, 2), scientific = FALSE, trim = TRUE, nsmall = 2)
  )]
  
  g[]
}

# ---- calculate overlap area between two rectangles on sphere ----
# for conservative regridding: weight = overlap_area / target_area
calc_overlap_area <- function(a_lat_min, a_lat_max, a_lon_min, a_lon_max,
                              b_lat_min, b_lat_max, b_lon_min, b_lon_max) {
  overlap_lat_min <- pmax(a_lat_min, b_lat_min)
  overlap_lat_max <- pmin(a_lat_max, b_lat_max)
  overlap_lon_min <- pmax(a_lon_min, b_lon_min)
  overlap_lon_max <- pmin(a_lon_max, b_lon_max)
  
  has_overlap <- (overlap_lat_min < overlap_lat_max) &
    (overlap_lon_min < overlap_lon_max)
  
  lat_extent <- pmax(0, overlap_lat_max - overlap_lat_min)
  lon_extent <- pmax(0, overlap_lon_max - overlap_lon_min)
  
  # latitude scaling for area on sphere: area ~ cos(latitude)
  lat_center <- (overlap_lat_min + overlap_lat_max) / 2
  area <- cos(lat_center * pi / 180) * lat_extent * lon_extent
  
  ifelse(has_overlap, area, 0)
}

# ---- build conservative regridding weights (fine -> coarse) ----
build_conservative_weights <- function(source_nc_file, target_nc_file,
                                       source_lat_step, source_lon_step,
                                       target_lat_step, target_lon_step,
                                       out_csv) {
  if (!(source_lat_step < target_lat_step && source_lon_step < target_lon_step))
    warning("expected source grid to be finer than target grid")
  
  nc_s <- nc_open(source_nc_file); on.exit(nc_close(nc_s), add = TRUE)
  nc_t <- nc_open(target_nc_file); on.exit(nc_close(nc_t), add = TRUE)
  
  src <- make_grid_cells(
    get_nc_coord(nc_s, "lat"), get_nc_coord(nc_s, "lon"),
    source_lat_step, source_lon_step
  )
  tgt <- make_grid_cells(
    get_nc_coord(nc_t, "lat"), get_nc_coord(nc_t, "lon"),
    target_lat_step, target_lon_step
  )
  
  cat("source grid:", nrow(src), "cells\n")
  cat("target grid:", nrow(tgt), "cells\n")
  
  weights_list <- vector("list", nrow(tgt))
  
  for (i in seq_len(nrow(tgt))) {
    ti <- tgt[i]
    
    cand <- src[
      lat_max > ti$lat_min & lat_min < ti$lat_max &
        lon_max > ti$lon_min & lon_min < ti$lon_max
    ]
    
    if (nrow(cand) == 0L) {
      weights_list[[i]] <- NULL
      next
    }
    
    overlap_areas <- calc_overlap_area(
      cand$lat_min, cand$lat_max, cand$lon_min, cand$lon_max,
      ti$lat_min,   ti$lat_max,   ti$lon_min,   ti$lon_max
    )
    
    target_area <- cos(ti$lat * pi / 180) * target_lat_step * target_lon_step
    weights     <- overlap_areas / target_area
    weights     <- weights / sum(weights)  # normalize to handle edge cases
    
    weights_list[[i]] <- data.table(
      target_id = ti$id,
      source_id = cand$id,
      weight    = weights
    )
  }
  
  w <- rbindlist(weights_list, fill = TRUE)
  if (!nrow(w))
    stop("no overlaps found between target and source grids", call. = FALSE)
  
  # verify weights sum to ~1 for each target cell
  weight_sums <- w[, .(weight_sum = sum(weight)), by = target_id]
  if (any(abs(weight_sums$weight_sum - 1) > 0.01))
    warning("some target cells have weights that don't sum to 1")
  
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  fwrite(w, out_csv)
  cat("wrote weights:", out_csv, "\n")
  cat("  rows:", nrow(w), "\n")
  cat("  unique targets:", length(unique(w$target_id)), "\n")
  cat("  unique sources:", length(unique(w$source_id)), "\n")
  
  invisible(w)
}

# ---- build equal-overlap regridding weights (coarse -> fine) ----
build_equal_overlap_weights <- function(source_nc_file, target_nc_file,
                                        source_lat_step, source_lon_step,
                                        target_lat_step, target_lon_step,
                                        out_csv) {
  if (!(source_lat_step > target_lat_step && source_lon_step > target_lon_step))
    warning("expected source grid cells to be larger than target grid cells")
  
  nc_s <- nc_open(source_nc_file); on.exit(nc_close(nc_s), add = TRUE)
  nc_t <- nc_open(target_nc_file); on.exit(nc_close(nc_t), add = TRUE)
  
  src <- make_grid_cells(
    get_nc_coord(nc_s, "lat"), get_nc_coord(nc_s, "lon"),
    source_lat_step, source_lon_step
  )
  tgt <- make_grid_cells(
    get_nc_coord(nc_t, "lat"), get_nc_coord(nc_t, "lon"),
    target_lat_step, target_lon_step
  )
  
  cat("source grid:", nrow(src), "cells\n")
  cat("target grid:", nrow(tgt), "cells\n")
  
  setkey(src, lat_min, lat_max)
  weights_list <- vector("list", nrow(tgt))
  
  for (i in seq_len(nrow(tgt))) {
    ti <- tgt[i]
    
    cand <- src[
      lat_max > ti$lat_min & lat_min < ti$lat_max &
        lon_max > ti$lon_min & lon_min < ti$lon_max
    ]
    
    if (nrow(cand) == 0L) {
      weights_list[[i]] <- NULL
      next
    }
    
    k <- nrow(cand)
    weights_list[[i]] <- data.table(
      target_id = ti$id,
      source_id = cand$id,
      weight    = rep(1 / k, k)
    )
  }
  
  w <- rbindlist(weights_list, fill = TRUE)
  if (!nrow(w))
    stop("no overlaps found between target and source grids", call. = FALSE)
  
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  fwrite(w, out_csv)
  cat("wrote weights:", out_csv, "\n")
  cat("  rows:", nrow(w), "\n")
  cat("  unique targets:", length(unique(w$target_id)), "\n")
  cat("  unique sources:", length(unique(w$source_id)), "\n")
  
  invisible(w)
}

# ---- dispatch: build weights from a loaded spec ----
build_weights_from_spec <- function(spec) {
  out_csv <- spec$out_csv
  stopifnot(!is.null(out_csv))
  
  # skip if weights file already exists
  if (file.exists(out_csv)) {
    cat("weights file already exists, skipping:", out_csv, "\n")
    return(invisible(NULL))
  }
  
  # pick representative nc files from folders
  source_nc_file <- get_representative_nc(spec$source$nc_folder)
  target_nc_file <- get_representative_nc(spec$target$nc_folder)
  
  cat("source nc:", source_nc_file, "\n")
  cat("target nc:", target_nc_file, "\n")
  
  method <- spec$method
  if (is.null(method))
    stop("spec must specify 'method': 'conservative' or 'equal_overlap'", call. = FALSE)
  
  if (method == "conservative") {
    build_conservative_weights(
      source_nc_file  = source_nc_file,
      target_nc_file  = target_nc_file,
      source_lat_step = spec$source$lat_step,
      source_lon_step = spec$source$lon_step,
      target_lat_step = spec$target$lat_step,
      target_lon_step = spec$target$lon_step,
      out_csv         = out_csv
    )
  } else if (method == "equal_overlap") {
    build_equal_overlap_weights(
      source_nc_file  = source_nc_file,
      target_nc_file  = target_nc_file,
      source_lat_step = spec$source$lat_step,
      source_lon_step = spec$source$lon_step,
      target_lat_step = spec$target$lat_step,
      target_lon_step = spec$target$lon_step,
      out_csv         = out_csv
    )
  } else {
    stop("unknown method: '", method, "'. use 'conservative' or 'equal_overlap'", call. = FALSE)
  }
}