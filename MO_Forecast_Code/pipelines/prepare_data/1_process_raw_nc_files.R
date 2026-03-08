# ==============================================================================
# File: 1_process_raw_nc_files.R
# ==============================================================================
# Purpose
#   Script to process *raw NetCDF rainfall files* directly into
#   monsoon-onset summary outputs.
#
#   This script is the driver/wrapper around run_single_pipeline() (defined in
#   R/prepare_data/nc_utils.R). It:
#     1) Reads a YAML spec (specs/raw_data/<spec_id>.yml)
#     2) Configures parallelism (if enabled in the spec)
#     3) Iterates over NetCDF files year-by-year (year extracted from filenames)
#     4) Produces final outputs (appended across years) in the output directory
#
# Pipeline modes (controlled by spec$type)
#   - "rainfall_forecast"
#       * Reads forecast NetCDFs using nc_read_forecast_wide()
#       * Optionally applies a grid-cell linear-combination transform (weights)
#       * Computes onset indices per ensemble member and aggregates summaries
#       * Writes/append: <output.out_dir>/<output.basename>_wide.csv
#
#   - "ground_truth_rainfall"
#       * Reads observed NetCDFs using nc_read_groundtruth_long()
#       * Optionally applies the same grid-cell transform (weights)
#       * Computes deterministic onset date per (id, year) and annotates series
#       * Writes/appends:
#           - <output.out_dir>/<output.basename>_wide.csv  (onset table)
#           - <output.out_dir>/<output.basename>_long.csv  (annotated daily series)
#
# Inputs
#   - YAML spec:
#       specs/raw_data/<spec_id>.yml
#     Key sections used by downstream utilities:
#       * type
#       * input: nc_folder, file_regex, value_col, wide_prefix, wide_day_dim,
#                parallel, workers, etc.
#       * dimensions: rename map for NetCDF dimension/coord names -> standardized
#       * options: window, (min_day/max_day for forecasts), cutoff_month_day for
#                  ground truth, optional min_year/max_year
#       * thresholds / mok / optional cell transform settings
#
# Outputs
#   - Always: combined (all years) CSV(s) written to output$out_dir with stem
#     output$basename.
#     Exact filenames depend on spec$type:
#       * rainfall_forecast      -> <stem>_wide.csv
#       * ground_truth_rainfall  -> <stem>_wide.csv and <stem>_long.csv
#
# Parallelism
#   - Controlled by spec$input$parallel and spec$input$workers.
#   - The underlying run_single_pipeline() uses future/furrr when enabled.
#   - Year-by-year iteration is sequential (for memory safety); files *within a
#     year* may be processed in parallel depending on your spec.
#
# How to run
#   From an sh file, see "submit" folder for details
#   In RStudio:
#     Run the script (change spec_id in the interactive() part of the code below)
#
# Dependencies
#   - R/_shared/misc.R        (%||% helper, etc.)
#   - R/_shared/read_spec.R   (load_spec(), validate_spec())
#   - R/prepare_data/onset_utils.R (find_onset(), read_mok_dates(), read_thresholds())
#   - R/prepare_data/nc_utils.R    (run_single_pipeline() and supporting fns)
#
# Notes
#   - Year is extracted from filenames via regex \\b(19|20)\\d{2}\\b, so filenames
#     like "2000.nc" and "data_2000.nc" both work.
#   - Spatial grouping is by `id` (lat_lon), not by separate lat/lon keys.
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(future)
  library(furrr)
  library(tidyr)
  library(ncdf4)
  library(data.table)
  library(lubridate)
})

source("R/_shared/misc.R")
source("R/_shared/read_spec.R")
source("R/prepare_data/onset_utils.R") # find_onset(), etc.
source("R/prepare_data/nc_utils.R") # main functions for this script




# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--spec_id"), type = "character", default = NULL,
              help = "spec id (loads specs/raw_data/<id>.yml)", metavar = "character")
)

if (sys.nframe() == 0) {
  if (interactive()) {
    opt <- list(spec_id = "aifs")
  } else {
    opt <- parse_args(OptionParser(option_list = option_list))
    if (is.null(opt$spec_id)) stop("Please provide --spec_id <id>.", call. = FALSE)
    }
  
  run_single_pipeline(opt$spec_id)
}
