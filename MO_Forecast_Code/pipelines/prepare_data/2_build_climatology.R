# ==============================================================================
# Script: 2_build_climatology.R
#
# Purpose:
#   Fit per-cell (grid cell id) onset-day KDEs from ground truth (training years)
#   and produce issue-date probability forecasts over a within-season date window
#   (test/output years).
#
#   This script supports MULTIPLE climatology runs from a single YAML spec:
#     - spec$climatologies is a named list of run configurations
#     - each run can use different training/test year ranges, issue windows,
#       onset column, and horizon rules
#     - outputs are written once per run (CSV + RDS), keyed by out_stem
#
# Description:
#   Consumes the ground-truth output written by 1_process_raw_nc_files.R
#   (spec$type == "ground_truth_rainfall"), specifically:
#
#     <spec$output$out_dir>/<spec$id>_wide.csv
#
#   For each run in spec$climatologies and for each grid cell id, fits a KDE to
#   historical onset_day values in the run's training years. When "conditional" is set to False, 
#   this produces a static climatology forecast. When "conditional" is set to True, this produces
#   the Evolving-Expectations  conditional model of climatology. In other words:
#   For each issue date in the run's within-season window, the script produces conditional probabilities:
#
#     P(onset occurs on day d0(time)+k | onset has not occurred by d0(time)),
#       for k = 1..H(time)
#
#   where d0(time) is the day offset from the run's season_start_md and H(time)
#   is either a fixed forecast_window or a piecewise horizon schedule.
#
# Inputs:
#   - YAML spec (specs/<type>/<id>.yml) containing:
#       - output.out_dir (basename derived from spec$id)
#       - paths.climatology_out_dir (optional)
#       - climatologies: named map of run configs, each with:
#           train_year_min/max, test_year_min/max, season_start_md, issue_end_md,
#           onset_col (optional), and either forecast_window or horizons
#   - Stage-2 *_wide.csv with columns: id, year, <onset_col> (per run)
#
# Outputs (per run):
#   - <paths$climatology_out_dir>/<out_stem>.csv : wide table of issue-date
#     probabilities with columns: time, year, id, model, predicted_prob_day_1..N
#   - <paths$climatology_out_dir>/<out_stem>.rds : list(spec, run, options,
#     forecasts, kdes)
#
# Notes:
#   - Training and test years may overlap. In this case forecasts will be produced for the test years by cross-validation (ie training on all other years)
#   - Forecast probabilities do not necessarily sum to 1 over the horizon, since
#     the horizon may not cover the full remaining tail of the onset distribution. The combine data script will add a "29plus" 
#     day to account for this
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(here)
})

source("R/_shared/misc.R")
source("R/_shared/read_spec.R")
source("R/prepare_data/climatology_utils.R")

option_list <- list(
  make_option(c("-i", "--spec_id"), type = "character", default = NULL,
              help = "spec id (loads specs/<type>/<id>.yml)", metavar = "character"),
  make_option(c("-r", "--run"), type = "character", default = NULL,
              help = "Optional: run only one climatology key (e.g. 'clim_pre_2000')", metavar = "character")
)

if (interactive()) {
  opt_cli <- list(spec_id = "imd",  run = NULL)
} else {
  opt_cli <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt_cli$spec_id)) stop("Please provide --spec_id <id>.", call. = FALSE)
}

spec <- load_spec(opt_cli$spec_id, type = "raw_data")
spec$id <- opt_cli$spec_id
# Validate we have multiple climatology runs
spec <- validate_spec(
  spec,
  required_paths = c("output.out_dir", "climatologies"),
  checks = list(function(s) {
    if (is.null(s$climatologies) || length(s$climatologies) == 0L)
      stop("Missing or empty 'climatologies' in spec.", call. = FALSE)
    
    # basic validation per run
    required <- c("train_year_min","train_year_max","test_year_min","test_year_max",
                  "season_start_md","issue_end_md")
    for (nm in names(s$climatologies)) {
      co <- s$climatologies[[nm]]
      miss <- required[!required %in% names(co)]
      if (length(miss) > 0L) {
        stop("Run '", nm, "' missing fields: ", paste(miss, collapse = ", "), call. = FALSE)
      }
      has_fw <- !is.null(co$forecast_window)
      has_hz <- !is.null(co$horizons)
      if (!has_fw && !has_hz)
        stop("Run '", nm, "': provide forecast_window or horizons.", call. = FALSE)
      if (has_fw && has_hz)
        stop("Run '", nm, "': provide only one of forecast_window or horizons.", call. = FALSE)
    }
  })
)

# Shared paths (GT path + output dir). Output stem becomes per-run.
paths <- get_paths_clim(spec)

# Read ground truth once
# (onset_col may differ per run, so we read the full wide file once and then pick columns per run)
gt_wide_path <- paths$gt_path
if (!file.exists(gt_wide_path)) stop("Ground-truth file not found: ", gt_wide_path, call. = FALSE)
gt_wide <- readRDS(gt_wide_path)

# Optionally subset which runs to execute
run_keys <- names(spec$climatologies)
if (!is.null(opt_cli$run)) {
  if (!opt_cli$run %in% run_keys) {
    stop("Requested --run '", opt_cli$run, "' not found. Available: ",
         paste(run_keys, collapse = ", "), call. = FALSE)
  }
  run_keys <- opt_cli$run
}

# Run loop
results <- purrr::imap(run_keys, function(run_key, idx_unused) {
  co <- spec$climatologies[[run_key]]
  opt <- get_climatology_options_from_run(co)
  
  # Read GT onset for this run's onset_col (from the already-loaded wide table)
  gt <- read_gt_onset_from_tbl(gt_wide, onset_col = opt$onset_col)
  
  gt_train <- filter_gt_training(gt, opt$train_year_min, opt$train_year_max)
  
  issue_grid <- build_issue_grid(
    opt$test_year_min,
    opt$test_year_max,
    opt$season_start_md,
    opt$issue_end_md
  )
  
  out <- compute_all_forecasts(
    gt_train        = gt_train,
    issue_grid      = issue_grid,
    season_start_md = opt$season_start_md,
    forecast_window = opt$forecast_window,
    horizons        = opt$horizons,
    conditional     = opt$conditional,
    cv_by_year      = opt$cv_by_year
  )
  
  forecast_tbl <- out$forecasts
  
  out_stem <- co$out_stem %||% paste0(paths$out_stem, "_", run_key)
  out_rds <- file.path(paths$out_dir, paste0(out_stem, ".rds"))

  saveRDS(forecast_tbl, out_rds)

  message("[", run_key, "] Wrote RDS: ", out_rds)

  list(run = run_key, out_rds = out_rds)
})

