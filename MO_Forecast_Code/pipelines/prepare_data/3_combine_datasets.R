# ==============================================================================
# Script: 3_combine_datasets.R
#
# Purpose:
#   Combine climatology + forecast families  + ground truth into one combined dataset for modeling
#   / downstream processing.
#
# Inputs:
#   - Ground truth: (lat, lon, year, true_onset_day, true_onset_date)
#       read via read_ground_truth_wide()
#   - Climatology: (lat, lon, time, year, clim_p_onset_day_<k>, [clim plus bin])
#       read via read_and_format_climatology_wide()
#   - Forecast families: produced by format_forecast_family() for each entry in spec$forecasts
#       * daily:      (lat, lon, time, year, <family>_<out>_day_<k>[plus])
#       * constants:  (lat, lon, time, year, <family>_<out>)
#
# What the script does:
#   1) Loads a "combine" spec by --spec_id (YAML-driven).
#   2) Reads ground truth, climatology WIDE, and each forecast family.
#   3) Merges all DAILY WIDE tables on (lat, lon, time, year).
#      - If spec$options$join == "full": uses full outer join (union of keys)
#      - Else: uses inner-style behavior via merge(all = FALSE) in the Reduce
#   4) Ensures `year` reflects issue-date year (year(time)), then merges ground truth
#      on (lat, lon, year).
#   5) Optional trimming:
#      - If spec$options$trim_forecasts_after_true_onset is TRUE, drops rows where
#        time > true_onset_date (keeping NA true_onset_date rows).
#   6) Merges all constants WIDE tables (full outer among constants) and then left-joins
#      them into the combined daily table on (lat, lon, time, year).
#   7) Writes the combined WIDE CSV to:
#        file.path(spec$output$out_dir, paste0(spec$id, "_combined_wide.csv"))
#
# Outputs:
#   - One CSV containing keys (lat, lon, time, year), climatology day columns,
#     all forecast daily day columns (and plus bins where configured), forecast
#     constants, and ground-truth onset columns.
#
# Dependencies:
#   - optparse, yaml, data.table, lubridate
#   - project helpers:
#       R/_shared/misc.R
#       R/_shared/read_spec.R
#       R/prepare_data/combine_forecasts_utils.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(data.table)
  library(lubridate)
})

source("R/_shared/misc.R")
source("R/_shared/read_spec.R")
source("R/prepare_data/combine_forecasts_utils.R")

# ----------------------------
# CLI
# ----------------------------
if(interactive())
{
  spec_id = "combine_template_clim_mok_date_2025"  
}else{
  option_list <- list(
    make_option(c("-i", "--spec_id"), type = "character", default = NULL)
  )
  
  opt <- parse_args(OptionParser(option_list = option_list))
  spec_id = opt$spec_id
  if (is.null(opt$spec_id))
    stop("Provide --spec_id", call. = FALSE)
}
spec <- load_spec(spec_id, "combine")


# ----------------------------
# Read inputs 
# ----------------------------

truth <- read_ground_truth_wide(spec$input$ground_truth_wide_rds)

clim_list <- list()

if (!is.null(spec$input$climatologies) && length(spec$input$climatologies) > 0L) {
  # named list/map preferred
  clim_list <- lapply(names(spec$input$climatologies), function(nm) {
    x <- spec$input$climatologies[[nm]]
    read_and_format_climatology_wide(
      path       = x$rds,
      out_prefix = x$out_prefix %||% paste0(nm, "_p_onset")
    )
  })
  names(clim_list) <- names(spec$input$climatologies)
  
} else {
  # backward compatible single-climatology config
  clim_list <- list(
    read_and_format_climatology_wide(
      path       = spec$input$climatology_rds,
      out_prefix = "clim_p_onset"
    )
  )
}
# ---------------------------------------------------------------------------

forecast_parts <- lapply(
  names(spec$forecasts),
  function(nm) format_forecast_family(nm, spec$forecasts[[nm]])
)
names(forecast_parts) <- names(spec$forecasts)

forecast_daily <- lapply(forecast_parts, `[[`, "daily")
forecast_const <- lapply(forecast_parts, `[[`, "constants")

# combine DAILY: now c(clim_list, forecast_daily)
daily_wide <- Reduce(
  function(a, b)
    merge(
      a, b,
      by = c("id", "time", "year"),
      all = (spec$options$join == "full")
    ),
  c(clim_list, forecast_daily)
)

# ----------------------------
# Exclude forecasts made after true onset 
# ----------------------------

daily_wide[, year := year(time)]

daily_wide <- merge(
  daily_wide,
  truth,
  by = c("id", "year"),
  all.x = TRUE
)

if (isTRUE(spec$options$trim_forecasts_after_true_onset)) {
  daily_wide <- daily_wide[
    is.na(true_onset_date) | time <= true_onset_date
  ]
}

# ----------------------------
# Merge constants 
# ----------------------------

if (length(forecast_const)) {
  const_all <- Reduce(
    function(a, b)
      merge(a, b, by = c("id", "time", "year"), all = TRUE),
    forecast_const
  )
  
  daily_wide <- merge(
    daily_wide,
    const_all,
    by = c("id", "time", "year"),
    all.x = TRUE
  )
}

# ----------------------------
# Write output
# ----------------------------

out_rds <- file.path(spec$output$out_dir,
                     paste0(spec_id, "_combined_wide.rds"))

saveRDS(daily_wide, out_rds)
message("Wrote combined wide dataset: ", out_rds)
