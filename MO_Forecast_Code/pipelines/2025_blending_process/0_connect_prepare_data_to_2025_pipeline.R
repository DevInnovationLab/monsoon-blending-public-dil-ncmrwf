# ==============================================================================
# File: 0_connect_prepare_data_to_2025_pipeline.R
# ==============================================================================
# Purpose
#   Convert day-level wide RDS from prepare_data pipeline into a weekly-bin
#   RDS that 1_blend_evaluation.R expects. Aggregates daily onset probabilities
#   and rain forecasts into weekly bins, computes logit-scale climatology
#   features, and derives rain-based predictors (diff_*, min_*_10day).
#
# Inputs
#   - YAML spec in specs/2025_blend/ (e.g., connect_mok.yml)
#   - RDS path specified in spec$input_rds
#
# Outputs
#   - RDS path specified in spec$output_rds
#
# Usage
#   Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_mok
#
# Dependencies
#   - R/2025_blending_process/connect_utils.R
#   - R/_shared/misc.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(tidyr)
  library(optparse)
  library(yaml)
})

source("R/2025_blending_process/connect_utils.R")
source("R/_shared/misc.R")

# ---- CLI argument parsing ----
option_list <- list(
  make_option("--spec_id", type = "character", default = "connect_mok",
              help = "Spec file name (without .yml) in specs/2025_blend/ [default %default]")
)
if (!interactive()) {
  opt <- parse_args(OptionParser(option_list = option_list))
  spec_id <- opt$spec_id
} else {
  spec_id <- "connect_no_mok_filter"
}

# ---- Load spec ----
spec_path <- file.path("specs", "2025_blend", paste0(spec_id, ".yml"))
spec <- yaml::read_yaml(spec_path)

# ---- Run conversion ----
wide_df_created <- make_cv_rds_from_daylevel(spec = spec)

cat("Wrote:", spec$output_rds, "\n")
cat("Rows:", nrow(wide_df_created), " Cols:", ncol(wide_df_created), "\n")
print(names(wide_df_created))
