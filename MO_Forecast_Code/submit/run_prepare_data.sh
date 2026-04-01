#!/bin/bash
# ==============================================================================
# run_prepare_data.sh
# ==============================================================================
# Run full prepare_data pipeline for ALL required specs:
#   - IMD ground truth (mok, clim_mok_date, no_mok_filter variants)
#   - NGCM and AIFS forecasts
#   - Combine datasets (standard, clim_mok_date, no_mok_filter)
#
# Usage:
#   bash submit/run_prepare_data.sh
# ==============================================================================

set -euo pipefail

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "Stage 1+2: IMD ground truth"
Rscript pipelines/prepare_data/1_process_raw_nc_files.R --spec_id imd
log "  1_process_raw_nc_files (imd) done"
Rscript pipelines/prepare_data/2_build_climatology.R --spec_id imd
log "  2_build_climatology (imd) done"


log "Stage 1: NCUM forecasts"
Rscript pipelines/prepare_data/1_process_raw_nc_files.R --spec_id ncmrwf
log "  1_process_raw_nc_files (ncum) done"

log "Stage 1: AIFS forecasts"
Rscript pipelines/prepare_data/1_process_raw_nc_files.R --spec_id aifs
log "  1_process_raw_nc_files (aifs) done"

log "Stage 3: Combine datasets"
Rscript pipelines/prepare_data/3_combine_datasets.R --spec_id combine_template_2025_inc_ncum
log "  3_combine (combine_template_2025) done"
log "Prepare data complete"
