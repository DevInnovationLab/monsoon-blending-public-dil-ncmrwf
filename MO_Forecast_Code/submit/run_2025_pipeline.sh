#!/bin/bash
# ==============================================================================
# run_2025_pipeline.sh
# ==============================================================================
# End-to-end 2025 blending pipeline (weekly-bin multinomial via nnet).
# Assumes prepare_data has already been run and
# Monsoon_Data/Processed_Data/Combined/*.rds exists.
#
# Runs:
#   - Connect/prepare weekly data (mok, clim_mok_date, no_mok_filter)
#   - CV evaluation (mok + variants)
#   - Hindcast evaluation (mok only)
#   - 2025 out-of-sample evaluation
#   - Figures
#
# Usage:
#   bash submit/run_2025_pipeline.sh
# ==============================================================================

set -euo pipefail

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "Step 0: Convert daily data -> weekly RDS (all modes)"
Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_mok_inc_ncum
log "  0_connect (mok) done"


log "Step 1: CV evaluation (mok)"
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id cv_models_inc_ncum
log "  CV 2000-2015 (mok) done"

log "2025 pipeline complete"