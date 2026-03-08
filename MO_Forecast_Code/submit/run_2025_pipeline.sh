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
Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_mok
log "  0_connect (mok) done"
Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_clim_mok_date
log "  0_connect (clim_mok_date) done"
Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_no_mok_filter
log "  0_connect (no_mok_filter) done"

log "Step 1: CV evaluation (mok)"
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id cv_models
log "  CV 2000-2024 (mok) done"

log "Step 1b: Hindcast evaluation (mok)"
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id hindcast_1965_1978
log "  Hindcast 1965-1978 (mok) done"

log "Step 1c: CV evaluation (clim_mok_date)"
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id cv_models_clim_mok_date
log "  CV 2000-2024 (clim_mok_date) done"

log "Step 1d: Hindcast evaluation (clim_mok_date)"
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id hindcast_1965_1978_clim_mok_date
log "  Hindcast 1965-1978 (clim_mok_date) done"

log "Step 1e: CV evaluation (no_mok_filter)"
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id cv_models_no_mok_filter
log "  CV 2000-2024 (no_mok_filter) done"

log "Step 1f: Hindcast evaluation (no_mok_filter)"
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id hindcast_1965_1978_no_mok_filter
log "  Hindcast 1965-1978 (no_mok_filter) done"


log "Step 2: 2025 out-of-sample evaluation"
Rscript pipelines/2025_blending_process/2_2025_evaluation.R
log "  2_2025_evaluation done"

log "Step 3: Produce figures"
Rscript pipelines/2025_blending_process/3_produce_figures.R
log "  3_produce_figures done"

log "2025 pipeline complete"