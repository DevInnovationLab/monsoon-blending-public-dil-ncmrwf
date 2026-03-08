# 2025 Weekly-Bin Blending Pipeline

Cross-validated weekly-bin multinomial onset blending using `nnet`. Combines climatology, NeuralGCM, and AIFS forecasts into optimized multi-model ensembles, evaluates skill, and produces publication figures.

## Stages

### Stage 0: `0_connect_prepare_data_to_2025_pipeline.R`

Converts day-level wide RDS from the prepare_data pipeline into weekly-bin format. Aggregates daily onset probabilities into 4 weekly bins plus a "later" bin, computes logit-scale climatology features, and derives rain-based predictors.

Driven by YAML specs in `specs/2025_blend/connect_*.yml`:
- `--spec_id connect_mok`: Uses observed MOK date per year
- `--spec_id connect_clim_mok_date`: Uses fixed climatological MOK date (June 1)
- `--spec_id connect_no_mok_filter`: No MOK-based filter

Each connect spec defines `mode`, `input_rds`, `output_rds`, `forecast_models` (with rain predictor types), and `climatology` (with optional window tags).

### Stage 1: `1_blend_evaluation.R`

Main cross-validation engine. Fits weekly-bin multinomial models (via `nnet`) using formulas defined in YAML specs. Supports multiple CV methods (global, local, neighbors, cluster). Computes metrics (Brier, RPS, AUC), reliability plots, and optimized multi-model ensemble (MME) weights. The summary_models outputs in the results folder store the primary metrics of interest. 

### Stage 2: `2_2025_evaluation.R`

Out-of-sample evaluation of 2025 monsoon onset forecasts. Scores the blended model, climatologies, and Platt-calibrated NGCM forecasts against three ground-truth variants. Stage 2 is only designed to run on the paper's dataset.

### Stage 3: `3_produce_figures.R`

Reads pre-computed model summaries and produces publication-ready figures: overall metrics by period, climatology training window variation, model skill comparisons, weekly performance, yearly time series, reliability diagrams, and spatial skill maps. Stage 3 is only designed to run on the paper's dataset. 

## Running

All commands run from the repository root (`MO_Forecast_Code/`):

```bash
# Stage 0: Convert daily -> weekly (all modes)
Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_mok
Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_clim_mok_date
Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_no_mok_filter

# Stage 1: CV evaluation
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id cv_models
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id hindcast_1965_1978

# Stage 1 variants (clim_mok_date, no_mok_filter)
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id cv_models_clim_mok_date
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id cv_models_no_mok_filter

# Stage 2: 2025 out-of-sample evaluation
Rscript pipelines/2025_blending_process/2_2025_evaluation.R

# Stage 3: Publication figures
Rscript pipelines/2025_blending_process/3_produce_figures.R

# Or run everything via shell wrapper:
bash submit/run_2025_pipeline.sh
```

## Spec Files

Located in `specs/2025_blend/*.yml`. Key fields:

- `run.cutoff_mode`: MOK filter mode (mok / clim_mok_date / no_mok_filter)
- `run.training_years`, `run.true_holdout_years`: Year splits
- `cv.methods`: CV strategy (global, local, neighbors, cluster2, cluster3)
- `models`: Named list of model formulas (each maps to an `nnet` multinomial logistic regression)
- `extras.clim_logits`: Climatology-based baseline model definitions
- `extras.forecasts`: Raw/calibrated forecast model definitions
- `mme`: Multi-model ensemble optimization settings (variants, `blend_models` list, weight optimization)

## Inputs

- `Monsoon_Data/Processed_Data/Combined/*.rds` (from prepare_data pipeline)
- `Monsoon_Data/dissemination_cells.csv`
- `Monsoon_Data/evaluation_2025/*.csv` (for Stage 2)

## Outputs

- `Monsoon_Data/results/2025_model_evaluation/` -- CV metrics, reliability plots, blend weights, map inputs
- `Monsoon_Data/results/2025_model_evaluation/evaluation/` -- 2025 out-of-sample metrics
- `figures/` -- Publication-ready PDF/SVG plots

## Adding a New Model Formula

1. Open the relevant spec in `specs/2025_blend/` (e.g., `cv_models.yml`).
2. Add a new entry under `models:` with a name and formula string.
3. The formula uses column names from the weekly-bin wide data (e.g., `ngcm_p_onset_week1`, `clim_logit_week1`).
4. Re-run Stage 1.

## Adding Window Variants

Configure in the spec under `extras.clim_logits`:
- Set `window_start_years` to a list of start years
- Set `window_end_year` to the end year
- Each start year produces a separate model variant (e.g., `clim_raw_1960_2024`)

## Running Hindcast Evaluation

Use a separate spec that sets `run.true_holdout_years` to the hindcast period (e.g., 1965-1978) and `mme.weights_year_tag` to read pre-computed MME weights from the main CV run:

```bash
Rscript pipelines/2025_blending_process/1_blend_evaluation.R --spec_id hindcast_1965_1978
```
