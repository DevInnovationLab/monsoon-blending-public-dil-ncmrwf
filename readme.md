# Replication Package: Designing probabilistic AI monsoon forecasts to inform agricultural decision-making
This is a fork of the monsoon-blending-public repository, modified to perform an NCUM/AIFS blnd instead of AIFS/NeuralGCM. Please contact Colin Aitken (caitken@uchicago.edu) for data to run this code if needed. 

Replication code for probabilistic Indian monsoon onset forecasts that blend rainfall observations (IMD) with AI weather prediction model forecasts (NeuralGCM, AIFS) through multinomial blending models. This package reproduces all results from data preparation through cross-validated evaluation to publication figures. See ``Running the Pipeline'' below for more information.

The data used in the original paper can be found at https://doi.org/10.5281/zenodo.18894299. Before running this pipeline, all three folders found in this dataset need to be added to Monsoon_Data/raw_nc.

## Adapting This Code to Your Own Data

Most of this pipeline is designed to be reusable beyond the specific datasets used in the paper. All pipeline behavior is controlled by YAML spec files rather than hardcoded values, so you can point the pipeline at your own data sources without modifying any R code. The code that produces figures is currently designed to produce the figures appearing in the paper, so will need further modification. The summary_global files in the results folder contain all of our primary metrics. 

- **Different ground truth rainfall**: Create a new spec in `specs/raw_data/` modeled after `imd.yml`. Point `input.nc_folder` at your NetCDF files and configure the variable name, dimension mappings, and onset thresholds for your grid.
- **Different AI forecast models**: Create a new spec in `specs/raw_data/` modeled after `ngcm.yml` (for ensembles) or `aifs.yml`. The pipeline handles ensemble forecasts in NetCDF format with configurable dimension names and variable mappings.
- **Different grid resolutions**: The pipeline supports regridding via weighted linear combinations. Set `cell_transform_enabled: true` in your raw data spec and provide a weights file mapping source grid cells to target grid cells.
- **Different onset definitions or time periods**: Edit the `climatologies` section of your ground truth spec to change training/test year ranges, onset column names, and season boundaries.
- **Different blending model formulas**: Edit the `models.formulas` section of the `specs/2025_blend/cv_models*.yml` specs to change which predictors enter the multinomial logistic regression.
- **Different forecast systems in the MME**: Edit the `mme.blend_models` list in `specs/2025_blend/cv_models*.yml` to add, remove, or reconfigure which calibrated forecast systems are blended.

The fully commented spec files in `specs/` serve as the primary reference for all configurable options. See [Spec-Driven Design](#spec-driven-design) below for details on the architecture.

If you change any of the inputs, you will need to add new "combine" and "2025_blend" spec files to ensure they are included in the combined climatology+ground truth + weather model dataset, and that their outputs are included in the blending code. 

## Repository Layout

```
MO_Forecast_Code/
├── R/
│   ├── _shared/                  Core utilities (spec parsing, scoring, misc)
│   │   ├── misc.R                  Null-coalescing, softmax, logit helpers
│   │   └── read_spec.R             Spec loading and validation
│   ├── prepare_data/             Data preparation helpers
│   │   ├── nc_utils.R              NetCDF reading, regridding, wide-table construction
│   │   ├── onset_utils.R           Onset detection, MOK dates, threshold loading
│   │   ├── climatology_utils.R     KDE fitting, climatological forecasts
│   │   └── combine_forecasts_utils.R  Merging climatology + forecast families
│   └── 2025_blending_process/    Blending pipeline helpers
│       ├── connect_utils.R         Day-to-week aggregation, logit transforms
│       ├── blend_evaluation_utils.R  CV evaluation, nnet multinomial, Platt calibration
│       ├── evaluation_2025_utils.R   Out-of-sample scoring utilities
│       └── blend_figure_utils.R      Figure generation utilities
├── pipelines/                     Each pipelines folder has its own README.md file for more information 
│   ├── prepare_data/             Stages 1-3: NetCDF -> onset tables -> climatology -> combined. 
│   └── 2025_blending_process/    Stages 0-3: weekly-bin blending + evaluation + figures
├── specs/
│   ├── raw_data/                 NetCDF input specs (aifs, ngcm, imd variants)
│   ├── combine/                  Data combination templates (for combining different weather model + climatology combinations)
│   └── 2025_blend/               Blended model specs (formulas, MME config, connect specs)
├── submit/                       Shell wrappers and PBS job scripts
├── Monsoon_Data/                 Data directory (not tracked in git)
│   ├── raw_nc/                     Raw NetCDF inputs (IMD, NGCM, AIFS)
│   ├── reference/                  Onset thresholds, MOK dates
│   ├── Processed_Data/
│   │   ├── Models/                   Per-system onset tables
│   │   ├── Climatology/              KDE climatology forecasts
│   │   ├── Combined/                 Merged modeling-ready wide tables
│   │   └── 2025_pipeline_input/      Weekly-bin data for blending
│   ├── results/
│   │   └── 2025_model_evaluation/    Model metrics, blend weights, figures
│   ├── evaluation_2025/            Out-of-sample forecast + ground truth files
│   └── maps/                       India boundary data for maps
└── figures/                      Figures (PDF/SVG). Note that PDF figures will only appear if Cairo is installed (for Arial fonts)
```

## Prerequisites

### Required Input Data

The pipeline expects the following data to be in place before running:

| Path | Description |
|------|-------------|
| `Monsoon_Data/raw_nc/IMD_2by2/` | IMD gridded rainfall NetCDF files (`data_YYYY.nc`), one per year |
| `Monsoon_Data/raw_nc/ngcm/` | NeuralGCM ensemble forecast NetCDF files, one per year |
| `Monsoon_Data/raw_nc/aifs/` | AIFS ensemble forecast NetCDF files, one per year |
| `Monsoon_Data/reference/thresholds_df.csv` | Per-grid-cell onset rainfall thresholds (`lat`, `lon`, `onset_threshold`) |
| `Monsoon_Data/reference/MOK Onset.csv` | Monsoon Onset Kerala dates by year (`Year`, `MOK`) |
| `Monsoon_Data/maps/` | India boundary shapefiles (for map figures) |
| `Monsoon_Data/evaluation_2025/` | Out-of-sample forecast and ground truth CSVs (for stage 2) |

### R Dependencies

```r
install.packages(c(
  # Core
  "optparse", "yaml", "data.table", "dplyr", "tidyr", "purrr", "stringr",
  "lubridate", "tibble", "tidyverse", "conflicted",
  # Parallel
  "furrr", "future", "future.apply", "parallel",
  # Data I/O
  "ncdf4", "readr",
  # Modeling
  "nnet",
  # Scoring
  "pROC", "ModelMetrics", "matrixStats",
  # Visualization
  "ggplot2", "sf", "forcats", "colorspace", "patchwork", "maps",
  # Utilities
  "tictoc", "here"
))
# Optional (used only if available):
# install.packages(c("svglite", "R.matlab"))
```

## Running the Pipeline

All scripts must be run from the repository root (`MO_Forecast_Code/`). Scripts use relative `source("R/...")` paths that break from other working directories. On a PBS server, the entire pipeline on the original dataset can be run using  qsub submit/run_full_2025.pbs from the MO_Forecast_Code directory. The two .sh files in that folder will run Stage 1 and Stage 2 separately.  The spec_ids called will need to be modified if you use new data sources.


### Stage 1: Prepare Data

Processes raw NetCDF files into onset tables, fits KDE climatologies, and combines everything into modeling-ready wide tables.

```bash
# Process raw NetCDF files into per-system onset tables
Rscript pipelines/prepare_data/1_process_raw_nc_files.R --spec_id imd
Rscript pipelines/prepare_data/1_process_raw_nc_files.R --spec_id ngcm
Rscript pipelines/prepare_data/1_process_raw_nc_files.R --spec_id aifs

# Build KDE climatology forecasts from IMD onset dates
Rscript pipelines/prepare_data/2_build_climatology.R --spec_id imd

# Combine ground truth, model forecasts, and climatology into wide tables
Rscript pipelines/prepare_data/3_combine_datasets.R --spec_id combine_template_2025
```

**Outputs**: `Monsoon_Data/Processed_Data/Combined/combine_template_2025_combined_wide.rds`

### Stage 2: Blending Pipeline

Converts daily probabilities to weekly bins, runs cross-validated multinomial model evaluation and multi-model ensemble (MME) weight optimization, scores against out-of-sample 2025 data, and produces publication figures.

```bash
# Convert daily data to weekly bins (one call per onset filter variant)
Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_mok
Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_clim_mok_date
Rscript pipelines/2025_blending_process/0_connect_prepare_data_to_2025_pipeline.R --spec_id connect_no_mok_filter

# Cross-validated model evaluation + MME weight optimization
Rscript pipelines/2025_blending_process/1_blend_evaluation.R

#If you have modified the pipeline to use new models, stop here and look in the Monsoon_Data/results folder! The next two scripts
# will not function correctly

# Out-of-sample 2025 evaluation
Rscript pipelines/2025_blending_process/2_2025_evaluation.R

# Publication figures
Rscript pipelines/2025_blending_process/3_produce_figures.R
```

**Outputs**: `Monsoon_Data/results/2025_model_evaluation/` (metrics CSVs, blend weights, reliability plots, publication figures)

### End-to-End via Shell Wrappers

```bash
bash submit/run_prepare_data.sh
bash submit/run_2025_pipeline.sh
```

### PBS Cluster Submission

```bash
qsub submit/run_full_2025.pbs
```

## Data Flow

```
Raw NetCDF (Monsoon_Data/raw_nc/)
    │
    ▼
┌──────────────────────────────────┐
│  1_process_raw_nc_files.R        │  specs/raw_data/*.yml
│  → Processed_Data/Models/*.rds   │
└──────────┬───────────────────────┘
           │
           ├──▶ 2_build_climatology.R
           │    → Processed_Data/Climatology/*.rds
           │                 │
           ▼                 ▼
┌──────────────────────────────────┐
│  3_combine_datasets.R            │  specs/combine/*.yml
│  → Processed_Data/Combined/*.rds │
└──────────┬───────────────────────┘
           │
           ▼
┌──────────────────────────────────────────┐
│  0_connect (day → week bins)             │  specs/2025_blend/connect_*.yml
│  → Processed_Data/2025_pipeline_input/   │
└──────────┬───────────────────────────────┘
           │
           ▼
┌──────────────────────────────────────────┐
│  1_blend_evaluation.R                    │  specs/2025_blend/cv_models*.yml
│  Cross-validated nnet multinomial        │  specs/2025_blend/hindcast_*.yml
│  Platt calibration, MME optimization     │
│  → results/2025_model_evaluation/        │
└──────────┬───────────────────────────────┘
           │
           ├──▶ 2_2025_evaluation.R
           │    Out-of-sample scoring (Brier, RPS, AUC)
           │                 │
           ▼                 ▼
┌──────────────────────────────────────────┐
│  3_produce_figures.R                     │
│  → figures/ (PDF/SVG)                    │
└──────────────────────────────────────────┘
```

## Onset Filter Variants

Three onset filter variants control which data points are included in training and evaluation, depending on how the monsoon onset date is defined relative to the Monsoon Onset Kerala (MOK):

| Variant | Spec suffix | Description |
|---------|-------------|-------------|
| **mok** | _(default)_ | Only issue dates after the observed MOK date each year |
| **clim_mok_date** | `_clim_mok_date` | Only issue dates after a fixed climatological MOK date (June 1) |
| **no_mok_filter** | `_no_mok_filter` | No MOK-based filtering; all issue dates from May 1 onward |

Each variant has its own set of connect, CV, and hindcast specs in `specs/2025_blend/`.

## Spec-Driven Design

All pipeline behavior is controlled by YAML specs in `specs/`. Specs define input paths, variable selection, modeling options, and output configuration. Output basenames are derived from the spec filename (`spec_id`), not from a field in the YAML. Pipeline scripts are thin orchestration layers: parse args -> load spec -> call helpers -> write artifacts.

### Spec directories

| Directory | Purpose | Used by |
|-----------|---------|---------|
| `specs/raw_data/` | NetCDF input configuration (paths, variables, thresholds, MOK) | `1_process_raw_nc_files.R`, `2_build_climatology.R` |
| `specs/combine/` | Which processed datasets to merge into wide tables | `3_combine_datasets.R` |
| `specs/2025_blend/connect_*.yml` | Day-to-week conversion settings (forecast models, rain predictors) | `0_connect_prepare_data_to_2025_pipeline.R` |
| `specs/2025_blend/cv_models*.yml` | Model formulas, MME configuration, forecast calibration | `1_blend_evaluation.R` |
| `specs/2025_blend/hindcast_*.yml` | Out-of-sample hindcast evaluation configuration | `1_blend_evaluation.R` |

### Key spec sections (cv_models*.yml)

- **`models.formulas`**: Named multinomial logistic regression formulas using `nnet::multinom`
- **`mme`**: Multi-model ensemble configuration — which calibrated models to blend, weight optimization settings
- **`mme.blend_models`**: List of models entering the MME, each with a `name`, `source` (clim or forecast), and optional `cal_variant`
- **`extras.forecasts`**: Per-system Platt calibration and raw/calibrated scoring options
- **`extras.clim_logits`**: Climatology baseline configurations with optional training-window variants

## Conventions

- Spatial key: `id = paste0(lat, "_", lon)`
- Time columns: `time` (POSIXct/Date), `year`
- Outcome categories: `week1` through `week4` plus `later` (5 weekly bins)
- Forecast probabilities stored in wide format with system-specific prefixes
- All scripts must be run from the repository root (`MO_Forecast_Code/`)
