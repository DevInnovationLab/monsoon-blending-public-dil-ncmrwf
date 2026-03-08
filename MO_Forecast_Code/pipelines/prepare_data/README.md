# prepare_data Pipeline

Converts raw NetCDF rainfall data (IMD ground truth, NeuralGCM and AIFS forecasts) into modeling-ready wide tables with onset probabilities, climatology forecasts, and combined features.

## Stages

### Stage 1: `1_process_raw_nc_files.R`

Reads raw NetCDF files, computes rolling-window onset probabilities per grid cell, and writes wide-format RDS tables.

- **Forecast systems** (e.g., `--spec_id ngcm`): Produces per-ensemble-member onset probabilities across lead days.
- **Ground truth** (e.g., `--spec_id imd`): Produces observed onset dates and long/wide rainfall tables.

### Stage 2: `2_build_climatology.R`

Fits per-cell KDE climatology models from historical IMD onset dates and produces issue-date probability forecasts over lead days. Supports multiple training windows (e.g., 1900-2024, 2000-2024) defined in the spec. One convention to note is that for unconditional climatology ("unc_clim_raw") a "day 0" forecast is used to store the probability that onset occurred before the forecast date.  

### Stage 3: `3_combine_datasets.R`

Joins ground truth, climatology, and forecast system outputs into a single wide table per combination template. Aligns time/space grids, adds "plus" remainder bins, and optionally trims post-onset forecasts.

## Running

All commands run from the repository root (`MO_Forecast_Code/`):

```bash
# Ground truth
Rscript pipelines/prepare_data/1_process_raw_nc_files.R --spec_id imd
Rscript pipelines/prepare_data/2_build_climatology.R --spec_id imd

# Forecast systems
Rscript pipelines/prepare_data/1_process_raw_nc_files.R --spec_id ngcm
Rscript pipelines/prepare_data/1_process_raw_nc_files.R --spec_id aifs

# Combine
Rscript pipelines/prepare_data/3_combine_datasets.R --spec_id combine_template_2025

# Or run everything via shell wrapper:
bash submit/run_prepare_data.sh
```

## Spec Files

- **`specs/raw_data/*.yml`**: One per data source. Key fields: `type` (ground_truth_rainfall / rainfall_forecast), `input.nc_folder`, `output.out_dir`, `thresholds`, `mok`, `paths.climatology_out_dir`, `climatologies`.
- **`specs/combine/*.yml`**: One per combination template. Key fields: `input.ground_truth_wide_rds`, `input.climatologies`, `forecasts`, `output.out_dir`.

## Inputs

- Raw NetCDF files in `Monsoon_Data/raw_nc/` (paths configured in `specs/raw_data/*.yml`):
  - `IMD_2by2/` — IMD ground truth at 2-degree (files: `data_YYYY.nc`)
  - `IMD_25by25/` — IMD ground truth at 0.25-degree
  - `ngcm/` — NeuralGCM forecast NetCDF files
  - `aifs/` — AIFS forecast NetCDF files
- Reference files in `Monsoon_Data/reference/`:
  - `thresholds_df.csv` — per-cell onset thresholds (2-degree)
  - `thresholds_quarter_degree_df.csv` — per-cell onset thresholds (0.25-degree)
  - `MOK Onset.csv` — observed MOK dates
  - `MOK Onset June.csv` — fixed June 1 MOK dates
  - `MOK Onset Any.csv` — no MOK filter variant
  - `allowed_cells.csv` — cells used for evaluation

## Outputs

- `Monsoon_Data/Processed_Data/Models/*.rds` -- Per-system wide onset tables
- `Monsoon_Data/Processed_Data/Climatology/*.rds` -- KDE climatology issue-date forecasts
- `Monsoon_Data/Processed_Data/Combined/*.rds` -- Merged modeling-ready wide tables

## NetCDF Input Requirements

### General Requirements (all files)

- **Format**: NetCDF-3 or NetCDF-4 (`.nc` / `.nc4`), readable by `ncdf4::nc_open()`
- **Filename convention**: Must contain a 4-digit year matching `(19|20)\d{2}` (e.g., `data_2000.nc`, `2000.nc`). The year is extracted from the filename, not from the file contents.
- **Coordinate dimensions**: Must include latitude, longitude, and time dimensions. Common name variants are handled automatically:
  - Latitude: `lat`, `latitude`, `LATITUDE`, `LAT`
  - Longitude: `lon`, `longitude`, `LONGITUDE`, `LON`, `long`, `LONG`
  - Time: `time`, `TIME`, `Time`, `t`, `T`
  - Any remaining mismatches can be fixed via `dimensions.rename` in the spec
- **Time encoding**: Numeric time values must have a CF-compliant `units` attribute (e.g., `"days since 1900-01-01"`). Supported units: seconds, minutes, hours, days.
- **Rainfall variable**: A single variable specified by `input.value_col` in the spec (e.g., `"tp"`, `"RAINFALL"`, `"precipitation_cumulative_mean"`). Must be a numeric array dimensioned by the coordinate dimensions.

### Forecast NetCDF Requirements (`type: "rainfall_forecast"`)

In addition to the general requirements:

- **Lead-day dimension**: A dimension for forecast lead days (e.g., `day`), specified via `input.wide_day_dim` in the spec. Values should be positive integers representing lead day numbers.
- **Ensemble dimension** (optional): A `number` dimension for ensemble members. If present, onset probabilities are computed per member and then aggregated (mean, sd, quantiles). The spec can optionally filter to a maximum number of members via `filter.max_number`.
- **Day range**: The spec's `options.min_day` and `options.max_day` control which lead days are kept. The file should contain at least days in `[min_day, max_day + window - 1]`.

### Ground Truth NetCDF Requirements (`type: "ground_truth_rainfall"`)

In addition to the general requirements:

- **No lead-day dimension**: The rainfall variable should be dimensioned by `(lat, lon, time)` only -- one value per grid cell per calendar day.
- **Daily resolution**: Time steps must be daily. The pipeline computes rolling-window rainfall sums for onset detection.
- **Coverage**: Should cover the full monsoon season (at least from `options.cutoff_month_day` through the end of the season) for each year.

### Optional: Cell Transform (Regridding)

When `options.cell_transform_enabled: true`, a CSV weights file is required with columns:
- `target_id`: Grid cell ID on the target grid (format `"lat_lon"`)
- `source_id`: Grid cell ID on the source grid
- `weight`: Linear combination weight (weights for a given target_id should sum to 1)

This is applied before onset detection, so thresholds must match the target grid.

## Data Dictionary

### Stage 1 Outputs (`Processed_Data/Models/`)

#### Forecast wide table (`<spec_id>_wide.rds`)

One row per (grid cell, initialization date, year). Ensemble members are aggregated.

| Column | Type | Description |
|--------|------|-------------|
| `id` | character | Grid cell ID (`"lat_lon"`, e.g., `"22_78"`) |
| `time` | Date | Forecast initialization date |
| `year` | integer | Year (from filename) |
| `onset_thresh` | numeric | Onset threshold for this grid cell (mm) |
| `mok_date` | Date | Monsoon Onset Kerala date for this year (or NA) |
| `forecast_rain_day_<k>` | numeric | Ensemble mean rainfall on lead day k |
| `forecast_rain_sd_day_<k>` | numeric | Ensemble std dev of rainfall on lead day k |
| `frac_raining_day_<k>` | numeric | Fraction of members with rainfall > 1mm on day k |
| `predicted_prob_day_<k>` | numeric | Onset probability on day k (raw, no start restriction) |
| `predicted_prob_clim_mok_date_day_<k>` | numeric | Onset probability on day k (onset restricted to after June 2) |
| `predicted_prob_mok_day_<k>` | numeric | Onset probability on day k (onset restricted to after MOK date) |

#### Ground truth wide table (`<spec_id>_wide.rds`)

One row per (grid cell, year).

| Column | Type | Description |
|--------|------|-------------|
| `id` | character | Grid cell ID (`"lat_lon"`) |
| `year` | integer | Year |
| `mr_onset_idx` | numeric | Index (position in daily series) of Moron-Robertson onset |
| `mr_onset_date` | Date | Calendar date of onset |
| `mr_onset_day` | numeric | Days from `cutoff_month_day` to onset (e.g., days since May 1) |
| `cutoff_date` | Date | Season start date for that year (from `options.cutoff_month_day`) |

#### Ground truth long table (`<spec_id>_long.rds`)

One row per (grid cell, day). Annotated daily rainfall series.

| Column | Type | Description |
|--------|------|-------------|
| `id` | character | Grid cell ID |
| `time` | Date | Calendar date |
| `year` | integer | Year |
| `<value_col>` | numeric | Daily rainfall (mm), variable name from spec (e.g., `rainfall`) |
| `onset_thresh` | numeric | Onset threshold for this cell |
| `mok_date` | Date | MOK date for this year |
| `mr_onset_date` | Date | Onset date for this cell-year |
| `mr_onset_flag` | logical | TRUE on the onset date, FALSE otherwise |

### Stage 2 Outputs (`Processed_Data/Climatology/`)

#### Climatology forecast table (`<out_stem>.rds`)

One row per (grid cell, issue date, lead day). Produced per climatology training window.

| Column | Type | Description |
|--------|------|-------------|
| `lat`, `lon` | numeric | Grid cell coordinates |
| `time` | Date | Issue (forecast) date |
| `day` | integer | Lead day (1 to `forecast_window`) |
| `predicted_prob` | numeric | P(onset on lead day k \| onset hasn't occurred by issue date), or unconditional mass if `conditional: false` |
| `model` | character | Climatology model label (e.g., `"kde"`) |

### Stage 3 Outputs (`Processed_Data/Combined/`)

#### Combined wide table (`<spec_id>_combined_wide.rds`)

One row per (grid cell, initialization date, year). All systems merged.

| Column | Type | Description |
|--------|------|-------------|
| `lat`, `lon` | numeric | Grid cell coordinates |
| `id` | character | Grid cell ID |
| `time` | Date | Forecast initialization date |
| `year` | integer | Year |
| `true_onset_day` | numeric | Days from season start to observed onset |
| `true_onset_date` | Date | Observed onset date |
| `clim_p_onset_day_<k>` | numeric | Conditional climatology (evolving-expectations model) onset probability for day k |
| `clim_unc_p_onset_day_<k>` | numeric | Unconditional climatology onset probability for day k |
| `<system>_p_onset_day_<k>` | numeric | Forecast system onset probability for day k (e.g., `ngcm_p_onset_day_1`) |
| `<system>_rain_mean_day_<k>` | numeric | Forecast system mean rainfall for day k |
| `<system>_onset_thresh` | numeric | Onset threshold (constant per cell) |
| `<system>_mok_date` | Date | MOK date (constant per year) |
| `*_day_<max+1>plus` | numeric | Remainder bin probability (1 - sum of day bins) |

Where `<system>` is `ngcm`, `aifs`, etc. and `<k>` ranges from 1 to `max_day` (typically 28 or 40).

## Adding a New Data Source

1. Create a new spec in `specs/raw_data/` (copy an existing forecast spec as template).
2. Set `input.nc_folder`, `value_col`, `dimensions.rename` for your NetCDF structure.
3. Run Stage 1 with `--spec_id <your_spec>`.
4. Add the new source to the relevant `specs/combine/*.yml` under `forecasts:`.
5. Re-run Stage 3.

## Adding a New Onset Filter Variant

1. Create a new IMD spec variant in `specs/raw_data/` (e.g., `imd_new_variant.yml`) with the appropriate `mok:` configuration.
2. Run Stages 1-2 with the new spec.
3. Create a matching combine template in `specs/combine/`.
4. Run Stage 3 with the new combine spec.
