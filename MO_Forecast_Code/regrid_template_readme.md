# Cell Transform Weights File

The pipeline can regrid rainfall from one spatial grid to another using a user-supplied CSV of linear combination weights. This document describes the required file structure.

## Overview

The pipeline applies the transform as:

```
rainfall(target_cell) = sum over source cells( weight_i * rainfall(source_cell_i) )
```

The weights file defines which source cells contribute to each target cell and with what coefficient. The regridding method itself (conservative overlap, nearest-neighbor, bilinear, etc.) is up to the user -- the pipeline only requires the final weights in the format below.

## File location

Save to the path specified in the YAML spec's `options.cell_transform_file` field. For the quarter-degree pipeline this is:

```
Monsoon_Data/weights/equal_overlap_weights_2_to_25.csv
```

## Required columns

The CSV must have exactly three columns:

| Column      | Type    | Description                                              |
|-------------|---------|----------------------------------------------------------|
| `target_id` | string  | Grid cell ID on the output (target) grid                 |
| `source_id` | string  | Grid cell ID on the input (source) grid                  |
| `weight`    | numeric | Coefficient for this source cell's contribution          |

No other columns should be present.

## ID format

Grid cell IDs throughout the pipeline are constructed as `paste0(lat, "_", lon)`, e.g., `"8.25_72.75"` or `"9_73"`.

- **`source_id`** values must exactly match the lat/lon coordinates in the source NetCDF files (after any `dimensions.rename` mapping in the spec is applied).
- **`target_id`** values must exactly match the IDs used in the thresholds file for the target grid.

Be careful with formatting details like trailing zeros (`8.0` vs `8`) or decimal precision (`22.5` vs `22.50`) -- these must match exactly.

## Row structure

One row per source-target pair. A target cell that receives contributions from multiple source cells will have multiple rows (one per contributing source). A target cell that maps to a single source cell has one row.

Example:

```csv
target_id,source_id,weight
8.125_72.125,9_73,1.0
8.125_72.375,9_73,1.0
8.375_72.125,9_73,1.0
10.125_74.125,11_75,0.5
10.125_74.125,11_73,0.5
```

In this example, the first three target cells each map entirely to one source cell (`weight = 1.0`), while the last target cell straddles two source cells (each contributing `weight = 0.5`).

## Weight constraint

**For each `target_id`, the weights must sum to 1.0.** This ensures the transformation preserves rainfall values in expectation.

## Where this file is referenced

- `specs/raw_data/ngcm_quarter_degree.yml` -- regrid NGCM forecasts (2-degree to 0.25-degree)
- `specs/raw_data/aifs_quarter_degree.yml` -- regrid AIFS forecasts (2-degree to 0.25-degree)
- IMD ground truth (`imd_quarter_degree.yml`) does **not** use this file because IMD data is already at 0.25-degree native resolution.

The relevant YAML fields are:

```yaml
options:
  cell_transform_enabled: true
  cell_transform_file: "Monsoon_Data/weights/equal_overlap_weights_2_to_25.csv"
```

## Code reference

The transform is applied by `transform_forecast_wide()` and `transform_groundtruth_long()` in `R/prepare_data/nc_utils.R`. The weights file is read by `read_cell_transform()` in the same file.
