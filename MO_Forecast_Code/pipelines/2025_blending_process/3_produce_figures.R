# ==============================================================================
# File: 3_produce_figures.R
# ==============================================================================
# Purpose
#   Reads pre-computed model evaluation summaries (RDS) and calibration chart
#   data (CSV) and produces publication-ready figures (PDF + SVG) for the
#   blended monsoon onset forecast paper.
#
# Inputs
#   - Monsoon_Data/results/2025_model_evaluation/summary_models_*.rds
#   - Monsoon_Data/results/2025_model_evaluation/yearly_metrics_global_*.rds
#   - Monsoon_Data/results/2025_model_evaluation/evaluation/model_metrics_*.rds
#   - Monsoon_Data/results/2025_model_evaluation/calibration plots/*.csv
#   - Monsoon_Data/results/2025_model_evaluation/cell_metrics_*.rds
#   - Monsoon_Data/maps/india_boundary.csv
#
# Outputs
#   - figures/*.pdf, figures/*.svg (overall, skill, weekly, yearly, reliability,
#     maps)
#
# Usage
#   Rscript pipelines/2025_blending_process/3_produce_figures.R
#
# Dependencies
#   - R/2025_blending_process/blend_figure_utils.R
# ==============================================================================

# ---------------------- 0. USER SETTINGS ----------------------

# Root folder for outputs (does not change working directory)
OUTPUT_DIR <- file.path(getwd(), "figures")

# Plot theme defaults
FONT_SIZE <- 16

# Vector formats to save by default
# cairo_pdf requires the Cairo graphics device; fall back to SVG-only on systems without it
HAS_CAIRO <- capabilities("cairo")
if (!HAS_CAIRO) {
  message("Note: Cairo not available — PDF output disabled, producing SVG only.")
  VECTOR_FORMATS_DEFAULT <- "svg"
} else {
  VECTOR_FORMATS_DEFAULT <- c("pdf", "svg")
}

# Optional raster format (set to NULL to disable)
RASTER_FORMAT_DEFAULT <- NULL  # e.g., "png"
RASTER_DPI_DEFAULT <- 300


# ---------------------- 1. PACKAGES ----------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(readr)
})


# ---------------------- 1b. HELPER FUNCTIONS ----------------------
source("R/2025_blending_process/blend_figure_utils.R")
source("R/2025_blending_process/blend_evaluation_utils.R")

# ---------------------- 2. FILE LOCATIONS ----------------------
# Update these paths for your environment.
PATHS <- list(
  # Summary-model (global) performance
  summary_1965_1978      = "Monsoon_Data/results/2025_model_evaluation/summary_models_1965_1978.rds",
  summary_1965_1978_clim_mok_date = "Monsoon_Data/results/2025_model_evaluation/summary_models_1965_1978_clim_mok_date.rds",
  summary_1965_1978_no_mok_filter = "Monsoon_Data/results/2025_model_evaluation/summary_models_1965_1978_no_mok_filter.rds",
  summary_2000_2024      = "Monsoon_Data/results/2025_model_evaluation/summary_models_2000_2024.rds",
  summary_2000_2024_clim_mok_date = "Monsoon_Data/results/2025_model_evaluation/summary_models_2000_2024_clim_mok_date.rds",
  summary_2000_2024_no_mok_filter = "Monsoon_Data/results/2025_model_evaluation/summary_models_2000_2024_no_mok_filter.rds",
  # 2025 out-of-sample evaluation (produced by 2_2025_evaluation.R)
  summary_2025_mr  = "Monsoon_Data/results/2025_model_evaluation/evaluation/model_metrics_sent_vs_mr_mok_year.rds",
  
  # Yearly time-series metrics
  yearly_1965_1978 = "Monsoon_Data/results/2025_model_evaluation/yearly_metrics_global_1965_1978.rds",
  yearly_1965_1978_clim_mok_date = "Monsoon_Data/results/2025_model_evaluation/yearly_metrics_global_clim_mok_date_1965_1978.rds",
  yearly_1965_1978_no_mok_filter = "Monsoon_Data/results/2025_model_evaluation/yearly_metrics_global_no_mok_filter_1965_1978.rds",
  yearly_2000_2024 = "Monsoon_Data/results/2025_model_evaluation/yearly_metrics_global_2000_2024.rds",
  yearly_2000_2024_clim_mok_date = "Monsoon_Data/results/2025_model_evaluation/yearly_metrics_global_clim_mok_date_2000_2024.rds",
  yearly_2000_2024_no_mok_filter = "Monsoon_Data/results/2025_model_evaluation/yearly_metrics_global_no_mok_filter_2000_2024.rds",
  
  # 2025 yearly metrics (produced by 2_2025_evaluation.R)
  yearly_2025     = "Monsoon_Data/results/2025_model_evaluation/yearly_metrics_2025.rds",
  yearly_2025_clim_mok_date = "Monsoon_Data/results/2025_model_evaluation/yearly_metrics_2025_clim_mok_date.rds",
  yearly_2025_no_mok_filter = "Monsoon_Data/results/2025_model_evaluation/yearly_metrics_2025_no_mok_filter.rds",
  
  # Reliability chartdata directory
  reliability_dir = "Monsoon_Data/results/2025_model_evaluation/calibration plots",
  
  # Input for producing skill maps.
  cell_metrics_2000_2024 = "Monsoon_Data/results/2025_model_evaluation/cell_metrics_2000_2024.rds",
  
  # India boundary data used for plotting (RDS or RDA or geojson; see loader below)
  # Important to use this one as the R default's northern border doesn't match the official Indian border
  india_boundary_path = "Monsoon_Data/maps/india_boundary.csv"
)


# ---------------------- 3. LABELS & COLORS ----------------------

# Display names for models (internal_name = "Pretty Name")
MODEL_LABELS <- c(
  "clim_raw"                    = "Evolving Expectations Model",
  "unc_clim_raw"                = "Static Climatology",
  "blended_model"                 = "Blended Model",
  "ngcm_raw"                    = "NGCM (raw)",
  "aifs_calibrated"               = "AIFS (calibrated)",
  "aifs_calibrated_clim_mok_date" = "AIFS (calibrated)",
  "ngcm_calibrated"               = "NGCM (calibrated)",
  "ngcm_blend"                  = "NGCM : EE",
  "int_all"                     = "NGCM : AIFS : EE",
  "ngcm_clim_mok_date_raw"               = "NGCM (raw)",
  "ngcm_calibrated_clim_mok_date"          = "NGCM (calibrated)",
  "mme_clim_mok_date_clim_raw_opt_rps" = "Best MME",
  "mme_no_mok_filter_clim_raw_opt_rps"       = "Best MME"
)

METRIC_LABELS <- c(
  "brier_skill" = "BSS",
  "rps_skill"   = "RPSS",
  "auc"         = "AUC"
)

# Period labels used in titles / axes
PERIOD_LABELS <- c(
  "2000-2024"      = "2000-2024",
  "1965-1978"      = "1965-1978",
  "2000-2024-clim_mok_date" = "2000-2024 (Climatological MOK Date)",
  "2025_MR"        = "2025"
)

# Okabe-Ito palette (colorblind-friendly)
OKABE_ITO <- c(
  black  = "#000000", orange = "#E69F00", sky = "#56B4E9", bluish_green = "#009E73",
  yellow = "#F0E442", blue   = "#0072B2", vermillion = "#D55E00", reddish_purple = "#CC79A7"
)

MODEL_COLORS_CORE <- c(
  "unc_clim_raw"                = "#237192",   # dark blue
  "clim_raw"                    = "#70c8d4",   # light cyan
  "ngcm_calibrated_clim_mok_date" = "#eb7900",   # orange
  "ngcm_calibrated_mok"           = "#eb7900",   # orange (known MOK variant)
  "ngcm_calibrated"               = "#eb7900",   # orange (no_mok_filter variant)
  "blended_model"               = "#ca1b00"    # red
)

METRIC_COLORS <- c(
  "brier_skill" = "#eb5e71",    # coral pink
  "rps_skill"   = "#93c1dc",    # light steel blue
  "auc"         = "#2b7d87"     # dark teal
)

# Colors for the three reliability-diagram histogram panels
RELIABILITY_COLORS <- c(
  "ngcm_clim_mok_date_raw"       = "#c49a2c",   # dark gold     (panel A)
  "ngcm_raw"                     = "#c49a2c",   # dark gold     (panel A, no_mok_filter variant)
  "ngcm_calibrated_clim_mok_date"  = "#eb7900",   # orange        (panel B)
  "ngcm_calibrated"                = "#eb7900",   # orange        (panel B, no_mok_filter variant)
  "blended_model"                = "#ca1b00"    # red           (panel C)
)

MODEL_COLORS_WEEKLY <- c(
  "unc_clim_raw"                = "#237192",    # dark blue (exact from fig_3)
  "clim_raw"                    = "#70c8d4",    # light cyan (exact from fig_3)
  "ngcm_calibrated_clim_mok_date" = "#eb7900",   # orange
  "ngcm_calibrated_mok"           = "#eb7900",    # orange
  "blended_model"               = "#ca1b00"     # red (exact from fig_3)
)


# ---------------------- 4. READ DATA ----------------------

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Summary datasets
summ_1965_1978      <- safe_read_rds(PATHS$summary_1965_1978)
summ_1965_1978_clim_mok_date <- optional_read_rds(PATHS$summary_1965_1978_clim_mok_date)
summ_1965_1978_no_mok_filter <- optional_read_rds(PATHS$summary_1965_1978_no_mok_filter)
summ_2000_2024      <- safe_read_rds(PATHS$summary_2000_2024)
summ_2000_2024_clim_mok_date <- safe_read_rds(PATHS$summary_2000_2024_clim_mok_date)
summ_2000_2024_no_mok_filter <- safe_read_rds(PATHS$summary_2000_2024_no_mok_filter)

summ_2025_mr  <- optional_read_rds(PATHS$summary_2025_mr)

# Yearly series
yearly_1965_1978 <- safe_read_rds(PATHS$yearly_1965_1978)
yearly_1965_1978_clim_mok_date <- optional_read_rds(PATHS$yearly_1965_1978_clim_mok_date)
yearly_1965_1978_no_mok_filter <- optional_read_rds(PATHS$yearly_1965_1978_no_mok_filter)
yearly_2000_2024 <- safe_read_rds(PATHS$yearly_2000_2024)
yearly_2000_2024_clim_mok_date <- safe_read_rds(PATHS$yearly_2000_2024_clim_mok_date)
yearly_2000_2024_no_mok_filter <- safe_read_rds(PATHS$yearly_2000_2024_no_mok_filter)

# Append 2025 to yearly series if available (optional; produced by 2_2025_evaluation.R)
yearly_2025     <- optional_read_rds(PATHS$yearly_2025)
yearly_2025_clim_mok_date <- optional_read_rds(PATHS$yearly_2025_clim_mok_date)
yearly_2025_no_mok_filter <- optional_read_rds(PATHS$yearly_2025_no_mok_filter)

if (!is.null(yearly_2025)) {
  yearly_2000_2024 <- dplyr::bind_rows(yearly_2000_2024, yearly_2025)
}
if (!is.null(yearly_2025_clim_mok_date) && !is.null(yearly_2000_2024_clim_mok_date)) {
  yearly_2000_2024_clim_mok_date <- dplyr::bind_rows(yearly_2000_2024_clim_mok_date, yearly_2025_clim_mok_date)
}
if (!is.null(yearly_2025_no_mok_filter) && !is.null(yearly_2000_2024_no_mok_filter)) {
  yearly_2000_2024_no_mok_filter <- dplyr::bind_rows(yearly_2000_2024_no_mok_filter, yearly_2025_no_mok_filter)
}

# Minimal column checks (guardrails for replication)
required_cols <- c("model", "brier_skill", "rps_skill", "auc")
check_required_cols(summ_1965_1978, required_cols, "summ_1965_1978")
check_required_cols(summ_2000_2024, required_cols, "summ_2000_2024")
check_required_cols(summ_2000_2024_clim_mok_date, required_cols, "summ_2000_2024_clim_mok_date")
check_required_cols(summ_2000_2024_no_mok_filter, required_cols, "summ_2000_2024_no_mok_filter")
if (!is.null(summ_1965_1978_clim_mok_date)) check_required_cols(summ_1965_1978_clim_mok_date, required_cols, "summ_1965_1978_clim_mok_date")
if (!is.null(summ_1965_1978_no_mok_filter)) check_required_cols(summ_1965_1978_no_mok_filter, required_cols, "summ_1965_1978_no_mok_filter")


# ---------------------- 5. OVERALL METRICS BY PERIOD ----------------------

# Assemble the pool for the "overall metrics" figure
first_pool_parts <- list(
  add_period_tag(summ_2000_2024, "2000-2024"),
  add_period_tag(summ_1965_1978, "1965-1978")
)
if (!is.null(summ_2025_mr)) first_pool_parts <- c(first_pool_parts, list(add_period_tag(summ_2025_mr, "2025_MR")))
first_pool <- bind_rows(first_pool_parts)

# Harmonize all NGCM calibrated variants to a single canonical name.
# Specs produce multiple NGCM variants; keep the most-preferred per period,
# but preserve sole entries (e.g. 2025 only has ngcm_calibrated).
ngcm_variants <- c("ngcm_calibrated_clim_mok_date",  "ngcm_calibrated")
ngcm_pref <- c("ngcm_calibrated_clim_mok_date" = 1L,  "ngcm_calibrated" = 2L)
first_pool <- first_pool %>%
  mutate(
    .ngcm_rank = ifelse(model %in% names(ngcm_pref), ngcm_pref[model], 0L),
    model = if_else(model %in% ngcm_variants, "ngcm_cal", model)
  ) %>%
  group_by(model, period) %>%
  arrange(.ngcm_rank) %>%
  slice(1) %>%
  ungroup() %>%
  select(-.ngcm_rank)

models_first <- c("unc_clim_raw", "clim_raw", "ngcm_cal", "blended_model")
models_first <- models_first[models_first %in% unique(first_pool$model)]
period_levels_raw <- if (!is.null(summ_2025_mr)) c("2000-2024", "1965-1978", "2025_MR") else c("2000-2024", "1965-1978")
period_levels_pretty <- unname(PERIOD_LABELS[period_levels_raw])

# Add label/color for canonical NGCM name
MODEL_LABELS[["ngcm_cal"]] <- "NGCM (calibrated)"
MODEL_COLORS_CORE[["ngcm_cal"]] <- "#eb7900"

first_long <- first_pool %>%
  filter(model %in% models_first) %>%
  select(model, period, brier_skill, rps_skill, auc) %>%
  pivot_longer(cols = c(brier_skill, rps_skill, auc), names_to = "metric", values_to = "value") %>%
  mutate(
    period = factor(period, levels = period_levels_raw),
    period_display = factor(recode(as.character(period), !!!PERIOD_LABELS), levels = period_levels_pretty),
    model = factor(model, levels = models_first),
    model_display = recode(as.character(model), !!!MODEL_LABELS, .default = as.character(model)),
    metric_display = factor(recode(metric, !!!METRIC_LABELS))
  )

auc_label <- METRIC_LABELS[["auc"]]

# Shared theme for overall metrics panels
.overall_theme <- function() {
  theme_minimal(base_size = 12, base_family = "Arial") +
    theme(
      text               = element_text(colour = "#000000"),
      legend.position    = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#ebebeb", linewidth = 0.25),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(fill = NA, colour = "#1a1a1a", linewidth = 0.5),
      axis.text          = element_text(size = 10, colour = "#000000"),
      axis.title         = element_text(size = 11, colour = "#000000")
    )
}

# Explicit metric ordering: Brier → RPS → AUC (top to bottom)
metric_order <- c(METRIC_LABELS[["brier_skill"]], METRIC_LABELS[["rps_skill"]], METRIC_LABELS[["auc"]])

# Build one panel per metric in Brier → RPS → AUC order
overall_panels <- lapply(metric_order, function(ml) {
  d <- first_long %>% filter(metric_display == ml)
  p <- ggplot(d, aes(x = period_display, y = value, fill = model)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    scale_fill_manual(
      values = MODEL_COLORS_CORE[models_first],
      limits = models_first,
      labels = unname(MODEL_LABELS[models_first]),
      name = NULL
    ) +
    labs(x = NULL, y = ml, title = NULL) +
    .overall_theme()
  # AUC panel: zoom to start at 0.5
  if (ml == auc_label) {
    p <- p + coord_cartesian(ylim = c(0.5, NA))
  }
  p
})

# Legend inside top-left of first panel (Brier)
overall_panels[[1]] <- overall_panels[[1]] +
  theme(
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.text          = element_text(size = 8),
    legend.key.size      = unit(0.35, "cm"),
    legend.margin        = margin(2, 4, 2, 4)
  ) +
  guides(fill = guide_legend(nrow = 2))

p_overall <- patchwork::wrap_plots(overall_panels, ncol = 1)

save_plot(
  file_stem = "overall_metrics",
  plot = p_overall,
  out_dir = OUTPUT_DIR,
  width = 7,
  height = 6,
  vector_formats = VECTOR_FORMATS_DEFAULT,
  raster_format = RASTER_FORMAT_DEFAULT,
  raster_dpi = RASTER_DPI_DEFAULT
)

# ---- Alternate: horizontal 3-panel layout ----
# Metric label inside top-center of each panel using annotate (won't clip)
overall_panels_h <- lapply(metric_order, function(ml) {
  d <- first_long %>% filter(metric_display == ml)
  y_max <- max(d$value, na.rm = TRUE)
  p <- ggplot(d, aes(x = period_display, y = value, fill = model)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.55) +
    annotate("text", x = mean(seq_along(levels(d$period_display))),
             y = y_max * 1.12, label = ml,
             fontface = "bold", size = 3.8, colour = "#000000", family = "Arial") +
    scale_fill_manual(
      values = MODEL_COLORS_CORE[models_first],
      limits = models_first,
      labels = unname(MODEL_LABELS[models_first]),
      name = NULL
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(x = NULL, y = NULL, title = NULL) +
    .overall_theme() +
    theme(axis.text.x = element_text(size = 8, colour = "#000000"))
  # AUC panel: start at 0.5
  if (ml == auc_label) {
    p <- p + coord_cartesian(ylim = c(0.5, y_max * 1.15))
  }
  p
})

# Legend inside first panel (Brier)
overall_panels_h[[1]] <- overall_panels_h[[1]] +
  theme(
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.text          = element_text(size = 7),
    legend.key.size      = unit(0.3, "cm"),
    legend.margin        = margin(1, 3, 1, 3)
  ) +
  guides(fill = guide_legend(nrow = 2))

p_overall_h <- patchwork::wrap_plots(overall_panels_h, ncol = 3)

save_plot(
  file_stem = "overall_metrics_horizontal",
  plot = p_overall_h,
  out_dir = OUTPUT_DIR,
  width = 14,
  height = 4.5,
  vector_formats = VECTOR_FORMATS_DEFAULT,
  raster_format = RASTER_FORMAT_DEFAULT,
  raster_dpi = RASTER_DPI_DEFAULT
)

# ---- Helper: build overall_metrics for a given variant ----
build_overall_variant <- function(summ_main, summ_hist, summ_2025, variant_suffix, ngcm_variants_extra = NULL,
                                  ngcm_pref_order = c("ngcm_calibrated_clim_mok_date" = 1L,"ngcm_calibrated" = 2L)) {
  pool_parts <- list(
    add_period_tag(summ_main, "2000-2024"),
    add_period_tag(summ_hist, "1965-1978")
  )
  if (!is.null(summ_2025)) pool_parts <- c(pool_parts, list(add_period_tag(summ_2025, "2025_MR")))
  pool <- bind_rows(pool_parts)
  
  # Harmonize NGCM variants — keep most-preferred per period, preserve sole entries
  all_ngcm <- c("ngcm_calibrated_clim_mok_date", "ngcm_calibrated", ngcm_variants_extra)
  pref_lookup <- ngcm_pref_order[intersect(names(ngcm_pref_order), all_ngcm)]
  pool <- pool %>%
    mutate(
      .ngcm_rank = ifelse(model %in% names(pref_lookup), pref_lookup[model], 0L),
      model = if_else(model %in% all_ngcm, "ngcm_cal", model)
    ) %>%
    group_by(model, period) %>%
    arrange(.ngcm_rank) %>%
    slice(1) %>%
    ungroup() %>%
    select(-.ngcm_rank)
  
  mf <- c("unc_clim_raw", "clim_raw", "ngcm_cal", "blended_model")
  mf <- mf[mf %in% unique(pool$model)]
  pl_raw <- if (!is.null(summ_2025)) c("2000-2024", "1965-1978", "2025_MR") else c("2000-2024", "1965-1978")
  pl_pretty <- unname(PERIOD_LABELS[pl_raw])
  
  vlong <- pool %>%
    filter(model %in% mf) %>%
    select(model, period, brier_skill, rps_skill, auc) %>%
    pivot_longer(cols = c(brier_skill, rps_skill, auc), names_to = "metric", values_to = "value") %>%
    mutate(
      period = factor(period, levels = pl_raw),
      period_display = factor(recode(as.character(period), !!!PERIOD_LABELS), levels = pl_pretty),
      model = factor(model, levels = mf),
      metric_display = factor(recode(metric, !!!METRIC_LABELS))
    )
  
  ml <- metric_order  # Brier → RPS → AUC
  
  # Vertical
  panels_v <- lapply(ml, function(m) {
    d <- vlong %>% filter(metric_display == m)
    p <- ggplot(d, aes(x = period_display, y = value, fill = model)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.65) +
      scale_fill_manual(values = MODEL_COLORS_CORE[mf], limits = mf, labels = unname(MODEL_LABELS[mf]), name = NULL) +
      labs(x = NULL, y = m, title = NULL) + .overall_theme()
    if (m == auc_label) p <- p + coord_cartesian(ylim = c(0.5, NA))
    p
  })
  panels_v[[1]] <- panels_v[[1]] +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.text = element_text(size = 8), legend.key.size = unit(0.35, "cm"),
          legend.margin = margin(2, 4, 2, 4)) +
    guides(fill = guide_legend(nrow = 2))
  p_v <- patchwork::wrap_plots(panels_v, ncol = 1)
  save_plot(paste0("overall_metrics", variant_suffix), p_v, OUTPUT_DIR, 7, 6, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
  
  # Horizontal
  panels_h <- lapply(ml, function(m) {
    d <- vlong %>% filter(metric_display == m)
    ym <- max(d$value, na.rm = TRUE)
    p <- ggplot(d, aes(x = period_display, y = value, fill = model)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.55) +
      annotate("text", x = mean(seq_along(levels(d$period_display))), y = ym * 1.12, label = m,
               fontface = "bold", size = 3.8, colour = "#000000", family = "Arial") +
      scale_fill_manual(values = MODEL_COLORS_CORE[mf], limits = mf, labels = unname(MODEL_LABELS[mf]), name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
      labs(x = NULL, y = NULL, title = NULL) + .overall_theme() + theme(axis.text.x = element_text(size = 8))
    if (m == auc_label) p <- p + coord_cartesian(ylim = c(0.5, ym * 1.15))
    p
  })
  panels_h[[1]] <- panels_h[[1]] +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.text = element_text(size = 7), legend.key.size = unit(0.3, "cm"),
          legend.margin = margin(1, 3, 1, 3)) +
    guides(fill = guide_legend(nrow = 2))
  p_h <- patchwork::wrap_plots(panels_h, ncol = 3)
  save_plot(paste0("overall_metrics_horizontal", variant_suffix), p_h, OUTPUT_DIR, 14, 4.5, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
}

# clim_mok_date variant
if (!is.null(summ_2000_2024_clim_mok_date)) {
  build_overall_variant(
    summ_2000_2024_clim_mok_date,
    if (!is.null(summ_1965_1978_clim_mok_date)) summ_1965_1978_clim_mok_date else summ_1965_1978,
    summ_2025_mr,
    "_clim_mok_date"
  )
}

# no_mok_filter variant
if (!is.null(summ_2000_2024_no_mok_filter)) {
  build_overall_variant(
    summ_2000_2024_no_mok_filter,
    if (!is.null(summ_1965_1978_no_mok_filter)) summ_1965_1978_no_mok_filter else summ_1965_1978,
    summ_2025_mr,
    "_no_mok_filter",
    ngcm_pref_order = c("ngcm_calibrated" = 1L)
  )
}

# ---------------------- 6. CLIMATOLOGY TRAINING WINDOW VARIATION ----------------------

# Identify versions of clim_raw / blended_model that exist in BOTH 1965-1978 and 2000-2024
models_in_both <- intersect(unique(summ_2000_2024$model), unique(summ_1965_1978$model))
models_keep <- models_in_both[str_detect(models_in_both, "^(clim_raw|blended_model)(_|$)")]

model_map <- tibble(model = models_keep) %>%
  mutate(
    model_type = case_when(
      str_detect(model, "^clim_raw")    ~ "Evolving Expectations Model",
      str_detect(model, "^blended_model") ~ "Blended Model",
      TRUE                               ~ "Other"
    ),
    clim_train = vapply(model, infer_train_window, character(1))
  ) %>%
  filter(!is.na(clim_train))

if (nrow(model_map) > 0) {
  
  train_levels <- sort(unique(model_map$clim_train))
  train_cols <- grDevices::hcl.colors(length(train_levels), palette = "Batlow", alpha = 0.95)
  clim_train_colors <- setNames(train_cols, train_levels)
  
  auc_scale <- 4
  auc_min <- summ_2000_2024 %>% filter(model == "unc_clim_raw") %>% pull(auc)
  
  six_long <- bind_rows(
    summ_1965_1978 %>% mutate(test_period = "1965-1978"),
    summ_2000_2024 %>% mutate(test_period = "2000-2024")
  ) %>%
    filter(model %in% model_map$model) %>%
    left_join(model_map, by = "model") %>%
    select(test_period, model_type, clim_train, brier_skill, rps_skill, auc) %>%
    pivot_longer(cols = c(brier_skill, rps_skill, auc), names_to = "metric", values_to = "value") %>%
    mutate(
      test_period = factor(test_period, levels = c("2000-2024", "1965-1978")),
      model_type  = factor(model_type,  levels = c("Evolving Expectations Model", "Blended Model")),
      clim_train  = factor(clim_train,  levels = train_levels),
      metric_display = recode(metric, !!!METRIC_LABELS),
      value_plot = if_else(metric == "auc", (value - auc_min) * auc_scale, value)
    )
  
  bars_df <- six_long %>% filter(test_period == "2000-2024")
  pts_df  <- six_long %>% filter(test_period == "1965-1978")
  
  # Shared theme for climatology panels
  .clim_theme <- function() {
    theme_minimal(base_size = 13, base_family = "Arial") +
      theme(
        text                 = element_text(colour = "#000000"),
        panel.grid.major.x   = element_blank(),
        panel.grid.major.y   = element_line(colour = "#ebebeb", linewidth = 0.25),
        panel.grid.minor     = element_blank(),
        panel.border         = element_rect(fill = NA, colour = "#1a1a1a", linewidth = 0.5),
        axis.text            = element_text(size = 11, colour = "#000000"),
        axis.title           = element_text(size = 12, colour = "#000000"),
        legend.position      = "none",
        plot.margin          = margin(5, 5, 5, 5)
      )
  }
  
  # Build one panel per metric using annotate() for corner labels
  # After coord_flip():
  #   x=-Inf → visual BOTTOM,  x=Inf → visual TOP
  #   y=-Inf → visual LEFT,    y=Inf → visual RIGHT
  clim_metrics <- unique(bars_df$metric_display)
  clim_panels <- lapply(clim_metrics, function(ml) {
    bd <- bars_df %>% dplyr::filter(metric_display == ml)
    pd <- pts_df  %>% dplyr::filter(metric_display == ml)
    # Compute horizontal midpoint for centering the title
    # Include 0 since geom_col anchors bars at 0, so the axis always includes it
    all_vals <- c(0, bd$value_plot, pd$value_plot)
    y_mid <- mean(range(all_vals, na.rm = TRUE))
    
    p <- ggplot(bd, aes(x = model_type, y = value_plot, fill = clim_train)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.72) +
      geom_point(
        data = pd,
        aes(x = model_type, y = value_plot, fill = clim_train),
        position = position_dodge(width = 0.8),
        shape = 1, size = 2.5, color = "grey25", stroke = 0.8,
        inherit.aes = FALSE, show.legend = FALSE
      ) +
      # METRIC NAME: top-center of box
      annotate("text", x = Inf, y = y_mid, label = ml,
               hjust = 0.5, vjust = 1.3,
               size = 5, fontface = "bold", colour = "#000000", family = "Arial") +
      scale_fill_manual(values = clim_train_colors, name = "Climatology training period") +
      scale_y_continuous(name = NULL) +
      labs(x = NULL, y = NULL) +
      .clim_theme() +
      coord_flip(clip = "off")
    
    p
  })
  
  # AUC panel: add bar/circle key in bottom-right corner
  auc_idx <- which(clim_metrics == METRIC_LABELS[["auc"]])
  if (length(auc_idx) > 0) {
    clim_panels[[auc_idx]] <- clim_panels[[auc_idx]] +
      annotate("text", x = -Inf, y = Inf,
               label = "bar = 2000-2024\n\u25cb = 1965-1978",
               hjust = 1.05, vjust = -0.3,
               size = 3.2, colour = "#000000", family = "Arial", lineheight = 0.9)
  }
  
  # Bottom panel: show legend
  clim_panels[[length(clim_panels)]] <- clim_panels[[length(clim_panels)]] +
    theme(
      legend.position      = "bottom",
      legend.text          = element_text(size = 9),
      legend.title         = element_text(size = 10),
      legend.key.size      = unit(0.4, "cm"),
      legend.key.width     = unit(0.6, "cm")
    ) +
    guides(fill = guide_legend(nrow = 2))
  
  p_train_window <- patchwork::wrap_plots(clim_panels, ncol = 1)
  
  save_plot(
    file_stem = "climatology_training_period_variation",
    plot = p_train_window,
    out_dir = OUTPUT_DIR,
    width = 12,
    height = 9,
    vector_formats = VECTOR_FORMATS_DEFAULT,
    raster_format = RASTER_FORMAT_DEFAULT,
    raster_dpi = RASTER_DPI_DEFAULT
  )
  
} else {
  message("Skipping climatology training window plot: no window variants found in both periods.")
}


# ---------------------- 7. MODEL COMPARISONS (FIG 4 STYLE — ALL 6 VARIANTS) ----------------------

# Model sets for each variant family
models_second_clim_mok_date_bars <- c(
  "clim_raw",
  "ngcm_clim_mok_date_raw",
  "aifs_calibrated_clim_mok_date",
  "ngcm_calibrated_clim_mok_date",
  "ngcm_blend",
  "int_all",
  "mme_clim_mok_date_clim_raw_opt_rps",
  "blended_model"
)

models_second_no_mok_filter <- c(
  "clim_raw",
  "ngcm_raw",
  "aifs_calibrated",
  "ngcm_calibrated",
  "ngcm_blend",
  "int_all",
  "mme_no_mok_filter_clim_raw_opt_rps",
  "blended_model"
)



# --- 2000-2024 variants ---
p_skill_2000_2024 <- make_fig4_variant(
  summ_2000_2024, models_second_clim_mok_date_bars,
  "2000-2024", "skill_comparison_2000_2024", OUTPUT_DIR
)

p_skill_2000_2024_clim_mok_date <- make_fig4_variant(
  summ_2000_2024_clim_mok_date, models_second_clim_mok_date_bars,
  "2000-2024, Climatological MOK Filter", "skill_comparison_2000_2024_clim_mok_date_filter", OUTPUT_DIR
)

p_skill_2000_2024_no_mok_filter <- make_fig4_variant(
  summ_2000_2024_no_mok_filter, models_second_no_mok_filter,
  "2000-2024, No MOK Filter", "skill_comparison_2000_2024_no_mok_filter", OUTPUT_DIR
)

# --- 1965-1978 variants ---
p_skill_1965_1978 <- make_fig4_variant(
  summ_1965_1978, models_second_clim_mok_date_bars,
  "1965-1978", "skill_comparison_1965_1978", OUTPUT_DIR
)

p_skill_1965_1978_clim_mok_date <- make_fig4_variant(
  summ_1965_1978_clim_mok_date, models_second_clim_mok_date_bars,
  "1965-1978, Climatological MOK Filter", "skill_comparison_1965_1978_clim_mok_date_filter", OUTPUT_DIR
)

p_skill_1965_1978_no_mok_filter <- make_fig4_variant(
  summ_1965_1978_no_mok_filter, models_second_no_mok_filter,
  "1965-1978, No MOK Filter", "skill_comparison_1965_1978_no_mok_filter", OUTPUT_DIR
)


# ---------------------- 8. WEEKLY PERFORMANCE ----------------------

p_weekly_1965_1978 <- make_weekly_bins_plot(
  df = summ_1965_1978,
  title_suffix = "1965-1978",
  variant = "standard",
  include_later = FALSE,
  model_labels = MODEL_LABELS,
  metric_labels = METRIC_LABELS,
  model_colors = MODEL_COLORS_WEEKLY,
  font_size = FONT_SIZE
)

p_weekly_2000_2024 <- make_weekly_bins_plot(
  df = summ_2000_2024,
  title_suffix = "2000-2024",
  variant = "standard",
  include_later = FALSE,
  model_labels = MODEL_LABELS,
  metric_labels = METRIC_LABELS,
  model_colors = MODEL_COLORS_WEEKLY,
  font_size = FONT_SIZE
)

p_weekly_2000_2024_clim_mok_date <- make_weekly_bins_plot(
  df = summ_2000_2024_clim_mok_date,
  title_suffix = "2000-2024 (Climatological MOK Filter)",
  variant = "clim_mok_date",
  include_later = FALSE,
  model_labels = MODEL_LABELS,
  metric_labels = METRIC_LABELS,
  model_colors = MODEL_COLORS_WEEKLY,
  font_size = FONT_SIZE
)

p_weekly_2000_2024_no_mok_filter <- make_weekly_bins_plot(
  df = summ_2000_2024_no_mok_filter,
  title_suffix = "2000-2024 (No MOK Filter)",
  variant = "no_mok_filter",
  include_later = FALSE,
  model_labels = MODEL_LABELS,
  metric_labels = METRIC_LABELS,
  model_colors = MODEL_COLORS_WEEKLY,
  font_size = FONT_SIZE
)

save_plot("skill_by_week_1965_1978", p_weekly_1965_1978, OUTPUT_DIR, 12, 7.5, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
save_plot("skill_by_week_2000_2024", p_weekly_2000_2024, OUTPUT_DIR, 12, 7.5, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
save_plot("skill_by_week_2000_2024_clim_mok_date_filter", p_weekly_2000_2024_clim_mok_date, OUTPUT_DIR, 12, 7.5, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
save_plot("skill_by_week_2000_2024_no_mok_filter", p_weekly_2000_2024_no_mok_filter, OUTPUT_DIR, 12, 7.5, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)


# ---------------------- 9. YEARLY PERFORMANCE ----------------------

p_yearly_1965_1978 <- make_yearly_plot(
  yearly_df = yearly_1965_1978,
  title_suffix = "1965-1978",
  model_labels = MODEL_LABELS,
  metric_labels = METRIC_LABELS,
  okabe_ito = OKABE_ITO,
  font_size = FONT_SIZE
)

p_yearly_2000_2024 <- make_yearly_plot(
  yearly_df = yearly_2000_2024,
  title_suffix = if (2025L %in% yearly_2000_2024$year) "2000-2025" else "2000-2024",
  model_labels = MODEL_LABELS,
  metric_labels = METRIC_LABELS,
  okabe_ito = OKABE_ITO,
  font_size = FONT_SIZE
)

p_yearly_2000_2024_clim_mok_date <- if (!is.null(yearly_2000_2024_clim_mok_date)) {
  make_yearly_plot(
    yearly_df = yearly_2000_2024_clim_mok_date,
    title_suffix = if (2025L %in% yearly_2000_2024_clim_mok_date$year) "2000-2025 (Climatological MOK Filter)" else "2000-2024 (Climatological MOK Filter)",
    model_labels = MODEL_LABELS,
    metric_labels = METRIC_LABELS,
    okabe_ito = OKABE_ITO,
    font_size = FONT_SIZE
  )
}

p_yearly_2000_2024_no_mok_filter <- if (!is.null(yearly_2000_2024_no_mok_filter)) {
  make_yearly_plot(
    yearly_df = yearly_2000_2024_no_mok_filter,
    title_suffix = if (2025L %in% yearly_2000_2024_no_mok_filter$year) "2000-2025 (No MOK Filter)" else "2000-2024 (No MOK Filter)",
    model_labels = MODEL_LABELS,
    metric_labels = METRIC_LABELS,
    okabe_ito = OKABE_ITO,
    font_size = FONT_SIZE,
    ngcm_model = "ngcm_calibrated"
  )
}

save_plot("yearly_metrics_by_model_1965_1978", p_yearly_1965_1978, OUTPUT_DIR, 12, 8, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
save_plot("yearly_metrics_by_model_2000_2024", p_yearly_2000_2024, OUTPUT_DIR, 12, 8, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
if (!is.null(p_yearly_2000_2024_clim_mok_date)) {
  save_plot("yearly_metrics_by_model_2000_2024_clim_mok_date", p_yearly_2000_2024_clim_mok_date, OUTPUT_DIR, 12, 8, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
}
if (!is.null(p_yearly_2000_2024_no_mok_filter)) {
  save_plot("yearly_metrics_by_model_2000_2024_no_mok_filter", p_yearly_2000_2024_no_mok_filter, OUTPUT_DIR, 12, 8, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
}


# ---------------------- 9b. COMBINED FIG 2 (WEEKLY + YEARLY) ----------------------

# Helper to generate a fig_3 combined plot for a given variant
make_fig2 <- function(summ_df, yearly_df, variant, file_stem, model_labels) {
  if (is.null(summ_df) || is.null(yearly_df)) return(NULL)
  tryCatch({
    p <- plot_fig2_combined(
      summ_df      = summ_df,
      yearly_df    = yearly_df,
      variant      = variant,
      model_labels = model_labels,
      model_colors = FIG2_MODEL_COLORS
    )
    
    save_plot(file_stem, p, OUTPUT_DIR, 6.5, 5.8,
              VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
    p
  }, error = function(e) { message("Skipping fig_3 (", file_stem, "): ", e$message); NULL })
}

p_fig2_2000_2024 <- make_fig2(
  summ_2000_2024, yearly_2000_2024,
  "standard", "fig2_combined_2000_2024", MODEL_LABELS
)

p_fig2_2000_2024_clim_mok <- make_fig2(
  summ_2000_2024_clim_mok_date, yearly_2000_2024_clim_mok_date,
  "clim_mok_date", "fig2_combined_2000_2024_clim_mok_date", MODEL_LABELS
)

p_fig2_2000_2024_no_mok <- make_fig2(
  summ_2000_2024_no_mok_filter, yearly_2000_2024_no_mok_filter,
  "no_mok_filter", "fig2_combined_2000_2024_no_mok_filter", MODEL_LABELS
)

p_fig2_1965_1978 <- make_fig2(
  summ_1965_1978, yearly_1965_1978,
  "standard", "fig2_combined_1965_1978", MODEL_LABELS
)

p_fig2_1965_1978_clim_mok <- make_fig2(
  summ_1965_1978_clim_mok_date, yearly_1965_1978_clim_mok_date,
  "clim_mok_date", "fig2_combined_1965_1978_clim_mok_date", MODEL_LABELS
)

p_fig2_1965_1978_no_mok <- make_fig2(
  summ_1965_1978_no_mok_filter, yearly_1965_1978_no_mok_filter,
  "no_mok_filter", "fig2_combined_1965_1978_no_mok_filter", MODEL_LABELS
)


# ---------------------- 10. RELIABILITY DIAGRAMS (FIG 3 STYLE — ALL 6 VARIANTS) ----------------------

rel_out_dir <- file.path(OUTPUT_DIR, "reliability")
dir.create(rel_out_dir, showWarnings = FALSE, recursive = TRUE)

# Define the 6 dataset configurations
# Each entry: period_tag (for file lookup), label (for output stem), and which
# model triplet to use (standard vs no_mok_filter model names).
reliability_configs <- list(
  list(
    period_tag  = "2000_2024",
    label       = "2000_2024",
    models      = c("ngcm_clim_mok_date_raw", "ngcm_calibrated_clim_mok_date", "blended_model")
  ),
  list(
    period_tag  = "clim_mok_date_2000_2024",
    label       = "2000_2024_clim_mok_date",
    models      = c("ngcm_clim_mok_date_raw", "ngcm_calibrated_clim_mok_date", "blended_model")
  ),
  list(
    period_tag  = "no_mok_filter_2000_2024",
    label       = "2000_2024_no_mok_filter",
    models      = c("ngcm_raw", "ngcm_calibrated", "blended_model")
  ),
  list(
    period_tag  = "1965_1978",
    label       = "1965_1978",
    models      = c("ngcm_clim_mok_date_raw", "ngcm_calibrated_clim_mok_date", "blended_model")
  ),
  list(
    period_tag  = "clim_mok_date_1965_1978",
    label       = "1965_1978_clim_mok_date",
    models      = c("ngcm_clim_mok_date_raw", "ngcm_calibrated_clim_mok_date", "blended_model")
  ),
  list(
    period_tag  = "no_mok_filter_1965_1978",
    label       = "1965_1978_no_mok_filter",
    models      = c("ngcm_raw", "ngcm_calibrated", "blended_model")
  )
)

for (cfg in reliability_configs) {
  # Build file-lookup grid for this variant
  rel_grid <- tibble::tibble(model = cfg$models) %>%
    dplyr::mutate(
      period_tag   = cfg$period_tag,
      path         = purrr::map_chr(model, ~ make_chartdata_path(PATHS$reliability_dir, .x, cfg$period_tag, bins_tag = "10bins")),
      model_pretty = pretty_model_name(model, MODEL_LABELS)
    ) %>%
    dplyr::filter(file.exists(path))
  
  if (nrow(rel_grid) == 0) {
    message("Skipping reliability 3-panel for ", cfg$label, ": no chart-data files found.")
    next
  }
  
  # Read all chart-data RDS files for this variant
  rel_data <- purrr::pmap_dfr(
    rel_grid %>% dplyr::select(path, model, model_pretty),
    function(path, model, model_pretty) {
      read_chartdata_one(path, model, model_pretty, period_pretty = cfg$label)
    }
  )
  
  # Produce the 3-panel reliability diagram
  p_rel3 <- plot_reliability_3panel(
    rel_data           = rel_data,
    model_order        = cfg$models,
    model_labels       = MODEL_LABELS,
    reliability_colors = RELIABILITY_COLORS,
    font_size          = FONT_SIZE
  )
  
  save_plot(
    file_stem      = paste0("reliability_3panel_", cfg$label),
    plot           = p_rel3,
    out_dir        = rel_out_dir,
    width          = 7.2,
    height         = 2.36,
    vector_formats = VECTOR_FORMATS_DEFAULT,
    raster_format  = RASTER_FORMAT_DEFAULT,
    raster_dpi     = RASTER_DPI_DEFAULT
  )
}

# ---------------------- 11. MAPS ----------------------

maps_out_dir <- file.path(OUTPUT_DIR, "maps_2000_2024")
dir.create(maps_out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read inputs ----

cell_2000_2024 <- safe_read_rds(PATHS$cell_metrics_2000_2024)

# Expectation: this is per-cell model output with (at least) lat, lon, model and metrics:
# brier, rps, auc (or already precomputed brier_skill/rps_skill/auc_diff).
# If your file uses different names, adjust in the "metric preparation" block below.
check_required_cols(cell_2000_2024, c("id", "lat", "lon", "model"), "cell_2000_2024")

# Remove the ALL summary row (id == "ALL", lat/lon are NA_character_) and ensure
# lat/lon are numeric for spatial operations (polygon building, maps).
cell_2000_2024 <- cell_2000_2024 %>%
  dplyr::filter(id != "ALL") %>%
  dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon))

# India boundary: support either sf object or a data.frame with long/lat/group
india_map <- read_india_boundary(PATHS$india_boundary_path)

# Filter to official India boundary (id == 253)
if ("id" %in% names(india_map)) {
  india_map <- india_map %>% dplyr::filter(id == 253)
}

# ---- Build polygons once ----

# Grid centers are unique id/lat/lon in the cell file; allowed_cells can be the same set.
grid_centers <- cell_2000_2024 %>% dplyr::distinct(id, lat, lon)
allowed_cells <- grid_centers

poly_obj <- build_polygons_for_mapping(grid_centers = grid_centers, allowed_cells = allowed_cells)
polygons_df <- poly_obj$polygons_df

# ---- Read dissemination cells (thicker border overlay) ----
dissem_cells_path <- "Monsoon_Data/dissemination_cells.csv"
dissem_poly_df <- NULL
if (file.exists(dissem_cells_path)) {
  dissem_raw <- readr::read_csv(dissem_cells_path, show_col_types = FALSE)
  # Match dissemination cells to grid centers by lat/lon
  dissem_matched <- grid_centers %>%
    dplyr::inner_join(dissem_raw, by = c("lat", "lon"))
  if (nrow(dissem_matched) > 0) {
    dissem_poly_obj <- build_polygons_for_mapping(
      grid_centers = dissem_matched,
      allowed_cells = dissem_matched
    )
    dissem_poly_df <- dissem_poly_obj$polygons_df
    message("Dissemination cells matched: ", dplyr::n_distinct(dissem_matched$id), " cells")
  }
}

# ---- Metric preparation (2000-2024 only) ----
# We will map "final vs clim baseline" per-cell skill/diffs.
# Uses your helper summarize_maps_compare(all_cells, method, clim_model, final_model).

# Pick which CV method to use if you have multiple; if you only have one, this is ignored.
# Common possibilities: cv_method column exists; otherwise we treat everything as one method.
method_levels <- if ("cv_method" %in% names(cell_2000_2024)) unique(cell_2000_2024$cv_method) else "default"

# Choose baseline + final models (change here if needed)
CLIM_MODEL  <- "clim_raw"
FINAL_MODEL <- "blended_model"

# Create per-cell comparison table (one row per lat/lon for the chosen method)
cell_comp_2000_2024 <- if ("cv_method" %in% names(cell_2000_2024)) {
  cell_2000_2024 %>%
    dplyr::filter(.data$cv_method %in% method_levels[1]) %>%  # take first method by default
    summarize_maps_compare(method = method_levels[1], clim_model = CLIM_MODEL, final_model = FINAL_MODEL)
} else {
  cell_2000_2024 %>%
    summarize_maps_compare(method = "default", clim_model = CLIM_MODEL, final_model = FINAL_MODEL)
}

# id already present from summarize_maps_compare output

# ---- Plot settings ----

# What to map
map_specs <- tibble::tribble(
  ~metric_col,   ~title_text,                                      ~file_stem,
  "brier_skill", "Brier Skill (Blended Model vs Evolving Expectations Model)", "map_brier_skill_2000_2024",
  "rps_skill",   "RPS Skill (Blended Model vs Evolving Expectations Model)",   "map_rps_skill_2000_2024",
  "auc_diff",    "AUC Difference (Blended Model - Evolving Expectations Model)","map_auc_diff_2000_2024"
)

# ---- Make and save maps ----
plots <- purrr::pmap(
  map_specs,
  function(metric_col, title_text, file_stem) {
    leg_title <- legend_title_for_metric(metric_col)
    
    p <- plot_metric_map(
      cell_metrics   = cell_comp_2000_2024,
      polygons_df    = polygons_df,
      metric_col     = metric_col,
      title_text     = title_text,
      india_map      = india_map,
      legend_title   = leg_title,
      dissem_poly_df = dissem_poly_df
    ) +
      ggplot2::theme_minimal(base_size = FONT_SIZE) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "right"
      )
    
    save_plot(
      file_stem = file_stem,
      plot = p,
      out_dir = maps_out_dir,
      width = 10,
      height = 10,
      vector_formats = VECTOR_FORMATS_DEFAULT,
      raster_format = RASTER_FORMAT_DEFAULT,
      raster_dpi = RASTER_DPI_DEFAULT
    )
    
    p
  }
)
