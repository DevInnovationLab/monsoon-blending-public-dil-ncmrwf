# ==============================================================================
# File: 3b_quick_figures.R
# ==============================================================================
# Purpose
#   Generic diagnostic figure script for a single run of script 1.
#   Auto-discovers all models present in the output data and produces four
#   types of plots: model comparison bars, weekly performance, yearly time
#   series, and metric maps.
#
# Inputs
#   - specs/2025_blend/<spec_id>.yml (same spec used by script 1)
#   - Monsoon_Data/results/2025_model_evaluation/summary_models_*.rds
#   - Monsoon_Data/results/2025_model_evaluation/yearly_metrics_global*.rds
#   - Monsoon_Data/results/2025_model_evaluation/cell_metrics_*.rds
#   - Monsoon_Data/results/2025_model_evaluation/map_inputs*.rds
#
# Outputs
#   - figures/<spec_id>/*.pdf, *.svg
#
# Usage
#   Rscript pipelines/2025_blending_process/3b_quick_figures.R --spec_id cv_models
#
# Dependencies
#   - R/2025_blending_process/blend_figure_utils.R
#   - R/2025_blending_process/blend_evaluation_utils.R
# ==============================================================================

# ---- 1. Packages + source utils ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(readr)
})

source("R/2025_blending_process/blend_figure_utils.R")
source("R/2025_blending_process/blend_evaluation_utils.R")

# ---- 2. Parse CLI args & read spec ----
if (!interactive()) {
  option_list <- list(
    optparse::make_option("--spec_id", type = "character", default = "cv_models",
                          help = "Spec file name (without .yml) in specs/2025_blend/")
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  spec_id <- opt$spec_id
} else {
  spec_id <- "cv_models"
}

spec_path <- file.path("specs", "2025_blend", paste0(spec_id, ".yml"))
spec <- yaml::read_yaml(spec_path)

cutoff_mode <- spec$run$cutoff_mode %||% "mok"
MR <- isTRUE(spec$run$MR %||% TRUE)
true_holdout_years <- as.integer(unlist(spec$run$true_holdout_years %||% integer(0)))
cv_holdout_years   <- as.integer(unlist(spec$run$cv_holdout_years %||% integer(0)))
holdout_years <- c(true_holdout_years, cv_holdout_years)

cutoff_tag <- make_cutoff_tag(cutoff_mode)
year_tag   <- make_year_tag(holdout_years)
output_tag <- paste0(cutoff_tag, year_tag)

path_box <- "Monsoon_Data/results/2025_model_evaluation"

# ---- Construct expected filenames (same logic as script 1) ----
if (length(holdout_years) == 1) {
  year_range <- as.character(holdout_years)
} else {
  year_range <- paste0(min(holdout_years), "_", max(holdout_years))
}

summary_file <- if (MR) {
  file.path(path_box, paste0("summary_models_", year_range, cutoff_tag, ".rds"))
} else {
  file.path(path_box, "summary_models_q.rds")
}
yearly_file     <- file.path(path_box, paste0("yearly_metrics_global", output_tag, ".rds"))
cell_file       <- file.path(path_box, paste0("cell_metrics_", year_range, cutoff_tag, ".rds"))
map_inputs_file <- file.path(path_box, paste0("map_inputs", output_tag, ".rds"))

# ---- Output directory ----
OUTPUT_DIR <- file.path("figures", spec_id)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Vector format defaults
HAS_CAIRO <- capabilities("cairo")
if (!HAS_CAIRO) {
  message("Note: Cairo not available - PDF output disabled, producing SVG only.")
  VECTOR_FORMATS <- "svg"
} else {
  VECTOR_FORMATS <- c("pdf", "svg")
}

# ---- 3. Read RDS files ----
cat("Looking for summary file:", summary_file, "\n")
summ <- safe_read_rds(summary_file)

yearly <- optional_read_rds(yearly_file)
if (is.null(yearly)) message("Yearly metrics file not found: ", yearly_file, " - skipping yearly plot.")

cell_metrics <- optional_read_rds(cell_file)
map_inputs   <- optional_read_rds(map_inputs_file)

# ---- 4. Auto-discover models, build palette + labels ----
all_models <- unique(summ$model)
# Exclude unc_clim_raw (static baseline) from bar models
bar_models <- setdiff(all_models, "unc_clim_raw")

cat("Discovered models:", paste(all_models, collapse = ", "), "\n")

# Known labels from 3_produce_figures.R
MODEL_LABELS <- c(
  "clim_raw"                        = "Evolving Expectations Model",
  "unc_clim_raw"                    = "Static Climatology",
  "blended_model"                   = "Blended Model",
  "ngcm_raw"                        = "NGCM (raw)",
  "aifs_calibrated"                 = "AIFS (calibrated)",
  "aifs_calibrated_clim_mok_date"   = "AIFS (calibrated)",
  "ngcm_calibrated"                 = "NGCM (calibrated)",
  "ngcm_blend"                      = "NGCM : EE",
  "int_all"                         = "NGCM : AIFS : EE",
  "ngcm_clim_mok_date_raw"          = "NGCM (raw)",
  "ngcm_calibrated_clim_mok_date"   = "NGCM (calibrated)",
  "mme_clim_mok_date_clim_raw_opt_rps" = "Best MME",
  "mme_no_mok_filter_clim_raw_opt_rps" = "Best MME"
)

METRIC_LABELS <- c(
  "brier_skill" = "BSS",
  "rps_skill"   = "RPSS",
  "auc"         = "AUC"
)

METRIC_COLORS <- c(
  "brier_skill" = "#eb5e71",
  "rps_skill"   = "#93c1dc",
  "auc"         = "#2b7d87"
)

# Build display labels: use known labels where available, else raw name.
# Disambiguate duplicates (e.g. multiple models mapping to "NGCM (calibrated)")
# by appending the raw model name in parentheses.
get_label <- function(m) {
  if (m %in% names(MODEL_LABELS)) MODEL_LABELS[[m]] else m
}
model_labels_vec <- setNames(vapply(all_models, get_label, character(1)), all_models)

dup_labels <- model_labels_vec[duplicated(model_labels_vec) | duplicated(model_labels_vec, fromLast = TRUE)]
for (nm in names(dup_labels)) {
  model_labels_vec[[nm]] <- paste0(dup_labels[[nm]], " (", nm, ")")
}

# Auto-generate color palette for all models
n_all <- length(all_models)
auto_colors <- setNames(
  grDevices::hcl.colors(max(n_all, 3), palette = "Dark 3")[seq_len(n_all)],
  all_models
)

# ---- 5. Plot 1: Model comparison bars (all models) ----
cat("Producing Plot 1: Model comparison bars...\n")

# Use long_format_for_fig4 + plot_2 from blend_figure_utils.R
auc_min_val <- if ("unc_clim_raw" %in% all_models) {
  summ %>% filter(model == "unc_clim_raw") %>% pull(auc)
} else {
  min(summ$auc, na.rm = TRUE)
}

if (length(auc_min_val) == 0) auc_min_val <- 0.5

# Determine reference model for clim_raw lines in plot_2
has_clim_raw <- "clim_raw" %in% all_models

# Order: clim_raw first (if present, used as ref lines), then bar_models
models_for_fig4 <- if (has_clim_raw) {
  c("clim_raw", setdiff(bar_models, "clim_raw"))
} else {
  bar_models
}

long_df <- long_format_for_fig4(
  df            = summ,
  models_keep   = models_for_fig4,
  auc_scale     = 4,
  auc_min       = auc_min_val,
  model_labels  = model_labels_vec,
  metric_labels = METRIC_LABELS
)

p_comparison <- plot_2(
  second_df       = long_df,
  year_label_text = spec_id,
  auc_scale       = 4,
  auc_min         = auc_min_val,
  metric_colors   = METRIC_COLORS,
  metric_labels   = METRIC_LABELS,
  font_size       = 16
) +
  theme(
    legend.position      = "bottom",
    legend.justification = "center"
  ) +
  guides(
    fill     = guide_legend(nrow = 2),
    linetype = guide_legend(nrow = 1)
  )

# Dynamic height: scale with number of models
n_bar <- length(models_for_fig4) - as.integer(has_clim_raw)
fig4_height <- max(3.5, 1.5 + n_bar * 0.4)

save_plot(
  file_stem      = "model_comparison",
  plot           = p_comparison,
  out_dir        = OUTPUT_DIR,
  width          = 7.2,
  height         = fig4_height,
  vector_formats = VECTOR_FORMATS
)

# ---- 6. Plot 2: Weekly skill bars (all models) ----
cat("Producing Plot 2: Weekly performance...\n")

# Detect weekly columns
weekly_brier_cols <- names(summ) %>% str_subset("^brier_week\\d+$")
weekly_auc_cols   <- names(summ) %>% str_subset("^auc_week\\d+$")

if (length(weekly_brier_cols) > 0 || length(weekly_auc_cols) > 0) {
  weekly_cols <- c(weekly_brier_cols, weekly_auc_cols)

  weekly_long <- summ %>%
    filter(model %in% all_models) %>%
    select(model, all_of(weekly_cols)) %>%
    pivot_longer(cols = -model, names_to = "metric_time", values_to = "value") %>%
    separate(metric_time, into = c("metric", "time"), sep = "_", extra = "merge") %>%
    mutate(
      week = as.integer(str_extract(time, "\\d+")),
      week_label = paste0("Week ", week)
    )

  # Compute Brier skill per week: 1 - (brier / baseline_brier)
  baseline_model <- if ("unc_clim_raw" %in% all_models) "unc_clim_raw" else all_models[1]

  brier_baseline <- weekly_long %>%
    filter(model == baseline_model, metric == "brier") %>%
    transmute(week_label, base_brier = value)

  weekly_skill <- weekly_long %>%
    left_join(brier_baseline, by = "week_label") %>%
    mutate(
      plot_value = if_else(metric == "brier", 1 - (value / base_brier), value),
      metric_display = dplyr::recode(metric,
        "brier" = METRIC_LABELS[["brier_skill"]],
        "auc"   = METRIC_LABELS[["auc"]]
      ),
      model = factor(model, levels = all_models),
      week_label = factor(week_label, levels = unique(week_label[order(week)]))
    )

  p_weekly <- ggplot(weekly_skill, aes(x = week_label, y = plot_value, fill = model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey25") +
    facet_grid(rows = vars(metric_display), scales = "free_y", switch = "y") +
    scale_fill_manual(
      values = auto_colors,
      labels = model_labels_vec[all_models],
      name   = "Model"
    ) +
    labs(x = "Lead time", y = NULL, title = NULL) +
    theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      text               = element_text(colour = "#000000"),
      legend.position    = "bottom",
      legend.text        = element_text(size = 7),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#ebebeb", linewidth = 0.25),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(fill = NA, colour = "#1a1a1a", linewidth = 0.5),
      strip.placement    = "outside",
      strip.text.y.left  = element_text(angle = 0, size = 9, colour = "#000000"),
      axis.text          = element_text(size = 7, colour = "#000000"),
      axis.title         = element_text(size = 9, colour = "#000000")
    ) +
    guides(fill = guide_legend(nrow = 3))

  save_plot("weekly_performance", p_weekly, OUTPUT_DIR, 12, 7.5, VECTOR_FORMATS)
} else {
  message("No weekly columns found in summary data - skipping weekly plot.")
}

# ---- 7. Plot 3: Yearly time series (all models) ----
if (!is.null(yearly)) {
  cat("Producing Plot 3: Yearly time series...\n")

  yearly_models <- intersect(all_models, unique(yearly$model))

  if (length(yearly_models) > 0) {
    # Determine which metric columns are available
    yearly_metric_cols <- intersect(c("brier_skill", "rpss", "auc"), names(yearly))

    yearly_long <- yearly %>%
      filter(model %in% yearly_models) %>%
      select(year, model, all_of(yearly_metric_cols)) %>%
      pivot_longer(cols = all_of(yearly_metric_cols), names_to = "metric", values_to = "value") %>%
      mutate(
        metric_display = dplyr::recode(metric,
          "brier_skill" = METRIC_LABELS[["brier_skill"]],
          "rpss"        = METRIC_LABELS[["rps_skill"]],
          "auc"         = METRIC_LABELS[["auc"]]
        ),
        model_display = model_labels_vec[model],
        model = factor(model, levels = yearly_models),
        model_display = factor(model_display, levels = model_labels_vec[yearly_models])
      )

    # Zero reference lines for skill panels
    skill_metrics <- c(METRIC_LABELS[["brier_skill"]], METRIC_LABELS[["rps_skill"]])
    zero_df <- tibble(
      metric_display = skill_metrics[skill_metrics %in% unique(yearly_long$metric_display)],
      yintercept = 0
    )

    n_yearly <- length(yearly_models)
    yearly_colors <- setNames(
      grDevices::hcl.colors(max(n_yearly, 3), palette = "Dark 3")[seq_len(n_yearly)],
      yearly_models
    )

    p_yearly <- ggplot(yearly_long, aes(x = year, y = value, color = model, group = model)) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.3, shape = 16) +
      geom_hline(
        data = zero_df,
        aes(yintercept = yintercept),
        linewidth = 0.25,
        color = "#000000"
      ) +
      facet_grid(rows = vars(metric_display), scales = "free_y", switch = "y") +
      scale_color_manual(
        values = yearly_colors,
        breaks = yearly_models,
        labels = model_labels_vec[yearly_models],
        name   = NULL
      ) +
      scale_x_continuous(breaks = pretty(unique(yearly_long$year))) +
      labs(x = "Year", y = NULL, title = NULL) +
      theme_minimal(base_size = 10, base_family = "Arial") +
      theme(
        text               = element_text(colour = "#000000"),
        legend.position    = "bottom",
        legend.text        = element_text(size = 7),
        panel.grid.major.x = element_line(colour = "#ebebeb", linewidth = 0.25),
        panel.grid.major.y = element_line(colour = "#ebebeb", linewidth = 0.25),
        panel.grid.minor   = element_blank(),
        panel.border       = element_rect(fill = NA, colour = "#1a1a1a", linewidth = 0.5),
        strip.placement    = "outside",
        strip.text.y.left  = element_text(angle = 0, size = 9, colour = "#000000"),
        axis.text          = element_text(size = 7, colour = "#000000"),
        axis.title         = element_text(size = 9, colour = "#000000")
      ) +
      guides(color = guide_legend(nrow = 3))

    save_plot("yearly_time_series", p_yearly, OUTPUT_DIR, 12, 8, VECTOR_FORMATS)
  } else {
    message("No matching models found in yearly data - skipping yearly plot.")
  }
} else {
  message("Yearly metrics not available - skipping yearly plot.")
}

# ---- 8. Plot 4: Metric maps (if data available) ----
if (!is.null(cell_metrics) && !is.null(map_inputs)) {
  cat("Producing Plot 4: Metric maps...\n")

  polygons_df    <- map_inputs$polygons_df
  india_boundary <- map_inputs$india_boundary

  # Determine which two models to compare
  clim_model  <- if ("clim_raw" %in% all_models) "clim_raw" else NULL
  final_model <- if ("blended_model" %in% all_models) "blended_model" else NULL

  # If either is missing, try first two available models
  available_cell_models <- unique(cell_metrics$model)
  available_cell_models <- available_cell_models[available_cell_models != "unc_clim_raw"]

  if (is.null(clim_model) || !clim_model %in% available_cell_models) {
    clim_model <- if (length(available_cell_models) >= 1) available_cell_models[1] else NULL
  }
  if (is.null(final_model) || !final_model %in% available_cell_models) {
    final_model <- if (length(available_cell_models) >= 2) available_cell_models[2] else NULL
  }

  if (!is.null(clim_model) && !is.null(final_model) && clim_model != final_model) {

    # Clean cell_metrics: remove ALL row, ensure numeric lat/lon
    cell_clean <- cell_metrics %>%
      filter(id != "ALL") %>%
      mutate(lat = as.numeric(lat), lon = as.numeric(lon))

    cell_comp <- summarize_maps_compare(
      all_cells   = cell_clean,
      method      = "global",
      clim_model  = clim_model,
      final_model = final_model
    )

    clim_label  <- get_label(clim_model)
    final_label <- get_label(final_model)

    map_specs <- tibble::tribble(
      ~metric_col,   ~title_text,                                                          ~file_stem,
      "brier_skill", paste0("Brier Skill (", final_label, " vs ", clim_label, ")"),         "map_brier_skill",
      "rps_skill",   paste0("RPS Skill (", final_label, " vs ", clim_label, ")"),            "map_rps_skill",
      "auc_diff",    paste0("AUC Difference (", final_label, " - ", clim_label, ")"),        "map_auc_diff"
    )

    maps_out_dir <- file.path(OUTPUT_DIR, "maps")
    dir.create(maps_out_dir, showWarnings = FALSE, recursive = TRUE)

    purrr::pwalk(map_specs, function(metric_col, title_text, file_stem) {
      leg_title <- legend_title_for_metric(metric_col)

      p <- plot_metric_map(
        cell_metrics = cell_comp,
        polygons_df  = polygons_df,
        metric_col   = metric_col,
        title_text   = title_text,
        india_map    = india_boundary,
        legend_title = leg_title
      ) +
        theme_minimal(base_size = 16) +
        theme(
          plot.title      = element_text(face = "bold"),
          legend.position = "right"
        )

      save_plot(
        file_stem      = file_stem,
        plot           = p,
        out_dir        = maps_out_dir,
        width          = 10,
        height         = 10,
        vector_formats = VECTOR_FORMATS
      )
    })
  } else {
    message("Need at least 2 distinct models in cell_metrics for maps - skipping.")
  }
} else {
  message("Cell metrics or map inputs not found - skipping maps.")
}

cat("\nDone! Figures saved to: ", OUTPUT_DIR, "\n")
