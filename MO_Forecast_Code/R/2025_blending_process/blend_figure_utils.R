# ==============================================================================
# File: blend_figure_utils.R
# ==============================================================================
# Purpose
#   Helper functions for 3_produce_figures.R. Provides safe RDS reading,
#   column checks, plot saving, model label inference, secondary axis helpers,
#   and specialized plot constructors for skill comparisons, weekly bins,
#   yearly time series, reliability diagrams, and metric maps.
#
## Function index
#   safe_read_rds(path)
#     Read RDS with informative error on failure.
#
#   optional_read_rds(path)
#     Read RDS if file exists, else return NULL.
#
#   check_required_cols(df, required_cols, name_for_error)
#     Stop with error if df is missing any required columns.
#
#   add_period_tag(df, period_label)
#     Add a 'period' column with a fixed label to a data frame.
#
#   fmt_num(x, digits)
#     Format a number to a fixed number of decimal places.
#
#   save_plot(file_stem, plot, out_dir, width, height, vector_formats, ...)
#     Save a ggplot to PDF/SVG/PNG with sensible defaults.
#
#   infer_train_window(model)
#     Parse training window label from model name string (e.g. "2000-2024").
#
#   sec_axis_inv(x, auc_scale, auc_min)
#     Inverse transformation for AUC secondary axis.
#
#   long_format_for_fig4(df, models_keep, auc_scale, auc_min, model_labels, ...)
#     Pivot model metrics to long format for grouped bar charts.
#
#   plot_fig4(second_df, year_label_text, auc_scale, auc_min, ...)
#     Horizontal bar chart of Brier skill, RPS skill, and AUC
#     (matching format of figure 4 in main text).
#      Reference lines for Evolving Expectations Model.
#
#   make_fig4_variant(summ_df, models_keep, label_text, file_stem, out_dir)
#     wrapper for long_format_for_fig4 and plot_fig4, creates a variant using either
#     different year or different onset condition 
#
#
#   make_weekly_bins_plot(df, title_suffix, variant, include_later, ...)
#     Faceted bar chart of per-week Brier skill and AUC by model.
#
#   make_yearly_plot(yearly_df, title_suffix, model_labels, ...)
#     Faceted line plot of yearly Brier skill, RPSS, and AUC by model.
#
#   pretty_model_name(m, model_labels)
#     Look up display name for a model string using model_labels.
#
#   make_chartdata_path(reliability_dir, model, period_tag, bins_tag)
#     Build file path for a reliability chart-data RDS file.
#
#   read_chartdata_one(path, model, model_pretty, period_pretty)
#     Read one reliability chart-data RDS and attach model/period metadata.
#
#   read_india_boundary(path)
#     Read India boundary from RDS (sf or data.frame) or CSV.
#     Returns data.frame with long/lat/group columns for geom_polygon.
#
#   legend_title_for_metric(metric_col)
#     Return "Skill" or "Difference" legend title based on metric column name.
#
#   plot_metric_map(cell_metrics, polygons_df, metric_col, title_text, ...)
#     Choropleth map of a per-cell metric with optional dissemination-cell overlay.
#
#   plot_reliability_3panel(rel_data, model_order, model_labels, ...)
#      3-panel reliability diagram (in the style of figure 3).
#
#   plot_fig2_combined(summ_df, yearly_df, variant, model_labels, ...)
#     Pixel-matched 4-panel combined figure: AUC bars, Brier bars,
#     AUC by year, Brier skill by year (in the style of figure 2) 
#
# Dependencies
#   - dplyr, tidyr, stringr, purrr, ggplot2, readr, sf (optional for maps)
# ==============================================================================

safe_read_rds <- function(path) {
  tryCatch(
    readRDS(path),
    error = function(e) stop(sprintf("Failed to read RDS at '%s': %s", path, e$message), call. = FALSE)
  )
}

# Read RDS if it exists, otherwise return NULL (for optional inputs)
optional_read_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

check_required_cols <- function(df, required_cols, name_for_error) {
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0) {
    stop(sprintf("Data frame '%s' is missing columns: %s", name_for_error, paste(miss, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

add_period_tag <- function(df, period_label) {
  df %>% mutate(period = period_label)
}

fmt_num <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

save_plot <- function(
    file_stem,
    plot,
    out_dir,
    width = 12,
    height = 8,
    vector_formats = c("pdf", "svg"),
    raster_format = NULL,
    raster_dpi = 300,
    bg = "white",
    print_plot = FALSE
) {
  stopifnot(inherits(plot, "ggplot"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(print_plot)) print(plot)

  # Vector outputs — use cairo_pdf/svglite to properly embed Arial
  for (fmt in vector_formats) {
    out_path <- file.path(out_dir, paste0(file_stem, ".", fmt))
    dev <- if (fmt == "pdf" && capabilities("cairo")) cairo_pdf
           else if (fmt == "pdf") grDevices::pdf
           else if (fmt == "svg" && requireNamespace("svglite", quietly = TRUE)) svglite::svglite
           else fmt
    ggsave(filename = out_path, plot = plot, width = width, height = height, bg = bg, device = dev)
    message("Saved: ", out_path)
  }

  # Optional raster
  if (!is.null(raster_format) && nzchar(raster_format)) {
    out_path <- file.path(out_dir, paste0(file_stem, ".", raster_format))
    ggsave(filename = out_path, plot = plot, width = width, height = height, dpi = raster_dpi, bg = bg, device = raster_format)
    message("Saved: ", out_path)
  }

  invisible(TRUE)
}

# Derive training window label from model name
# Examples:
# - clim_raw -> "1900-2024"
# - clim_raw_2000_2024 -> "2000-2024"
# - blend_clim_mok_date_clim_raw_opt_rps -> "1900-2024" (treated as baseline)
infer_train_window <- function(model) {
  m <- as.character(model)

  if (m %in% c("clim_raw", "blended_model", "blend_clim_mok_date_clim_raw_opt_rps")) return("1900-2024")

  mm <- stringr::str_match(m, ".*_(\\d{4})_(\\d{4})(?:_opt_rps)?$")
  if (is.na(mm[1, 1])) return(NA_character_)
  paste0(mm[1, 2], "-", mm[1, 3])
}

sec_axis_inv <- function(x, auc_scale, auc_min) {
  (x / auc_scale) + auc_min
}

long_format_for_fig4 <- function(df, models_keep, auc_scale, auc_min, model_labels, metric_labels) {
  df %>%
    filter(model %in% models_keep) %>%
    select(model, brier_skill, rps_skill, auc) %>%
    pivot_longer(cols = c(brier_skill, rps_skill, auc), names_to = "metric", values_to = "value") %>%
    mutate(
      metric_display = recode(metric, !!!metric_labels),
      value_centered = if_else(metric == "auc", (value - auc_min) * auc_scale, value),
      model = factor(model, levels = models_keep),
      model_display = recode(as.character(model), !!!model_labels, .default = as.character(model)),
      model_display = factor(
        model_display,
        levels = recode(models_keep, !!!model_labels, .default = models_keep)
      )
    )
}


# ==============================================================================
# Skill-comparison horizontal bar chart (fig_5 style — pixel-matched)
# ==============================================================================
# Fig 5 skill comparison color palette
#   Bar fills:        Brier=#eb5e71  RPS=#93c1dc  AUC=#2b7d87
#   Ref-line strokes: Brier=#a72429  RPS=#93c1dc  AUC=#2e67ab  (SOLID, thin)
#   ALL text:         #000000, ArialMT, NOT bold
#   Y-axis labels:    10pt, ArialMT, plain
#   Bottom ticks:     6pt;  Top ticks: 7pt
#   Axis titles:      10pt, black
#   Grid:             #d9d9d9, width 0.25
#   Zero line:        #333333, SOLID, width 1.0
#   Bottom breaks:    every 0.04
#   Top breaks:       every 0.01
#   Legend:           inside plot, bottom-right
#   Legend order:     Evolving Exp. (vline glyph), AUC, RPS, Brier (square glyphs)
# ------------------------------------------------------------------------------


plot_2 <- function(second_df,
                        year_label_text,
                        auc_scale,
                        auc_min,
                        metric_colors,
                        metric_labels,
                        font_size = 16) {
  
  bars_df <- second_df %>% dplyr::filter(model != "clim_raw")
  
  # Ensure metric factor order matches reference: AUC (top), RPS (mid), Brier (bottom)
  metric_order <- c("auc", "rps_skill", "brier_skill")
  bars_df$metric <- factor(bars_df$metric, levels = metric_order)
  
  vline_df <- second_df %>%
    dplyr::filter(model == "clim_raw") %>%
    dplyr::distinct(metric, value_centered)
  
  # ---- Reference-line colours (DIFFERENT from bar fills) ----
  vline_colors <- c(
    "brier_skill" = "#a72429",
    "rps_skill"   = "#93c1dc",
    "auc"         = "#2e67ab"
  )
  
  # ---- Axis breaks: every 0.04 (bottom), every 0.01 (top) ----
  all_vals <- c(bars_df$value_centered, vline_df$value_centered)
  x_lo <- floor(min(all_vals, na.rm = TRUE) / 0.04) * 0.04
  x_hi <- ceiling(max(all_vals, na.rm = TRUE) / 0.04) * 0.04
  bottom_breaks <- seq(x_lo, x_hi, by = 0.04)
  top_breaks    <- seq(x_lo / auc_scale, x_hi / auc_scale, by = 0.01)
  
  # "0" for zero; drop trailing zeros otherwise
  clean_label <- function(x) {
    ifelse(x == 0, "0", sub("0+$", "", sub("\\.$", "", sprintf("%.2f", x))))
  }
  
  # Dummy for "Evolving Exp." legend entry — solid vline key glyph
  dummy_line <- data.frame(lab = "Evolving Exp.")
  
  ggplot2::ggplot(bars_df,
                  ggplot2::aes(x = value_centered, y = model_display, fill = metric)) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(width = 0.8, reverse = TRUE),
      width = 0.75
    ) +
    # Coloured reference lines — SOLID, separate colours from bar fills
    ggplot2::geom_vline(
      data = vline_df,
      ggplot2::aes(xintercept = value_centered, color = metric),
      linetype = "solid",
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    # Invisible dummy segment for "Evolving Exp." legend entry
    ggplot2::geom_segment(
      data = dummy_line,
      ggplot2::aes(x = 0, xend = 0, y = 1, yend = 1, linetype = lab),
      inherit.aes = FALSE,
      linewidth = 0.5,
      alpha = 0,
      key_glyph = "vline"
    ) +
    # Zero line — SOLID
    ggplot2::geom_vline(xintercept = 0, linewidth = 1.0, color = "#333333") +
    # Fill scale (bars) — order: AUC, RPS, Brier top-to-bottom in legend
    ggplot2::scale_fill_manual(
      values = metric_colors,
      name   = NULL,
      breaks = c("auc", "rps_skill", "brier_skill"),
      labels = c("AUC", "RPSS", "BSS")
    ) +
    # Colour scale (ref lines — NOT in legend)
    ggplot2::scale_color_manual(values = vline_colors, guide = "none") +
    # Linetype — "Evolving Exp." with solid vline glyph
    ggplot2::scale_linetype_manual(
      values = c("Evolving Exp." = "solid"),
      name   = NULL,
      guide  = ggplot2::guide_legend(
        order = 1,
        override.aes = list(linewidth = 0.8, alpha = 1, colour = "#2e67ab")
      )
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::scale_x_continuous(
      breaks       = bottom_breaks,
      labels       = clean_label,
      minor_breaks = NULL,
      name         = "BSS/RPSS",
      sec.axis     = ggplot2::sec_axis(
        ~ . / auc_scale,
        name   = "AUC Difference",
        breaks = top_breaks,
        labels = clean_label
      )
    ) +
    ggplot2::coord_cartesian(xlim = c(x_lo, x_hi)) +
    ggplot2::scale_y_discrete(name = NULL) +
    ggplot2::theme_minimal(base_size = 10, base_family = "Arial") +
    ggplot2::theme(
      legend.position      = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background    = ggplot2::element_rect(fill = "white", colour = NA),
      legend.box           = "vertical",
      legend.text          = ggplot2::element_text(size = 10, colour = "#000000"),
      legend.key.size      = ggplot2::unit(0.45, "cm"),
      panel.grid.major.y   = ggplot2::element_blank(),
      panel.grid.major.x   = ggplot2::element_line(color = "#d9d9d9", linewidth = 0.25),
      panel.grid.minor     = ggplot2::element_blank(),
      axis.ticks.x         = ggplot2::element_line(color = "#000000", linewidth = 0.25),
      axis.text.y          = ggplot2::element_text(size = 10, colour = "#000000"),
      axis.text.x.bottom   = ggplot2::element_text(size = 6, colour = "#000000"),
      axis.text.x.top      = ggplot2::element_text(size = 7, colour = "#000000"),
      axis.title.x.bottom  = ggplot2::element_text(size = 10, colour = "#000000"),
      axis.title.x.top     = ggplot2::element_text(size = 10, colour = "#000000"),
      plot.title           = ggplot2::element_blank(),
      plot.margin          = ggplot2::margin(5, 5, 5, 5),
      legend.spacing.y = ggplot2::unit(-.225, "cm"),
    )
}


# Creates and saves a variant of figure 4, using plot_fig4 as a helper
make_fig4_variant <- function(summ_df, models_keep, label_text, file_stem, out_dir) {
  if (is.null(summ_df)) return(NULL)
  auc_min_val <- summ_df %>% filter(model == "unc_clim_raw") %>% pull(auc)
  if (length(auc_min_val) == 0) return(NULL)
  
  long_df <- long_format_for_fig4(
    summ_df,
    models_keep = models_keep,
    auc_scale   = 4,
    auc_min     = auc_min_val,
    model_labels = MODEL_LABELS,
    metric_labels = METRIC_LABELS
  )
  
  p <- plot_2(
    second_df       = long_df,
    year_label_text = label_text,
    auc_scale       = 4,
    auc_min         = auc_min_val,
    metric_colors   = METRIC_COLORS,
    metric_labels   = METRIC_LABELS,
    font_size       = FONT_SIZE
  )
  
  save_plot(file_stem, p, out_dir, 7.2, 3.69, VECTOR_FORMATS_DEFAULT, RASTER_FORMAT_DEFAULT, RASTER_DPI_DEFAULT)
  p
}



make_weekly_bins_plot <- function(
    df,
    title_suffix,
    variant = c("standard", "clim_mok_date", "no_mok_filter"),
    include_later = FALSE,
    model_labels,
    metric_labels,
    model_colors,
    font_size
) {
  variant <- match.arg(variant)

  # Bar models by variant
  if (variant %in% c("standard", "clim_mok_date")) {
    bar_levels <- c("unc_clim_raw", "clim_raw", "ngcm_calibrated_clim_mok_date", "blended_model")
  } else {
    bar_levels <- c("unc_clim_raw", "clim_raw", "ngcm_calibrated", "blended_model")
  }

  bar_levels <- bar_levels[bar_levels %in% unique(df$model)]

  weekly_cols_exist <- names(df) %>%
    stringr::str_subset(if (include_later) "^(brier|auc)_(week\\d+|later)$" else "^(brier|auc)_week\\d+$")

  if (length(weekly_cols_exist) == 0) {
    stop(sprintf("No weekly columns found in this dataset for %s.", title_suffix), call. = FALSE)
  }

  weekly_long <- df %>%
    filter(model %in% bar_levels) %>%
    select(model, all_of(weekly_cols_exist)) %>%
    pivot_longer(cols = -model, names_to = "metric_time", values_to = "value") %>%
    separate(metric_time, into = c("metric", "time"), sep = "_", remove = TRUE) %>%
    mutate(
      week = stringr::str_extract(time, "\\d+") %>% as.integer(),
      week = dplyr::coalesce(week, dplyr::if_else(time == "later", 99L, NA_integer_)),
      week_label = dplyr::if_else(time == "later", "Later", paste0("Week", week))
    )

  brier_baseline <- weekly_long %>%
    filter(model == "unc_clim_raw", metric == "brier") %>%
    transmute(week_label, base_brier = value)

  weekly_skill <- weekly_long %>%
    left_join(brier_baseline, by = "week_label") %>%
    mutate(
      plot_value = if_else(metric == "brier", 1 - (value / base_brier), value),
      metric_display = recode(metric, "brier" = metric_labels[["brier_skill"]], "auc" = metric_labels[["auc"]])
    )

  bars_df <- weekly_skill %>%
    mutate(
      model = factor(model, levels = bar_levels),
      week_label = factor(week_label, levels = unique(week_label))
    )

  fill_vals <- model_colors[bar_levels]
  names(fill_vals) <- bar_levels

  auc_label <- metric_labels[["auc"]]

  ggplot(bars_df, aes(x = week_label, y = plot_value, fill = model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey25") +
    facet_grid(rows = vars(metric_display), scales = "free_y", switch = "y") +
    scale_fill_manual(
      values = fill_vals,
      breaks = bar_levels,
      limits = bar_levels,
      drop = FALSE,
      labels = unname(model_labels[bar_levels]),
      name = "Model"
    ) +
    labs(x = "Lead time", y = NULL, title = NULL) +
    theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      text               = element_text(colour = "#000000"),
      legend.position    = "top",
      legend.text        = element_text(size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#ebebeb", linewidth = 0.25),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(fill = NA, colour = "#1a1a1a", linewidth = 0.5),
      strip.placement    = "outside",
      strip.text.y.left  = element_text(angle = 0, size = 9, colour = "#000000"),
      axis.text          = element_text(size = 7, colour = "#000000"),
      axis.title         = element_text(size = 9, colour = "#000000")
    )
}

make_yearly_plot <- function(yearly_df, title_suffix, model_labels, metric_labels, okabe_ito, font_size,
                             ngcm_model = "ngcm_calibrated_clim_mok_date") {
  # Lines/points series (main)
  models_main <- c("unc_clim_raw", "clim_raw", ngcm_model, "blended_model")
  models_needed <- models_main

  model_colors_yearly <- c(
    "unc_clim_raw"                = "#237192",
    "clim_raw"                    = "#70c8d4",
    "ngcm_calibrated_clim_mok_date" = "#eb7900",
    "blended_model"               = "#ca1b00",
    "ngcm_calibrated"               = "#eb7900"
  )
  
  
  yearly_long <- yearly_df %>%
    filter(model %in% models_needed) %>%
    select(year, model, brier_skill, rpss, auc) %>%
    pivot_longer(cols = c(brier_skill, rpss, auc), names_to = "metric", values_to = "value") %>%
    mutate(
      metric_display = recode(
        metric,
        "brier_skill" = metric_labels[["brier_skill"]],
        "rpss"        = metric_labels[["rps_skill"]],
        "auc"         = metric_labels[["auc"]]
      ),
      model_display = recode(model, !!!model_labels, .default = model)
    )

  yearly_lines <- yearly_long %>%
    filter(model %in% models_main) %>%
    mutate(
      model = factor(model, levels = models_main),
      model_display = factor(model_display, levels = recode(models_main, !!!model_labels, .default = models_main))
    )

  zero_df <- tibble::tibble(
    metric_display = c(metric_labels[["brier_skill"]], metric_labels[["rps_skill"]]),
    yintercept = 0
  )

  
  ggplot(yearly_lines, aes(x = year, y = value, color = model, group = model)) +
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
      values = model_colors_yearly[models_main],
      breaks = models_main,
      labels = unname(model_labels[models_main]),
      name = NULL
    ) +
    scale_x_continuous(breaks = pretty(unique(yearly_lines$year))) +
    labs(x = "Year", y = NULL, title = NULL) +
    theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      text               = element_text(colour = "#000000"),
      legend.position    = "top",
      legend.text        = element_text(size = 8),
      panel.grid.major.x = element_line(colour = "#ebebeb", linewidth = 0.25),
      panel.grid.major.y = element_line(colour = "#ebebeb", linewidth = 0.25),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(fill = NA, colour = "#1a1a1a", linewidth = 0.5),
      strip.placement    = "outside",
      strip.text.y.left  = element_text(angle = 0, size = 9, colour = "#000000"),
      axis.text          = element_text(size = 7, colour = "#000000"),
      axis.title         = element_text(size = 9, colour = "#000000")
    )
}

pretty_model_name <- function(m, model_labels) {
  dplyr::recode(m, !!!model_labels, .default = m)
}

make_chartdata_path <- function(reliability_dir, model, period_tag, bins_tag = "10bins") {
  fn <- paste0("reliability_quantiles_global_", model, "_", bins_tag, "_", period_tag, "_chartdata.rds")
  file.path(reliability_dir, fn)
}

read_chartdata_one <- function(path, model, model_pretty, period_pretty) {
  df <- readRDS(path) %>%
    mutate(
      model = model,
      model_pretty = model_pretty,
      period_pretty = period_pretty
    )

  bw <- df %>%
    filter(series == "hist") %>%
    mutate(bin_width = bin_right - bin_left) %>%
    summarise(bin_width = dplyr::first(bin_width)) %>%
    pull(bin_width)

  df %>% mutate(bin_width = bw)
}
# Read India boundary from .rds or .csv:
# - If it's an sf object, convert to a fortify-like data.frame for geom_polygon.
# - If it's already a data.frame with long/lat/group, just return it.
read_india_boundary <- function(path) {
  ext <- tolower(tools::file_ext(path))
  obj <- if (ext == "csv") readr::read_csv(path, show_col_types = FALSE) else safe_read_rds(path)

  if (inherits(obj, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Package 'sf' is required to read sf-format india_boundary. Install sf or supply a data.frame boundary.")
    # Cast to polygons then build a dataframe with long/lat/group
    obj <- sf::st_cast(obj, "MULTIPOLYGON", warn = FALSE)

    # Convert to coordinates dataframe
    coords <- sf::st_coordinates(obj)
    # coords columns include X, Y, L1, L2, L3 ... depending on geometry nesting
    # We'll define a stable group id using all available L* columns.
    l_cols <- grep("^L\\d+$", colnames(coords), value = TRUE)
    group_id <- apply(coords[, l_cols, drop = FALSE], 1, paste, collapse = "_")

    df <- tibble::tibble(
      long = coords[, "X"],
      lat  = coords[, "Y"],
      group = group_id
    )
    return(df)
  }

  if (is.data.frame(obj) && all(c("long", "lat", "group") %in% names(obj))) {
    return(obj)
  }

  stop("india_boundary_path must be an sf object or a data.frame with columns long, lat, group.", call. = FALSE)
}



#' Return a legend title string based on the metric column name.
legend_title_for_metric <- function(metric_col) {
  if (grepl("skill|rpss", metric_col, ignore.case = TRUE)) "Skill" else "Difference"
}

#' Plot a metric map for cell polygons with optional dissemination-cell highlighting.
plot_metric_map <- function(cell_metrics, polygons_df, metric_col, title_text, india_map,
                            legend_title = metric_col, dissem_poly_df = NULL) {
  poly_data <- polygons_df %>%
    dplyr::left_join(cell_metrics, by = "id") %>%
    dplyr::filter(is.na(.data[[metric_col]]) | .data[[metric_col]] >= -.2)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(x = lon, y = lat, group = id, fill = .data[[metric_col]]),
      color = "black", linewidth = 0.2
    )
  
  if (!is.null(dissem_poly_df) && nrow(dissem_poly_df) > 0) {
    p <- p +
      ggplot2::geom_polygon(
        data = dissem_poly_df,
        ggplot2::aes(x = lon, y = lat, group = id),
        fill = NA, color = "black", linewidth = 1.6
      )
  }
  
  p <- p +
    ggplot2::geom_polygon(
      data = india_map,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = NA, color = "blue"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title_text, fill = legend_title) +
    ggplot2::theme_minimal()
  
  if (grepl("diff|skill", metric_col)) {
    p <- p + ggplot2::scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0)
  }
  
  p
}



# ==============================================================================
# Reliability 3-panel plot (style of figure 3 in text)
# ==============================================================================
# Reliability diagram color palette
#   Histogram fills: A=#c49a2c  B=#eb7900  C=#ca1b00  (alpha = 1.0, FULL opacity)
#   Grid:           #ebebeb, width 0.25
#   Axis tick text: #231f20, ArialMT, 8pt
#   Axis titles:    #000000, ArialMT, 10pt
#   Model names:    #1a1a1a, ArialMT, 10pt, plain
#   Cal. line:      #000000, width ~1.9, solid
#   Cal. dots:      #000000, size ~3pt
#   Diagonal ref:   #000000, dashed [6 4], width 0.75
#   Panel border:   #000000, width 1.0
#   Badge letter:   MyriadPro-Bold 12pt, white on coloured square
#   Font:           Arial (ArialMT in the PDF)
# ------------------------------------------------------------------------------
plot_reliability_3panel <- function(rel_data,
                                    model_order,
                                    model_labels,
                                    reliability_colors,
                                    font_size = 16) {
  
  panel_ids <- setNames(seq_along(model_order), model_order)
  
  cal_df  <- rel_data %>% dplyr::filter(series == "calibration")
  hist_df <- rel_data %>% dplyr::filter(series == "hist")
  
  cal_df <- cal_df %>%
    dplyr::mutate(panel_id = factor(panel_ids[model], levels = seq_along(model_order)))
  hist_df <- hist_df %>%
    dplyr::mutate(panel_id = factor(panel_ids[model], levels = seq_along(model_order)))
  
  bw <- hist_df$bin_width[1]
  
  badge_df <- data.frame(
    panel_id   = factor(seq_along(model_order), levels = seq_along(model_order)),
    model      = model_order,
    letter     = LETTERS[seq_along(model_order)],
    model_name = unname(model_labels[model_order]),
    stringsAsFactors = FALSE
  )
  
  # "0" not "0.00"
  axis_fmt <- function(x) ifelse(x == 0, "0", sprintf("%.2f", x))
  
  ggplot2::ggplot() +
    # Histogram bars — FULL OPACITY (alpha = 1.0)
    ggplot2::geom_col(
      data = hist_df,
      ggplot2::aes(x = bin_mid, y = y, fill = model),
      width = bw,
      show.legend = FALSE
    ) +
    # Diagonal reference — dashed, width 0.75
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         linewidth = 0.75, colour = "#000000") +
    # Calibration line — thick solid black, width ~1.9
    ggplot2::geom_line(
      data = cal_df,
      ggplot2::aes(x = pred_mean, y = obs_frac),
      linewidth = 1.9, colour = "#000000"
    ) +
    # Calibration dots — larger (~3pt)
    ggplot2::geom_point(
      data = cal_df,
      ggplot2::aes(x = pred_mean, y = obs_frac),
      size = 3, colour = "#000000"
    ) +
    # Badge: coloured square with white bold letter
    ggplot2::geom_tile(
      data = badge_df,
      ggplot2::aes(x = 0.045, y = 0.955, fill = model),
      width = 0.065, height = 0.065,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = badge_df,
      ggplot2::aes(x = 0.045, y = 0.955, label = letter),
      color = "white", fontface = "bold",
      family = "Arial",
      size = 4.2   # ~12pt
    ) +
    # Badge: model name (plain face, #1a1a1a)
    ggplot2::geom_text(
      data = badge_df,
      ggplot2::aes(x = 0.13, y = 0.955, label = model_name),
      hjust = 0, vjust = 0.5,
      fontface = "plain", color = "#1a1a1a",
      family = "Arial",
      size = 3.5   # ~10pt
    ) +
    ggplot2::scale_fill_manual(values = reliability_colors) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = axis_fmt
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = axis_fmt
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::facet_wrap(~ panel_id, nrow = 1) +
    ggplot2::labs(x = "Forecast probability", y = "Observed frequency") +
    ggplot2::theme_minimal(base_size = 10, base_family = "Arial") +
    ggplot2::theme(
      panel.grid.major   = ggplot2::element_line(colour = "#ebebeb", linewidth = 0.25),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(fill = NA, colour = "#000000",
                                                 linewidth = 1.0),
      strip.text         = ggplot2::element_blank(),
      plot.title         = ggplot2::element_blank(),
      axis.text          = ggplot2::element_text(colour = "#231f20", size = 8),
      axis.title.x       = ggplot2::element_text(size = 10, colour = "#000000"),
      axis.title.y       = ggplot2::element_text(size = 10, colour = "#000000"),
      plot.margin        = ggplot2::margin(5, 5, 5, 5)
    )
}


# ==============================================================================
# Combined 4-panel figure (style of figure 2 in main text)
# ==============================================================================
# Fig 2 combined plot color palette
#   Model colors:     Static Clim=#237192  Evolving Exp=#70c8d4
#                     NGCM(cal)=#eb7900    Blended=#ca1b00
#   ALL text:         ArialMT, #000000
#   Legend text:       7pt
#   Panel labels:     Arial-BoldMT 10pt (A, B, C, D)
#   Axis ticks:       6pt
#   Axis titles:      8pt
#   Lines:            width 1.0, solid, filled circle markers ~2.6
#   Grid:             #ebebeb, width 0.25
#   Bars:             no outline stroke
#   Zero line (D):    #000000, width 0.25
#   Page:             342.99 x 304.37 pt -> aspect ~ 4.76 x 4.23 in
#   Layout:           top = A|B bars, middle = C lines, bottom = D lines
# ------------------------------------------------------------------------------

# Colours for the 4-model palette used in fig_3
FIG2_MODEL_COLORS <- c(
  "unc_clim_raw"                = "#237192",
  "clim_raw"                    = "#70c8d4",
  "ngcm_calibrated_clim_mok_date" = "#eb7900",
  "ngcm_calibrated"               = "#eb7900",
  "ngcm_calibrated_mok"           = "#eb7900",
  "blended_model"               = "#ca1b00"
)

# Shared theme for all 4 panels
.fig2_theme <- function() {
  ggplot2::theme_minimal(base_size = 8, base_family = "Arial") +
    ggplot2::theme(
      text                = ggplot2::element_text(colour = "#000000"),
      panel.grid.major    = ggplot2::element_line(colour = "#ebebeb", linewidth = 0.25),
      panel.grid.minor    = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_blank(),
      axis.text           = ggplot2::element_text(size = 6, colour = "#000000"),
      axis.title          = ggplot2::element_text(size = 8, colour = "#000000"),
      plot.title          = ggplot2::element_text(size = 8, colour = "#000000", face = "plain",
                                                  hjust = 0, margin = ggplot2::margin(0, 0, 2, 0)),
      plot.title.position = "plot",
      legend.position     = "none",
      plot.margin         = ggplot2::margin(2, 4, 2, 4)
    )
}

plot_fig2_combined <- function(summ_df,
                               yearly_df,
                               variant = c("standard", "clim_mok_date", "no_mok_filter"),
                               model_labels,
                               model_colors = FIG2_MODEL_COLORS) {
  
  variant <- match.arg(variant)
  
  # ---- Determine models ----
  if (variant %in% c("standard", "clim_mok_date")) {
    bar_models <- c("unc_clim_raw", "clim_raw", "ngcm_calibrated_clim_mok_date", "blended_model")
  } else {
    bar_models <- c("unc_clim_raw", "clim_raw", "ngcm_calibrated", "blended_model")
  }
  bar_models <- bar_models[bar_models %in% unique(summ_df$model)]
  line_models <- bar_models[bar_models %in% unique(yearly_df$model)]
  
  fill_vals <- model_colors[bar_models]
  col_vals  <- model_colors[line_models]
  model_labs <- unname(model_labels[bar_models])
  
  # ---- Shared theme for bar panels (A, B) ----
  bar_theme <- function() {
    .fig2_theme() +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(fill = NA, colour = "#1a1a1a", linewidth = 0.5),
        plot.margin  = ggplot2::margin(2, 6, 2, 4)
      )
  }
  
  # ---- Shared theme for line panels (C, D) ----
  line_theme <- function() {
    .fig2_theme() +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_line(colour = "#ebebeb", linewidth = 0.25),
        panel.border       = ggplot2::element_rect(fill = NA, colour = "#1a1a1a", linewidth = 0.5),
        plot.margin        = ggplot2::margin(2, 6, 2, 4)
      )
  }
  
  # ---- A. Extract weekly data ----
  week_cols <- names(summ_df) %>% stringr::str_subset("^(brier|auc)_week\\d+$")
  weekly_long <- summ_df %>%
    dplyr::filter(model %in% bar_models) %>%
    dplyr::select(model, dplyr::all_of(week_cols)) %>%
    tidyr::pivot_longer(-model, names_to = "mc", values_to = "value") %>%
    tidyr::separate(mc, into = c("metric", "time"), sep = "_", extra = "merge") %>%
    dplyr::mutate(
      week = as.integer(stringr::str_extract(time, "\\d+")),
      model = factor(model, levels = bar_models)
    )
  
  # Brier skill = 1 - brier/baseline_brier
  brier_base <- weekly_long %>%
    dplyr::filter(model == "unc_clim_raw", metric == "brier") %>%
    dplyr::transmute(week, base_brier = value)
  
  weekly_plot <- weekly_long %>%
    dplyr::left_join(brier_base, by = "week") %>%
    dplyr::mutate(
      plot_value = dplyr::if_else(metric == "brier", 1 - (value / base_brier), value),
      lead_label = factor(paste0(week, " weeks"), levels = paste0(4:1, " weeks"))
    )
  levels(weekly_plot$lead_label)[levels(weekly_plot$lead_label) == "1 weeks"] <- "1 week"
  
  # ---- Panel A: AUC bars ----
  auc_bars <- weekly_plot %>% dplyr::filter(metric == "auc")
  pA <- ggplot2::ggplot(auc_bars, ggplot2::aes(x = lead_label, y = plot_value, fill = model)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::scale_fill_manual(values = fill_vals, labels = model_labs, name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::coord_cartesian(ylim = c(0.5, NA)) +
    ggplot2::labs(x = "Lead time", y = NULL) +
    ggplot2::ggtitle(expression(bold("A") ~ " Area under ROC curve")) +
    bar_theme() +
    ggplot2::theme(legend.position = "none")
  
  # ---- Panel B: Brier skill bars + legend ----
  brier_bars <- weekly_plot %>% dplyr::filter(metric == "brier")
  pB <- ggplot2::ggplot(brier_bars, ggplot2::aes(x = lead_label, y = plot_value, fill = model)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::scale_fill_manual(values = fill_vals, labels = model_labs, name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = "Lead time", y = NULL) +
    ggplot2::ggtitle(expression(bold("B") ~ " Brier skill score")) +
    bar_theme() +
    ggplot2::theme(
      legend.position      = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background    = ggplot2::element_rect(fill = "white", colour = NA),
      legend.text          = ggplot2::element_text(size = 6, colour = "#000000"),
      legend.key.size      = ggplot2::unit(0.3, "cm"),
      legend.spacing.y     = ggplot2::unit(0.02, "cm"),
      legend.margin        = ggplot2::margin(1, 3, 1, 3)
    )
  
  # ---- C/D. yearly time series ----
  yearly_long <- yearly_df %>%
    dplyr::filter(model %in% line_models) %>%
    dplyr::select(year, model, brier_skill, auc) %>%
    dplyr::mutate(model = factor(model, levels = line_models))
  
  year_range <- range(yearly_long$year, na.rm = TRUE)
  yr_breaks  <- seq(floor(year_range[1] / 5) * 5, ceiling(year_range[2] / 5) * 5, by = 5)
  
  # ---- Panel C: AUC by year — legend inside top-right ----
  pC <- ggplot2::ggplot(yearly_long, ggplot2::aes(x = year, y = auc, colour = model)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.3, shape = 16) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "C", fontface = "bold",
                      size = 3.5, hjust = -0.3, vjust = 1.3, family = "Arial") +
    ggplot2::scale_color_manual(values = col_vals, labels = unname(model_labels[line_models]), name = NULL) +
    ggplot2::scale_x_continuous(breaks = yr_breaks) +
    ggplot2::scale_y_continuous(limits = c(0.6, 1.0)) +
    ggplot2::labs(x = NULL, y = "Area under ROC curve") +
    line_theme() +
    ggplot2::theme(legend.position = "none")
  
  # ---- Panel D: Brier skill by year — no legend ----
  pD <- ggplot2::ggplot(yearly_long, ggplot2::aes(x = year, y = brier_skill, colour = model)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "#000000") +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.3, shape = 16) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "D", fontface = "bold",
                      size = 3.5, hjust = -0.3, vjust = 1.3, family = "Arial") +
    ggplot2::scale_color_manual(values = col_vals, labels = unname(model_labels[line_models]), name = NULL) +
    ggplot2::scale_x_continuous(breaks = yr_breaks) +
    ggplot2::labs(x = "Year", y = "Brier skill score") +
    line_theme() +
    ggplot2::theme(legend.position = "none")
  
  # ---- Compose with patchwork ----
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for fig_3 combined plots. Install with: install.packages('patchwork')")
  }
  
  top_row <- (pA | pB)
  full    <- (top_row / pC / pD) +
    patchwork::plot_layout(heights = c(3, 2, 2.5))
  
  full
}


