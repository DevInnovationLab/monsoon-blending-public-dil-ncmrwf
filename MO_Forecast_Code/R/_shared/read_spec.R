# ---------------------------------------------------------------------------
# read_spec.R
#
# Purpose:
#   Read YAML “spec” files and provide lightweight validation utilities.
#   This file centralizes:
#     - locating and loading a spec by id/type
#     - validating presence of required fields (top-level or nested paths)
#     - small helpers that derive labels/column names from spec components
#
# Function summary:
#   get_arg(key, default = NULL)
#     - Retrieves a command-line argument of the form key=value from commandArgs().
#
#   load_spec(spec_id, type)
#     - Loads YAML from specs/<type>/<spec_id>.yml.
#
#   validate_spec(spec,
#                 required_top   = character(),
#                 required_paths = character(),
#                 checks         = list(),
#                 label          = "spec")
#     - Generic validator:
#         * required_top: top-level keys that must exist
#         * required_paths: dotted paths (e.g. "input.nc_folder") that must exist
#         * checks: list of functions(spec) run after structural checks
#     - Returns spec unchanged on success; errors with clear message on failure.
#
#   outcome_levels_from_spec(spec)
#     - Builds ordered outcome labels: day_1..day_max plus spec$targets$plus_label.
#
#   clim_cols_from_spec(spec)
#     - Constructs climatology column names from spec$climatology and target max_day.
#     - Important: “plus” column is last (downstream code depends on this ordering).
#
# ---------------------------------------------------------------------------



library(yaml)
library(dplyr)
library(purrr)

get_arg <- function(key, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

load_spec <- function(spec_id, type) {
  path <- file.path("specs",type,  paste0(spec_id, ".yml"))
  if (!file.exists(path)) {
    stop("Spec file not found: ", path)
  }
  yaml::read_yaml(path)
}

# validate_spec <- function(spec) {
#   required_top <- c("id", "data", "targets", "forecasts",
#                     "smoothing", "design", "stan_model", "eval")
#   missing_top <- setdiff(required_top, names(spec))
#   if (length(missing_top) > 0L) {
#     stop("Spec missing top-level fields: ", paste(missing_top, collapse = ", "))
#   }
#   
#   if (is.null(spec$data$path)) {
#     stop("spec$data$path is required.")
#   }
#   if (is.null(spec$data$year_min) || is.null(spec$data$year_max)) {
#     stop("spec$data$year_min and year_max are required.")
#   }
#   
#   spec
# }
# 
# validate spec: ensure that the spec contains required variables for the script 

validate_spec <- function(spec,
                          required_top = character(),
                          required_paths = character(),
                          checks = list(),
                          label = "spec") {
  
  # --- helpers ---
  .stop_missing <- function(msg) stop(msg, call. = FALSE)
  
  .has_path <- function(x, path) {
    parts <- strsplit(path, "\\.", fixed = FALSE)[[1]]
    cur <- x
    for (p in parts) {
      if (is.null(cur) || !is.list(cur) || is.null(cur[[p]])) return(FALSE)
      cur <- cur[[p]]
    }
    TRUE
  }
  
  .get_path <- function(x, path) {
    parts <- strsplit(path, "\\.", fixed = FALSE)[[1]]
    cur <- x
    for (p in parts) cur <- cur[[p]]
    cur
  }
  
  # --- 1) top-level fields ---
  if (length(required_top) > 0L) {
    miss <- setdiff(required_top, names(spec))
    if (length(miss) > 0L) {
      .stop_missing(paste0(label, " missing top-level fields: ", paste(miss, collapse = ", ")))
    }
  }
  
  # --- 2) nested via dotted paths: c("data.path", "data.year_min") ---
  if (length(required_paths) > 0L) {
    miss <- required_paths[!vapply(required_paths, function(p) .has_path(spec, p), logical(1))]
    if (length(miss) > 0L) {
      .stop_missing(paste0(label, " missing required fields: ", paste(miss, collapse = ", ")))
    }
  }
  
  # --- 3) extra checks (functions) ---
  if (length(checks) > 0L) {
    if (!is.list(checks) || any(!vapply(checks, is.function, logical(1)))) {
      .stop_missing("checks must be a list of functions, e.g. list(function(spec) {...}).")
    }
    for (fn in checks) fn(spec)
  }
  
  spec
}


outcome_levels_from_spec <- function(spec) {
  max_day <- spec$targets$max_day
  plus_label <- spec$targets$plus_label
  c(paste0("day_", seq_len(max_day)), plus_label)
}


# Important: "plus" must be the last column (downstream code depends on this ordering)
clim_cols_from_spec <- function(spec) {
  max_day <- spec$targets$max_day
  pref    <- spec$climatology$prefix
  plus    <- spec$climatology$plus_col
  
  cols_1_max <- paste0(pref, 1:max_day)
  c(cols_1_max, plus)
}