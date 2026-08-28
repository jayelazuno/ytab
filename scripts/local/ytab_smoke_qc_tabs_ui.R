#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: ytab_smoke_qc_tabs_ui.R <project.yaml>", call. = FALSE)
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1L]]),
                      winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/", mustWork = TRUE)

library(shiny)
library(bslib)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left

source(file.path(repo_root, "app/shiny/R/project_discovery.R"))
source(file.path(repo_root, "app/shiny/R/ui_helpers.R"))
source(file.path(repo_root, "app/shiny/R/ui_components.R"))
source(file.path(repo_root, "app/shiny/R/plot_customization_helpers.R"))
source(file.path(repo_root, "app/shiny/R/plot_display_helpers.R"))
source(file.path(repo_root, "app/shiny/R/table_display_helpers.R"))
source(file.path(repo_root, "app/shiny/R/qc_result_state.R"))
source(file.path(repo_root, "app/shiny/R/qc_plot_utils.R"))
source(file.path(repo_root, "app/shiny/R/qc_summary_library_plots.R"))
source(file.path(repo_root, "app/shiny/R/qc_library_diagnostics_plots.R"))
source(file.path(repo_root, "app/shiny/R/ui_qc.R"))

project <- read_project_summary(args[[1L]], repo_root)

render_plot <- function(input, id, expr, width = 1700, height = 1100) {
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = width, height = height, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  ytab_with_plot_display_options(input, id, expr)
  grDevices::dev.off()
  on.exit(NULL)
  stopifnot(file.exists(file), file.info(file)$size > 1000)
  invisible(file)
}

summary_data <- qc_summary_stats_plot_data(project)
stopifnot(nrow(summary_data) >= 4L)
summary_labels <- as.character(summary_data$sample)
stopifnot(any(nchar(summary_labels) >= 18L))

summary_horizontal <- list(
  summary_qc_text_size = "medium",
  summary_qc_label_mode = "full",
  summary_qc_label_angle = "0",
  summary_qc_grid = "show",
  summary_qc_bar_orientation = "horizontal",
  summary_qc_show_value_labels = FALSE
)
summary_vertical <- list(
  summary_qc_text_size = "medium",
  summary_qc_label_mode = "full",
  summary_qc_label_angle = "90",
  summary_qc_grid = "hide",
  summary_qc_bar_orientation = "vertical",
  summary_qc_show_value_labels = TRUE
)

render_plot(summary_horizontal, "summary_qc", plot_qc_summary_metric(project, "complexity"))
render_plot(summary_vertical, "summary_qc", plot_qc_summary_metric(project, "features"))
render_plot(summary_horizontal, "summary_qc", plot_qc_summary_metric(project, "feature_intergenic"))
if (nrow(qc_summary_pairwise_plot_data(project))) {
  render_plot(summary_horizontal, "summary_qc", plot_qc_summary_pairwise_correlations(project))
}

library_summary <- qc_library_summary_plot_data(project)
stopifnot(nrow(library_summary) >= 4L)

library_horizontal <- list(
  library_diagnostics_text_size = "medium",
  library_diagnostics_label_mode = "full",
  library_diagnostics_label_angle = "0",
  library_diagnostics_grid = "show",
  library_diagnostics_bar_orientation = "horizontal",
  library_diagnostics_show_value_labels = FALSE
)
library_vertical <- list(
  library_diagnostics_text_size = "small",
  library_diagnostics_label_mode = "full",
  library_diagnostics_label_angle = "90",
  library_diagnostics_grid = "hide",
  library_diagnostics_bar_orientation = "vertical",
  library_diagnostics_show_value_labels = TRUE
)

if (nrow(qc_library_midlc_plot_data(project))) {
  render_plot(library_horizontal, "library_diagnostics", plot_qc_library_midlc(project))
}
render_plot(library_horizontal, "library_diagnostics", plot_qc_library_jackpot_depth(project))
if (nrow(qc_library_sequence_bias_plot_data(project))) {
  render_plot(library_vertical, "library_diagnostics", plot_qc_library_sequence_bias(project))
}

summary_text <- paste(readLines(file.path(repo_root, "app/shiny/R/qc_summary_library_plots.R"), warn = FALSE),
                      collapse = "\n")
library_text <- paste(readLines(file.path(repo_root, "app/shiny/R/qc_library_diagnostics_plots.R"), warn = FALSE),
                      collapse = "\n")
ui_qc_text <- paste(readLines(file.path(repo_root, "app/shiny/R/ui_qc.R"), warn = FALSE),
                    collapse = "\n")
plot_controls_text <- paste(readLines(file.path(repo_root, "app/shiny/R/plot_customization_helpers.R"), warn = FALSE),
                            collapse = "\n")
app_text <- paste(readLines(file.path(repo_root, "app/shiny/app.R"), warn = FALSE),
                  collapse = "\n")

stopifnot(
  grepl("qc_plot_display_labels", summary_text, fixed = TRUE),
  grepl("qc_plot_label_margin_lines", summary_text, fixed = TRUE),
  grepl("qc_plot_draw_vertical_labels", summary_text, fixed = TRUE),
  grepl("qc_plot_draw_horizontal_labels", summary_text, fixed = TRUE),
  grepl("qc_plot_grid_enabled", summary_text, fixed = TRUE),
  grepl("qc_plot_show_value_labels", summary_text, fixed = TRUE),
  grepl("qc_plot_metric_key_row", summary_text, fixed = TRUE),
  grepl("qc_plot_metric_key_row", library_text, fixed = TRUE),
  grepl("qc_plot_begin_key_layout", library_text, fixed = TRUE),
  grepl("summary_qc_plot_controls", ui_qc_text, fixed = TRUE),
  grepl("ytab_plot_customization_controls\\(\"summary_qc\"", app_text),
  !grepl("default_bar_orientation = \"horizontal\"", ui_qc_text, fixed = TRUE),
  grepl("default_bar_orientation = \"vertical\"", plot_controls_text, fixed = TRUE),
  grepl("default_label_angle = \"30\"", plot_controls_text, fixed = TRUE),
  grepl("library_diagnostics_plot_controls", ui_qc_text, fixed = TRUE),
  grepl("ytab_plot_customization_controls\\(\"library_diagnostics\"", app_text),
  grepl("diagnostic_show_archived_static", ui_qc_text, fixed = TRUE),
  grepl("diagnostic_files_table", ui_qc_text, fixed = TRUE),
  grepl("!isTRUE\\(input\\$diagnostic_show_archived_static\\)", app_text),
  !grepl("diagnostic_plot_cards", ui_qc_text, fixed = TRUE),
  !grepl("diagnostic_gallery_page_size", app_text, fixed = TRUE),
  !grepl("App-rendered: Summary QC", ui_qc_text, fixed = TRUE)
)

cat("PASS\n")
