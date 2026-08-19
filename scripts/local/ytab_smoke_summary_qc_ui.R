#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: ytab_smoke_summary_qc_ui.R <project.yaml>", call. = FALSE)
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
source(file.path(repo_root, "app/shiny/R/qc_plot_utils.R"))
source(file.path(repo_root, "app/shiny/R/qc_summary_library_plots.R"))
source(file.path(repo_root, "app/shiny/R/ui_qc.R"))

project <- read_project_summary(args[[1L]], repo_root)

render_plot <- function(input, expr, width = 1600, height = 1000) {
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = width, height = height, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  ytab_with_plot_display_options(input, "summary_qc", expr)
  grDevices::dev.off()
  on.exit(NULL)
  stopifnot(file.exists(file), file.info(file)$size > 1000)
  invisible(file)
}

data <- qc_summary_combined_features_hit_data(project)
stopifnot(nrow(data) >= 1L)
choices <- qc_summary_combined_feature_group_choices(project)
choice_values <- unname(unlist(choices, use.names = FALSE))
stopifnot(length(choices) >= 1L)
stopifnot("control" %in% choice_values || "treated" %in% choice_values)

has_control <- "control" %in% as.character(data$group)
has_treated <- "treated" %in% as.character(data$group)
if (has_control) stopifnot("control" %in% choice_values)
if (has_treated) stopifnot("treated" %in% choice_values)
if (has_control && has_treated) stopifnot("both" %in% choice_values)

labels <- qc_summary_combined_feature_group_labels(project)
project_id <- as.character(project$project_id %||% "")
if (identical(project_id, "Zn_toxicity_screen")) {
  stopifnot(
    identical(labels$control, "Mock controls combined"),
    identical(labels$treated, "1.5 mM Zn-treated combined")
  )
}
if (identical(project_id, "H2O2_screen_v1")) {
  stopifnot(
    identical(labels$control, "Parents combined"),
    identical(labels$treated, "H2O2-treated combined")
  )
}

vertical_large <- list(
  summary_qc_text_size = "large",
  summary_qc_label_mode = "full",
  summary_qc_label_angle = "90",
  summary_qc_grid = "show",
  summary_qc_bar_orientation = "vertical",
  summary_qc_show_value_labels = TRUE
)
horizontal_small <- list(
  summary_qc_text_size = "small",
  summary_qc_label_mode = "full",
  summary_qc_label_angle = "0",
  summary_qc_grid = "hide",
  summary_qc_bar_orientation = "horizontal",
  summary_qc_show_value_labels = FALSE
)
vertical_45 <- list(
  summary_qc_text_size = "medium",
  summary_qc_label_mode = "compact",
  summary_qc_label_angle = "45",
  summary_qc_grid = "show",
  summary_qc_bar_orientation = "vertical",
  summary_qc_show_value_labels = TRUE
)

if (has_control) render_plot(vertical_large, plot_qc_summary_combined_features_hit(project, "control"))
if (has_treated) render_plot(horizontal_small, plot_qc_summary_combined_features_hit(project, "treated"))
if (has_control && has_treated) {
  render_plot(vertical_45, plot_qc_summary_combined_features_hit(project, "both"))
}

old <- options(ytab.plot.text_size = "small")
small_sizes <- qc_plot_text_sizes()
small_cex <- qc_plot_label_cex(1.05)
sample_labels <- as.character(qc_summary_stats_plot_data(project)$sample)
small_left_margin <- qc_plot_label_margin_lines(sample_labels, orientation = "horizontal")
small_bottom_margin <- qc_plot_label_margin_lines(sample_labels, angle = 90L, orientation = "vertical")
options(ytab.plot.text_size = "medium")
medium_sizes <- qc_plot_text_sizes()
options(ytab.plot.text_size = "large")
large_sizes <- qc_plot_text_sizes()
large_cex <- qc_plot_label_cex(1.05)
large_left_margin <- qc_plot_label_margin_lines(sample_labels, orientation = "horizontal")
large_bottom_margin <- qc_plot_label_margin_lines(sample_labels, angle = 90L, orientation = "vertical")
options(old)
stopifnot(
  large_cex > small_cex,
  medium_sizes$axis > small_sizes$axis,
  large_sizes$axis > medium_sizes$axis,
  medium_sizes$sample > small_sizes$sample,
  large_sizes$sample > medium_sizes$sample,
  medium_sizes$lab > small_sizes$lab,
  large_sizes$lab > medium_sizes$lab,
  large_left_margin >= small_left_margin,
  large_bottom_margin >= small_bottom_margin
)

summary_text <- paste(readLines(file.path(repo_root, "app/shiny/R/qc_summary_library_plots.R"), warn = FALSE),
                      collapse = "\n")
ui_text <- paste(readLines(file.path(repo_root, "app/shiny/R/ui_qc.R"), warn = FALSE),
                 collapse = "\n")
app_text <- paste(readLines(file.path(repo_root, "app/shiny/app.R"), warn = FALSE),
                  collapse = "\n")

stopifnot(
  grepl("qc_summary_combined_feature_group_choices", summary_text, fixed = TRUE),
  grepl("qc_summary_combined_feature_sample_roles", summary_text, fixed = TRUE),
  grepl("control_or_treated", summary_text, fixed = TRUE),
  grepl("library_role", summary_text, fixed = TRUE),
  grepl("rect\\(", summary_text),
  grepl("bar_half <- if \\(n == 1L\\) 0\\.22", summary_text),
  grepl("lwd = max\\(1\\.6, qc_plot_lwd\\)", summary_text),
  grepl("qc_plot_label_cex\\(1\\.05\\)", summary_text),
  grepl("qc_plot_value_cex", paste(readLines(file.path(repo_root, "app/shiny/R/qc_plot_utils.R"), warn = FALSE),
                                   collapse = "\n"), fixed = TRUE),
  grepl("qc_plot_text_sizes", paste(readLines(file.path(repo_root, "app/shiny/R/qc_plot_utils.R"), warn = FALSE),
                                    collapse = "\n"), fixed = TRUE),
  grepl("summary_qc_combined_group_selector", ui_text, fixed = TRUE),
  grepl("summary_qc_combined_group_selector", app_text, fixed = TRUE),
  grepl("input\\$summary_qc_combined_group", app_text),
  grepl("Combined group", app_text, fixed = TRUE)
)

cat("PASS\n")
