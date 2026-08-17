#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: ytab_smoke_mapping_qc_ui.R <project.yaml>", call. = FALSE)
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1L]]),
                      winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/", mustWork = TRUE)

library(shiny)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left

source(file.path(repo_root, "app/shiny/R/project_discovery.R"))
source(file.path(repo_root, "app/shiny/R/plot_customization_helpers.R"))
source(file.path(repo_root, "app/shiny/R/qc_plot_utils.R"))
source(file.path(repo_root, "app/shiny/R/ui_qc.R"))
source(file.path(repo_root, "app/shiny/R/qc_mapping_stats_plot.R"))

project <- read_project_summary(args[[1L]], repo_root)
data <- qc_mapping_stats_plot_data(project)
stopifnot(nrow(data) >= 8L)
labels <- as.character(data$sample)
stopifnot(any(grepl("H2O2-treated-facs-pool", labels, fixed = TRUE)))

default_input <- list(
  mapping_qc_text_size = "medium",
  mapping_qc_label_mode = "full",
  mapping_qc_label_angle = "0",
  mapping_qc_grid = "show",
  mapping_qc_bar_orientation = "horizontal",
  mapping_qc_show_value_labels = FALSE
)
compact_input <- list(
  mapping_qc_text_size = "small",
  mapping_qc_label_mode = "compact",
  mapping_qc_label_angle = "45",
  mapping_qc_grid = "hide",
  mapping_qc_bar_orientation = "vertical",
  mapping_qc_show_value_labels = TRUE
)

render_plot <- function(input) {
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 1600, height = 1000, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  ytab_with_plot_display_options(input, "mapping_qc", plot_qc_mapping_stats(project))
  file
}

default_png <- render_plot(default_input)
compact_png <- render_plot(compact_input)
stopifnot(file.exists(default_png), file.info(default_png)$size > 1000)
stopifnot(file.exists(compact_png), file.info(compact_png)$size > 1000)

stopifnot(
  identical(qc_plot_display_labels("yH298-H2O2-treated-facs-pool1"),
            "yH298-H2O2-treated-facs-pool1"),
  !any(grepl("App-rendered: Mapping summary",
             paste(readLines(file.path(repo_root, "app/shiny/R/ui_qc.R"), warn = FALSE),
                   collapse = "\n"),
             fixed = TRUE))
)

cat("PASS\n")

