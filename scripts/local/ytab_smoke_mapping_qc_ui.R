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
stopifnot(any(nchar(labels) >= 18L))
stopifnot("total_records" %in% names(data))
stopifnot("primary_mapped" %in% names(data))
stopifnot("percent_mapped" %in% names(data))
stopifnot("avg_mapq_mapped_primary" %in% names(data))

default_input <- list(
  mapping_qc_plot_choice = "read_counts",
  mapping_qc_text_size = "medium",
  mapping_qc_label_mode = "full",
  mapping_qc_label_angle = "0",
  mapping_qc_grid = "show",
  mapping_qc_bar_orientation = "horizontal",
  mapping_qc_show_value_labels = FALSE
)
compact_input <- list(
  mapping_qc_plot_choice = "percent_mapped",
  mapping_qc_text_size = "small",
  mapping_qc_label_mode = "compact",
  mapping_qc_label_angle = "90",
  mapping_qc_grid = "hide",
  mapping_qc_bar_orientation = "vertical",
  mapping_qc_show_value_labels = TRUE
)
full_vertical_input <- list(
  mapping_qc_plot_choice = "mapq",
  mapping_qc_text_size = "medium",
  mapping_qc_label_mode = "full",
  mapping_qc_label_angle = "90",
  mapping_qc_grid = "show",
  mapping_qc_bar_orientation = "vertical",
  mapping_qc_show_value_labels = FALSE
)

render_plot <- function(input) {
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 1600, height = 1000, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  ytab_with_plot_display_options(input, "mapping_qc",
                                 plot_qc_mapping_stats(project, input$mapping_qc_plot_choice))
  grDevices::dev.off()
  on.exit(NULL)
  file
}

default_png <- render_plot(default_input)
compact_png <- render_plot(compact_input)
full_vertical_png <- render_plot(full_vertical_input)
stopifnot(file.exists(default_png), file.info(default_png)$size > 1000)
stopifnot(file.exists(compact_png), file.info(compact_png)$size > 1000)
stopifnot(file.exists(full_vertical_png), file.info(full_vertical_png)$size > 1000)

options(ytab.plot.label_mode = "full")
stopifnot(identical(qc_plot_display_labels(labels[[1L]]), labels[[1L]]))
options(ytab.plot.label_mode = "compact")
stopifnot(!identical(qc_plot_display_labels(labels[[1L]]), ""))
options(ytab.plot.text_size = "small")
small_sizes <- qc_plot_text_sizes()
small_left_margin <- qc_plot_label_margin_lines(labels, orientation = "horizontal")
small_bottom_margin <- qc_plot_label_margin_lines(labels, angle = 90L, orientation = "vertical")
options(ytab.plot.text_size = "medium")
medium_sizes <- qc_plot_text_sizes()
options(ytab.plot.text_size = "large")
large_sizes <- qc_plot_text_sizes()
large_left_margin <- qc_plot_label_margin_lines(labels, orientation = "horizontal")
large_bottom_margin <- qc_plot_label_margin_lines(labels, angle = 90L, orientation = "vertical")
stopifnot(
  medium_sizes$axis > small_sizes$axis,
  large_sizes$axis > medium_sizes$axis,
  medium_sizes$sample > small_sizes$sample,
  large_sizes$sample > medium_sizes$sample,
  medium_sizes$lab > small_sizes$lab,
  large_sizes$lab > medium_sizes$lab,
  large_left_margin >= small_left_margin,
  large_bottom_margin >= small_bottom_margin
)

mapping_plot_text <- paste(readLines(file.path(repo_root, "app/shiny/R/qc_mapping_stats_plot.R"), warn = FALSE),
                           collapse = "\n")
qc_utils_text <- paste(readLines(file.path(repo_root, "app/shiny/R/qc_plot_utils.R"), warn = FALSE),
                       collapse = "\n")
ui_qc_text <- paste(readLines(file.path(repo_root, "app/shiny/R/ui_qc.R"), warn = FALSE),
                    collapse = "\n")
stopifnot(
  grepl("read_counts", ui_qc_text, fixed = TRUE),
  grepl("percent_mapped", ui_qc_text, fixed = TRUE),
  grepl("mapq", ui_qc_text, fixed = TRUE),
  grepl("mapping_qc_plot_key", paste(readLines(file.path(repo_root, "app/shiny/app.R"), warn = FALSE),
                                     collapse = "\n"), fixed = TRUE),
  grepl("ytab-inline-legend", paste(readLines(file.path(repo_root, "app/shiny/www/ytab_release_ui.css"),
                                             warn = FALSE), collapse = "\n"), fixed = TRUE),
  grepl("Total reads", mapping_plot_text, fixed = TRUE),
  grepl("Mapped reads", mapping_plot_text, fixed = TRUE),
  grepl("% reads mapped", mapping_plot_text, fixed = TRUE),
  grepl("Average MAPQ mapped primary", mapping_plot_text, fixed = TRUE),
  grepl("primary_mapped", mapping_plot_text, fixed = TRUE),
  grepl("avg_mapq_mapped_primary", mapping_plot_text, fixed = TRUE),
  !grepl("% l", mapping_plot_text, fixed = TRUE),
  !grepl("legend(\"topleft\"", mapping_plot_text, fixed = TRUE),
  !grepl("legend(\"right\"", mapping_plot_text, fixed = TRUE),
  grepl("legend(\"center\"", mapping_plot_text, fixed = TRUE),
  grepl("layout(matrix", mapping_plot_text, fixed = TRUE),
  grepl("qc_plot_label_margin_lines", mapping_plot_text, fixed = TRUE),
  grepl("qc_plot_text_sizes", qc_utils_text, fixed = TRUE),
  grepl("small = 0.95", qc_utils_text, fixed = TRUE),
  grepl("large = 1.55", qc_utils_text, fixed = TRUE),
  grepl("qc_plot_value_cex", mapping_plot_text, fixed = TRUE),
  grepl("qc_plot_key_cex", mapping_plot_text, fixed = TRUE),
  !grepl("App-rendered: Mapping summary", ui_qc_text, fixed = TRUE)
)

cat("PASS\n")
