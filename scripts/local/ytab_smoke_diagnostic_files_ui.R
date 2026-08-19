#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)

all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1]]),
                      winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/", mustWork = TRUE)

library(shiny)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
panel_card <- function(title, ...) tags$section(tags$h3(title), ...)
ytab_two_column_layout <- function(controls, main) tags$div(controls, main)
ytab_control_panel <- function(title, ...) tags$aside(tags$h4(title), ...)
ytab_plot_card <- function(title, content, description = NULL, ...) tags$section(tags$h4(title), content)
job_progress_ui <- function(id, compact = FALSE) tags$div(id = id)
ytab_plot_customization_controls <- function(...) tags$div()

source(file.path(root, "app/shiny/R/project_discovery.R"), local = TRUE)
source(file.path(root, "app/shiny/R/qc_plot_utils.R"), local = TRUE)
source(file.path(root, "app/shiny/R/ui_qc.R"), local = TRUE)

project <- read_project_summary(args[[1]], root)
inventory <- build_diagnostic_file_inventory(project)

stopifnot(
  nrow(inventory) > 0L,
  all(c("file", "filename", "extension", "display_type", "plot_type", "sample",
        "sample_set", "run_id", "size_bytes", "size_display", "modified",
        "relative_path", "viewable", "is_plot") %in% names(inventory)),
  !any(inventory$extension %in% c("csv", "tsv") & inventory$is_plot),
  diagnostic_plot_title("sample.centromere_bias.png") == "Centromere bias",
  diagnostic_plot_title("sample.metaplots.png") == "Feature metaplots",
  grepl("KB|MB|B", human_file_size(319 * 1024))
)

plots <- inventory[inventory$is_plot, , drop = FALSE]
tables <- inventory[!inventory$is_plot, , drop = FALSE]
filtered <- filter_diagnostic_inventory(inventory, sample = if (nrow(inventory)) inventory$sample[[1]] else "All")
default_visible <- inventory[!(inventory$is_plot & inventory$extension %in% c("png", "jpg", "jpeg", "svg")), , drop = FALSE]

stopifnot(
  nrow(filtered) <= nrow(inventory),
  nrow(default_visible) == nrow(tables),
  nrow(default_visible) <= nrow(inventory)
)

ui <- paste(readLines(file.path(root, "app/shiny/R/ui_qc.R"), warn = FALSE), collapse = "\n")
stopifnot(
  grepl("diagnostic_files_table", ui, fixed = TRUE),
  grepl("diagnostic_show_archived_static", ui, fixed = TRUE),
  !grepl("diagnostic_view_mode", ui, fixed = TRUE),
  !grepl("Plot gallery", ui, fixed = TRUE),
  !grepl("diagnostic_plot_cards", ui, fixed = TRUE),
  !grepl("base64", ui, ignore.case = TRUE)
)

cat(sprintf("Indexed files: %d; archived/static images indexed but hidden by default: %d\n",
            nrow(inventory), nrow(plots)))
cat("PASS\n")
