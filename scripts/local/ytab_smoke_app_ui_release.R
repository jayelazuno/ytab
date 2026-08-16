#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
repo_root <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

must_exist <- function(path) {
  full <- file.path(repo_root, path)
  if (!file.exists(full)) fail(paste("Missing expected file:", path))
  full
}

must_contain <- function(path, pattern, label = pattern) {
  text <- paste(readLines(file.path(repo_root, path), warn = FALSE), collapse = "\n")
  if (!grepl(pattern, text, perl = TRUE)) fail(paste("Missing", label, "in", path))
  invisible(TRUE)
}

must_not_contain <- function(paths, pattern, label = pattern) {
  for (path in paths) {
    if (!file.exists(path) || dir.exists(path)) next
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    if (grepl(pattern, text, perl = TRUE)) fail(paste("Unexpected", label, "in", path))
  }
}

library(shiny)
library(bslib)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
human_file_size <- function(bytes) {
  bytes <- suppressWarnings(as.numeric(bytes))
  ifelse(is.na(bytes), "", paste0(round(bytes), " B"))
}

invisible(parse(file = must_exist("app/shiny/app.R")))

helper_files <- c(
  "app/shiny/R/ui_components.R",
  "app/shiny/R/plot_display_helpers.R",
  "app/shiny/R/table_display_helpers.R",
  "app/shiny/R/plot_customization_helpers.R"
)
for (helper in helper_files) {
  full <- must_exist(helper)
  invisible(parse(file = full))
  source(full, local = TRUE)
}

invisible(must_exist("app/shiny/www/ytab_release_ui.css"))
must_contain("app/shiny/app.R", "ytab_release_ui\\.css", "release CSS reference")
must_contain("app/shiny/app.R", "ui_components\\.R", "shared UI helper source")
must_contain("app/shiny/app.R", "plot_customization_helpers\\.R", "plot customization helper source")

for (fn in c(
  "ytab_page_header", "ytab_result_card", "ytab_metric_card",
  "ytab_plot_card", "ytab_static_image_card", "ytab_table_card",
  "ytab_technical_details", "ytab_empty_state",
  "ytab_plot_output", "ytab_static_image_ui", "ytab_generated_file_gallery",
  "ytab_static_image_output_card", "ytab_plot_frame", "ytab_two_column_layout",
  "ytab_datatable", "ytab_file_table", "ytab_plot_customization_controls",
  "ytab_plot_customization_values", "ytab_plot_height_px",
  "ytab_with_plot_display_options"
)) {
  if (!exists(fn, mode = "function")) fail(paste("Missing helper function:", fn))
}

must_contain("app/shiny/R/plot_customization_helpers.R", "Display-only controls", "display-only control note")
must_contain("app/shiny/R/plot_customization_helpers.R", "plot_width", "plot width control")
must_contain("app/shiny/R/plot_customization_helpers.R", "label_angle", "label angle control")
must_contain("app/shiny/www/ytab_release_ui.css", "height: auto !important", "aspect-ratio preserving image CSS")
must_contain("app/shiny/www/ytab_release_ui.css", "ytab-two-column-layout", "two-column layout CSS")
must_contain("app/shiny/www/ytab_release_ui.css", "ytab-workspace-header\\.ytab-release-header", "readable header contrast CSS")
must_contain("app/shiny/www/ytab_release_ui.css", "ytab-diagnostic-gallery-grid", "diagnostic gallery density CSS")

workspace_text <- paste(readLines(must_exist("app/shiny/R/ui_workspace.R"), warn = FALSE), collapse = "\n")
for (tab in c("Quality Control", "Essentiality", "Fitness Screen",
              "Gene & Domain Insertion Explorer", "Comparative Species",
              "Results & Exports")) {
  if (!grepl(tab, workspace_text, fixed = TRUE)) fail(paste("Missing major tab:", tab))
}

must_contain("app/shiny/R/ui_qc.R", "summary_qc_plot_choice", "Summary QC plot selector")
must_contain("app/shiny/R/ui_qc.R", "ytab_plot_customization_controls\\(\"summary_qc\"", "Summary QC display controls")
must_contain("app/shiny/R/ui_qc.R", "ytab_two_column_layout", "Summary QC two-column layout")
must_contain("app/shiny/R/ui_essentiality.R", "ytab_two_column_layout", "Essentiality two-column layout")
must_contain("app/shiny/R/ui_essentiality.R", "ytab_plot_customization_controls\\(\"essentiality\"", "Essentiality display controls")
must_contain("app/shiny/R/ui_fitness.R", "fitness_visualization_selector", "Fitness plot selector")
must_contain("app/shiny/R/ui_fitness.R", "ytab_plot_customization_controls\\(\"fitness\"", "Fitness display controls")
must_contain("app/shiny/R/ui_fitness.R", "ytab_two_column_layout", "Fitness two-column layout")
must_contain("app/shiny/app.R", "ytab-static-image-card", "Fitness generated plot static card")
must_contain("app/shiny/app.R", "ytab_plot_frame", "App-rendered plot frame")
must_contain("app/shiny/R/ui_gene_domain_explorer.R", "gene_domain_track_preset", "Gene Explorer track preset")
must_contain("app/shiny/R/ui_gene_domain_explorer.R", "gene_domain_display_height", "Gene Explorer display height")
must_contain("app/shiny/R/gene_domain_explorer_server.R", "ytab_static_image_output_card", "Gene Explorer static image card")
must_contain("app/shiny/R/gene_domain_explorer_state.R", "gene_domain_preset_track_rows", "Gene Explorer preset logic")
must_contain("app/shiny/R/gene_domain_explorer_server.R", "Gene Explorer could not generate a figure", "Gene Explorer friendly failure state")
must_contain("app/shiny/R/gene_domain_explorer_server.R", "Search for a gene, choose tracks, then click Generate figure", "Gene Explorer friendly empty state")
must_contain("app/shiny/R/ui_comparative.R", "comparative_single_species", "Comparative species controls")
must_contain("app/shiny/R/ui_comparative.R", "comparative_ortholog_table", "Comparative orthology table")
must_contain("app/shiny/R/plot_display_helpers.R", "ytab_generated_file_gallery", "generated file gallery helper")
must_contain("app/shiny/R/fitness_generated_plots.R", "ytab-static-image-card", "Fitness generated plot static card")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "ytab-static-image-card", "Diagnostics generated plot static card")
must_contain("app/shiny/app.R", "Show path / filename", "Diagnostic gallery hides paths by default")
must_contain("app/shiny/R/qc_summary_library_plots.R", "qc_plot_bar_horizontal", "Summary QC long-label handling")
must_contain("docs/ui_reference_research_step14F1.md", "Runtime dependency confirmation", "research manifest")

runtime_files <- file.path(repo_root, c(
  "app/shiny/app.R",
  "app/shiny/R/ui_workspace.R",
  "app/shiny/R/ui_qc.R",
  "app/shiny/R/ui_fitness.R",
  "app/shiny/R/ui_gene_domain_explorer.R",
  "app/shiny/R/gene_domain_explorer_server.R",
  "app/shiny/R/ui_comparative.R",
  helper_files
))
runtime_files <- runtime_files[file.exists(runtime_files)]
reference_name <- paste0("RNA", "cross")
external_reference_dir <- paste0("codex/", "ui_", "references")
blocked <- paste0("codex/", reference_name, "|source\\(.*", reference_name,
                  "|readRDS\\(.*codex/", reference_name, "|", external_reference_dir)
must_not_contain(runtime_files, blocked, "runtime dependency on external UI reference")

cat("PASS\n")
