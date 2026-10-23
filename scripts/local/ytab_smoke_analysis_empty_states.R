#!/usr/bin/env Rscript
# A project with no analysis results yet: every output of the new panels (library correlation, normalization,
# essentiality, multi-condition, genetic interactions) must render its empty state or stay silent, never error.
# (Regression: req(nrow(x)) accepted 0 rows and indexed the first run -> "subscript out of bounds".)
args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
suppressPackageStartupMessages({library(shiny); library(bslib)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
for (name in c("job_manager.R", "project_discovery.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "analysis_runs.R", "essentiality_state.R", "qc_library_correlation.R", "normalization_module.R",
               "essentiality_module.R", "treated_vs_control.R", "comparisons_module.R"))
  source(file.path(root, "app/shiny/R", name))
config <- file.path(root, "output/projects/Zn_toxicity_screen/config/project.yaml")
project <- modifyList(read_project_summary(config, root), list(project_root = tempfile("empty_project_")))
dir.create(file.path(project$project_root, "config"), recursive = TRUE)
jobs <- list(job_is_running = function() FALSE, start_job = function(...) NULL)
active <- reactive(project); path <- function() config; python <- function() "python"
projects <- reactiveVal(data.frame(project_config = config, project_id = "Zn", species = "glabrata", included_count = 8L))
check <- function(server, ids) testServer(server, {
  session$flushReact()
  for (id in ids) {
    value <- tryCatch(output[[id]], shiny.silent.error = function(e) "silent", error = function(e) stop(id, ": ", conditionMessage(e), call. = FALSE))
    stopifnot(!is.null(value))
  }
})
check(function(input, output, session) qc_library_correlation_server(input, output, session, root, active, path, projects, jobs, python),
      c("libcor_summary", "libcor_result_selector", "libcor_pairs_table"))
check(function(input, output, session) normalization_server(input, output, session, root, active, path, jobs, python),
      c("norm_status", "norm_result_selector", "norm_plot_frame", "norm_native_panels", "norm_stats_table"))
check(function(input, output, session) essentiality_module_server(input, output, session, root, active, path, jobs, python),
      c("ess_run_status", "ess_summary", "ess_result_selector", "ess_chromosome_selector", "ess_plot_frame", "ess_genes_table"))
check(function(input, output, session) comparisons_server(input, output, session, root, active, path, jobs, python),
      c("mc_inputs", "mc_summary", "mc_result_selector", "mc_pair_selector", "mc_plot_frame", "mc_table",
        "gi_inputs", "gi_summary", "gi_result_selector", "gi_plot_frame", "gi_table"))
cat("PASS\n")
