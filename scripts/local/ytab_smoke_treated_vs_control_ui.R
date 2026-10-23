#!/usr/bin/env Rscript
# Treated vs Control: method choice in Run and Results, labels of the three methods, resampling/Gumbel result views.
# Usage: Rscript scripts/local/ytab_smoke_treated_vs_control_ui.R REPO_ROOT [RESULT_DIR ...]
# RESULT_DIR: optional `ytab treated-vs-control` outputs whose plots and tables are drawn.
args <- commandArgs(trailingOnly = TRUE); root <- normalizePath(args[[1L]])
suppressPackageStartupMessages({library(shiny); library(bslib)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
for (name in c("job_manager.R", "job_progress.R", "project_discovery.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "analysis_runs.R", "normalization_module.R", "ui_fitness.R", "treated_vs_control.R", "comparisons_module.R"))
  source(file.path(root, "app/shiny/R", name))
html <- htmltools::renderTags(comparisons_module_ui())$html
stopifnot(grepl("Treated vs Control: Run", html, fixed = TRUE), grepl("Treated vs Control: Results", html, fixed = TRUE),
          !grepl("Fitness: Run", html, fixed = TRUE),
          grepl("Local-abundance z-score", html, fixed = TRUE), grepl(">Resampling<", html, fixed = TRUE),
          grepl("Control-conditioned Gumbel (depletion only)", html, fixed = TRUE), !grepl("narrowly failed FDR benchmark)", html, fixed = TRUE),
          grepl("passed its pre-registered calibration gates", html, fixed = TRUE),   # the validation status stays in Methods
          grepl("fitness_normalization_selector", html, fixed = TRUE), grepl("tvc_result_selector", html, fixed = TRUE),
          !grepl("inside the R analysis", html, fixed = TRUE))
# TRANSIT is named only in Methods references.
stopifnot(!grepl("transit", gsub('(?s)<p class="ytab-reference">.*?</p>', "", html, perl = TRUE), ignore.case = TRUE))
stopifnot(identical(tvc_allowed_methods$gumbel, c("ttr", "totreads", "nzmean")), "betageom" %in% tvc_allowed_methods$resampling)
for (run in args[-1L]) {
  m <- ytab_read_manifest(run); genes <- ytab_read_analysis_tsv(file.path(run, "genes.tsv"))
  stopifnot(identical(m$kind, "treated_vs_control"), nrow(genes) > 0)
  p <- tvc_ggplot(genes[genes$group == genes$group[[1L]], , drop = FALSE], m$method)
  built <- ytab_ggplotly(p); stopifnot(length(built$x$data) >= 2L)
  g <- genes[genes$group == genes$group[[1L]], , drop = FALSE]
  sites <- utils::read.delim(file.path(run, g$group[[1L]], "sites.tsv.gz"), stringsAsFactors = FALSE)
  top <- tvc_tested(g); top <- top$ORF[order(top$q_bh)][[1L]]
  gs <- sites[vapply(strsplit(as.character(sites$genes), "[;,]"), function(x) top %in% x, logical(1)), , drop = FALSE]
  stopifnot(nrow(gs) > 0)
  for (choice in tvc_plot_choices(m$method)) {
    q <- tvc_plot(choice, g, m$method, gs, top); stopifnot(inherits(q, "ggplot"))
    invisible(ytab_ggplotly(q))
  }
  f <- tempfile(fileext = ".png"); ytab_save_ggplot(f, p); stopifnot(file.size(f) > 1000)
}
cat("PASS\n")
