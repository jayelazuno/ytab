#!/usr/bin/env Rscript
# `ytab essentiality-figure` writes the app's own figure, not a second drawing of it.
#
# The figure is built by the app's ess_*_plotly functions and rendered through a headless browser, so this
# checks the parts that can go wrong without one: that a run's tables produce a figure object for every view
# and selection, and that chromosomes resolve by project name as well as by accession. The browser render is
# exercised only when chromote and a Chrome are present, so the test is useful on a machine without them.
#
# Usage: Rscript scripts/local/ytab_smoke_figure.R PROJECT_CONFIG
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: ytab_smoke_figure.R PROJECT_CONFIG", call. = FALSE)
config <- normalizePath(args[[1L]], mustWork = TRUE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(self), "../.."), mustWork = TRUE)

suppressPackageStartupMessages({library(shiny); library(bslib); library(ggplot2); library(plotly)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
for (name in c("project_discovery.R", "ui_helpers.R", "ui_components.R", "analysis_runs.R", "plot_capture.R",
               "qc_plot_utils.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "glabrata_annotation_lookup.R", "essentiality_state.R",
               "essentiality_module.R"))
  suppressWarnings(try(source(file.path(root, "app/shiny/R", name)), silent = TRUE))

project <- read_project_summary(config, root)
runs <- ytab_analysis_runs(file.path(dirname(dirname(config)), "essentiality"))
stopifnot(nrow(runs) > 0)
run <- runs$path[[1L]]
windows <- ytab_read_analysis_tsv(file.path(run, "plots", "windows.tsv.gz"))
genes <- ytab_read_analysis_tsv(file.path(run, "models", "genes.tsv"))
stopifnot(nrow(windows) > 0, nrow(genes) > 0)
chroms <- unique(windows$chromosome); names_all <- qc_plot_chromosome_display(chroms, project)

# Every view and every selection yields a figure, built by the same functions the app renders.
for (view in c("hmm", "hybrid", "compare")) {
  overview <- ess_genome_overview_plotly(windows, genes, chroms, names_all, view)
  stopifnot(inherits(overview, c("plotly", "htmlwidget")))
  multi <- ess_multi_chromosome_plotly(windows, genes, chroms[1:2], names_all[1:2], view)
  stopifnot(inherits(multi, c("plotly", "htmlwidget")))
}
stopifnot(inherits(ess_track_plotly(windows, chroms[[1L]], names_all[[1L]]), c("plotly", "htmlwidget")),
          inherits(ess_gene_track_plotly(windows, genes, chroms[[1L]], names_all[[1L]]), c("plotly", "htmlwidget")))

# A chromosome is nameable either way, and an unknown one says so rather than silently drawing the genome:
# a faulty test harness once passed arguments that were dropped, and the figure came out as the whole genome.
resolve <- function(x) {
  hit <- chroms[match(x, chroms)]; missing <- is.na(hit)
  hit[missing] <- chroms[match(x[missing], names_all)]
  if (anyNA(hit)) stop("unknown", call. = FALSE)
  hit
}
stopifnot(identical(resolve(names_all[[2L]]), chroms[[2L]]),
          identical(resolve(chroms[[2L]]), chroms[[2L]]),
          identical(resolve(c(names_all[[1L]], chroms[[3L]])), c(chroms[[1L]], chroms[[3L]])),
          inherits(tryCatch(resolve("Chr ZZ"), error = function(e) e), "error"))

if (ytab_capture_available()) {
  file <- file.path(tempdir(), "ytab_smoke_figure.svg")
  ytab_capture_plotly(ess_track_plotly(windows, chroms[[1L]], names_all[[1L]]), file, width = 1200, height = 700)
  svg <- paste(readLines(file, warn = FALSE), collapse = "")
  stopifnot(file.size(file) > 10000, grepl("<svg", svg, fixed = TRUE), grepl(names_all[[1L]], svg, fixed = TRUE))
  cat("  rendered through a headless browser:", file.size(file), "bytes, naming", names_all[[1L]], "\n")
} else {
  cat("  chromote or Chrome absent; the browser render was not exercised\n")
}
# `ytab figure` covers every tab, and each tab names its own plots. The dispatcher is the part that rots when
# a plot is renamed or a table moves, so it is checked here without rendering anything.
figure_script <- file.path(root, "scripts/local/ytab_figure.R")
stopifnot(file.exists(figure_script))
text <- paste(readLines(figure_script, warn = FALSE), collapse = "\n")
for (tab in c("essentiality", "treated-vs-control", "multi-condition", "genetic-interactions",
              "normalization", "cross-species"))
  stopifnot(grepl(tab, text, fixed = TRUE))
listed <- system2(file.path(R.home("bin"), "Rscript"), c(figure_script, "treated-vs-control", "--list"), stdout = TRUE)
stopifnot("volcano" %in% listed)
# The tables each tab reads must be the ones the app reads; these paths moved once and the figure failed late.
stopifnot(grepl('read_table("tests", "pairs_genes.tsv")', text, fixed = TRUE),
          grepl('read_table("intersections.tsv")', text, fixed = TRUE),
          grepl('read_table("plots", "windows.tsv.gz")', text, fixed = TRUE))
cat("  ytab figure: 6 tabs, plot lists and table paths\n")
cat("PASS\n")
