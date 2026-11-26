#!/usr/bin/env Rscript
# A figure file from any analysis run, without opening the app.
#
# Each tab's figure is built by the same function the app calls -- tvc_plot, mc_ggplot, gi_ggplot,
# ess_*_plotly and so on -- reading the run's own tables, and is then rendered through a headless browser by
# the same Plotly.toImage the app's download button uses. So a figure produced here is the figure the app
# shows: there is one definition of each plot in this repository, not a print copy that can drift from it.
#
# Usage:
#   Rscript scripts/local/ytab_figure.R TAB RUN_DIR PROJECT_CONFIG [options]
#
#   TAB            essentiality | treated-vs-control | multi-condition | genetic-interactions |
#                  normalization | cross-species
#   --plot NAME    which plot of that tab (see --list)
#   --output FILE  .svg (default) or .png; defaults to RUN_DIR/figures/
#   --list         print the plots the tab offers and exit
#
#   Essentiality only:
#     --view hmm|hybrid|compare      which model the results are read through
#     --chromosome NAME              by project name ("Chr A") or accession; repeatable.
#                                    None = genome overview, one = detailed track, several = compared.
#   Multi-condition only:
#     --pair NAME                    which pair of conditions
#   Rendering:
#     --width N --height N  pixels (1600 x 900)      --scale N  raster multiplier for .png (3)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: ytab_figure.R TAB RUN_DIR PROJECT_CONFIG [--plot NAME] [--output FILE]", call. = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/", mustWork = TRUE)

one <- function(flag, default = NULL) { i <- which(args == flag); if (!length(i)) default else args[[i[[1L]] + 1L]] }
many <- function(flag) { i <- which(args == flag); if (!length(i)) character(0) else args[i + 1L] }
flag <- function(name) name %in% args

TABS <- c("essentiality", "treated-vs-control", "multi-condition", "genetic-interactions", "normalization", "cross-species")
PLOTS <- list(
  essentiality = c("track", "scores", "calls", "target", "roc"),
  `treated-vs-control` = c("volcano", "ma", "sites", "comet", "nullz"),
  `multi-condition` = c("volcano", "correlation", "heatmap"),
  `genetic-interactions` = c("scatter", "forest", "volcano"),
  normalization = c("pairs", "pca", "sites"),
  `cross-species` = c("sets", "agreement"))

tab <- args[[1L]]
if (!tab %in% TABS) stop("Unknown tab '", tab, "'. One of: ", paste(TABS, collapse = ", "), call. = FALSE)
if (flag("--list")) { cat(paste(PLOTS[[tab]], collapse = "\n"), "\n", sep = ""); quit(save = "no") }
if (length(args) < 3L) stop("Usage: ytab_figure.R ", tab, " RUN_DIR PROJECT_CONFIG [options]", call. = FALSE)
run_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
config  <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
plot_name <- one("--plot", PLOTS[[tab]][[1L]])
if (!plot_name %in% PLOTS[[tab]])
  stop("'", plot_name, "' is not a plot of ", tab, ". One of: ", paste(PLOTS[[tab]], collapse = ", "), call. = FALSE)
width <- as.numeric(one("--width", 1600)); height <- as.numeric(one("--height", 900)); scale <- as.numeric(one("--scale", 3))

suppressPackageStartupMessages({library(shiny); library(bslib); library(ggplot2); library(plotly)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
# ytab_read_manifest needs this, but it lives in job_manager.R, which starts a session poller when sourced and
# so blocks a plain script. One small function is cheaper than loading that file.
read_json_safely <- function(path) {
  if (!nzchar(path %||% "") || !file.exists(path)) return(NULL)
  tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
}
# Named rather than sourcing the whole directory: some of the app's files expect a running session and block
# when loaded outside one, which hung this script before the figure was ever built.
for (name in c("project_discovery.R", "ui_helpers.R", "ui_components.R", "analysis_runs.R", "plot_capture.R",
               "qc_plot_utils.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "glabrata_annotation_lookup.R", "reference_projects.R",
               "essentiality_state.R", "essentiality_module.R", "treated_vs_control.R", "comet.R",
               "comparisons_module.R", "normalization_module.R", "cross_species.R"))
  suppressWarnings(try(source(file.path(repo_root, "app/shiny/R", name)), silent = TRUE))
if (!ytab_capture_available())
  stop("Writing a figure needs the chromote package and a Chrome or Chromium installation.", call. = FALSE)

project <- read_project_summary(config, repo_root)
read_table <- function(...) ytab_read_analysis_tsv(file.path(run_dir, ...))
display <- ytab_plot_customization_values(list(), "cli")

figure <- switch(tab,

  essentiality = {
    windows <- read_table("plots", "windows.tsv.gz"); genes <- read_table("models", "genes.tsv")
    if (!nrow(windows)) stop("No windows table in ", run_dir, call. = FALSE)
    view <- match.arg(one("--view", "hmm"), c("hmm", "hybrid", "compare"))
    chroms <- unique(windows$chromosome); names_all <- qc_plot_chromosome_display(chroms, project)
    wanted <- many("--chromosome")
    chosen <- if (!length(wanted)) character(0) else {
      hit <- chroms[match(wanted, chroms)]; gap <- is.na(hit)
      hit[gap] <- chroms[match(wanted[gap], names_all)]
      if (anyNA(hit)) stop("Unknown chromosome: ", paste(wanted[is.na(hit)], collapse = ", "),
                           ". This run has: ", paste(names_all, collapse = ", "), call. = FALSE)
      hit
    }
    labels <- if (length(chosen)) qc_plot_chromosome_display(chosen, project) else names_all
    if (identical(plot_name, "track")) {
      if (!length(chosen)) ess_genome_overview_plotly(windows, genes, chroms, names_all, view)
      else if (length(chosen) > 1L) ess_multi_chromosome_plotly(windows, genes, chosen, labels, view)
      else if (identical(view, "hybrid")) ess_gene_track_plotly(windows, genes, chosen, labels)
      else ess_track_plotly(windows, chosen, labels)
    } else {
      manifest <- ytab_read_manifest(run_dir)
      ytab_ggplotly(ess_ggplot(plot_name, genes, manifest = manifest, run_dir = run_dir, view = view,
                               project_label = project_display_label(project)))
    }
  },

  `treated-vs-control` = {
    genes <- read_table("genes.tsv"); manifest <- ytab_read_manifest(run_dir)
    if (!nrow(genes)) stop("No genes table in ", run_dir, call. = FALSE)
    genes$Name <- ytab_display_gene_names(genes$ORF, genes$Name, repo_root)
    tvc_plotly(tvc_plot(plot_name, genes, manifest$method %||% "resampling", display,
                        list(unit = manifest$value_unit %||% "", label_n = as.integer(one("--label-n", 10)),
                             direction = one("--direction", "both"))))
  },

  `multi-condition` = {
    # The same files mc_data() reads in the app, under tests/ and plots/.
    label <- function(d) { if (nrow(d) && "ORF" %in% names(d)) d$Name <- ytab_display_gene_names(d$ORF, d$Name, repo_root); d }
    data <- list(manifest = ytab_read_manifest(run_dir),
                 anova = label(read_table("tests", "anova_genes.tsv")),
                 pairs = label(read_table("tests", "pairs_genes.tsv")),
                 correlation = read_table("plots", "response_correlation.tsv"),
                 heatmap = label(read_table("plots", "heatmap_matrix.tsv")))
    if (!nrow(data$pairs)) stop("No pairs table in ", run_dir, call. = FALSE)
    pair <- one("--pair", unique(data$pairs$pair)[[1L]])
    if (!pair %in% data$pairs$pair) stop("Unknown pair '", pair, "'. This run has: ",
                                         paste(unique(data$pairs$pair), collapse = ", "), call. = FALSE)
    tvc_plotly(mc_ggplot(plot_name, data, pair = pair, s = display))
  },

  `genetic-interactions` = {
    genes <- read_table("genes.tsv"); manifest <- ytab_read_manifest(run_dir)
    if (!nrow(genes)) stop("No genes table in ", run_dir, call. = FALSE)
    tvc_plotly(gi_ggplot(plot_name, genes, manifest$reference_strain %||% "A", manifest$query_strain %||% "B", s = display))
  },

  normalization = {
    data <- norm_plot_data(run_dir)
    ytab_ggplotly(norm_ggplot(data, choice = plot_name))
  },

  `cross-species` = {
    if (identical(plot_name, "sets")) ytab_ggplotly(xs_sets_ggplot(read_table("intersections.tsv"), project$species))
    else ytab_ggplotly(xs_agreement_ggplot(read_table("pairwise.tsv"), project$species))
  })

where <- if (identical(tab, "essentiality") && length(many("--chromosome")))
  paste0("_", paste(gsub("[^A-Za-z0-9]+", "", qc_plot_chromosome_display(many("--chromosome"), project)), collapse = "-")) else ""
file <- one("--output", file.path(run_dir, "figures",
            sprintf("%s_%s_%s%s.svg", project_display_label(project), gsub("-", "_", tab), plot_name, where)))
dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
ytab_capture_plotly(figure, file, width = width, height = height, scale = scale)
cat(jsonlite::toJSON(list(status = "success", tab = tab, plot = plot_name, figure = file,
                          bytes = file.size(file)), auto_unbox = TRUE), "\n", sep = "")
