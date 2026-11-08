#!/usr/bin/env Rscript
# Treated vs Control: method choice in Run and Results, labels of the three methods, resampling/Gumbel result views.
# Usage: Rscript scripts/local/ytab_smoke_treated_vs_control_ui.R REPO_ROOT [RESULT_DIR ...]
# RESULT_DIR: optional `ytab treated-vs-control` outputs whose plots and tables are drawn.
args <- commandArgs(trailingOnly = TRUE); root <- normalizePath(args[[1L]])
suppressPackageStartupMessages({library(shiny); library(bslib)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
for (name in c("job_manager.R", "job_progress.R", "project_discovery.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "analysis_runs.R", "design_support.R", "reference_projects.R", "normalization_module.R", "ui_fitness.R", "treated_vs_control.R", "comparisons_module.R"))
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
# Every plot of every method from the shared layer, with the display options applied.
s <- ytab_plot_customization_values(list(tvc_point_size = 2, tvc_point_opacity = .4, tvc_text_size = "large", tvc_grid = "show"), "tvc")
stopifnot(identical(names(tvc_plot_choices("zscore"))[[1L]], "Volcano"), identical(names(tvc_plot_choices("gumbel"))[[1L]], "Volcano"),
          "runs" %in% tvc_plot_choices("gumbel"), !"pvalues" %in% tvc_plot_choices("zscore"), "nullz" %in% tvc_plot_choices("zscore"))
for (run in args[-1L]) {
  zscore <- file.exists(file.path(run, "zscore_manifest.json"))
  if (zscore) { method <- "zscore"; genes <- tvc_zscore_genes(run) } else {
    m <- ytab_read_manifest(run); method <- m$method; genes <- ytab_read_analysis_tsv(file.path(run, "genes.tsv"))
    stopifnot(identical(m$kind, "treated_vs_control")) }
  stopifnot(nrow(genes) > 0)
  g <- genes[genes$group == genes$group[[1L]], , drop = FALSE]
  tested <- tvc_tested(g); top <- tested$ORF[order(-tvc_significance(tested, method))][[1L]]
  if (zscore) {
    v <- utils::read.csv(file.path(run, "tables", "gene_values.long.csv"), stringsAsFactors = FALSE)
    controls <- strsplit(g$controls[[1L]], ";")[[1L]]; v <- v[v$feature_id == top & v$sample %in% c(controls, strsplit(g$treated[[1L]], ";")[[1L]]), ]
    v$role <- ifelse(v$sample %in% controls, "control", "treated"); gs <- v
    nz <- utils::read.csv(file.path(run, "tables", "control_control_z_scores.csv")); nz <- nz[nz$background == g$background[[1L]], ]
    pr <- utils::read.csv(file.path(run, "tables", "treated_vs_parent.by_pool.log2fc_z.csv"))
    pr <- unique(pr[pr$background == g$background[[1L]], c("parent_sample", "treated_sample")])
  } else {
    sites <- utils::read.delim(file.path(run, g$group[[1L]], "sites.tsv.gz"), stringsAsFactors = FALSE)
    gs <- sites[vapply(strsplit(as.character(sites$genes), "[;,]"), function(x) top %in% x, logical(1)), , drop = FALSE]; nz <- NULL; pr <- NULL
  }
  stopifnot(nrow(gs) > 0)
  ids <- if (zscore) unique(gs$sample) else c(sprintf("control_%d", seq_along(unlist(m$groups[[g$group[[1L]]]]$controls))), sprintf("treatment_%d", seq_along(unlist(m$groups[[g$group[[1L]]]]$treatments))))
  labs <- stats::setNames(ytab_sample_labels(if (zscore) ids else c(unlist(m$groups[[g$group[[1L]]]]$controls), unlist(m$groups[[g$group[[1L]]]]$treatments)), "output/projects/Gale_2020"), ids)
  for (choice in tvc_plot_choices(method)) {
    q <- tvc_plot(choice, g, method, s, list(sites = gs, gene = top, null_z = nz, label_n = 5, label_genes = top, unit = "CPM", library_labels = labs, pairs = pr, pseudocount = .5))
    if (choice == "sites" && zscore) {   # paired dot plot: one point per library, condition means, one line per design pair
      stopifnot(nrow(q$data) == nrow(gs), all(abs(q$data$y - log2(q$data$value + .5)) < 1e-9),
                nrow(Filter(function(l) inherits(l$geom, "GeomSegment"), q$layers)[[1L]]$data) == nrow(pr)) }
    if (choice == "sites") { shown <- unique(q$data$library); stopifnot(!any(grepl("^SRR", shown))) }
    stopifnot(inherits(q, "ggplot")); built <- ytab_ggplotly(q); stopifnot(length(built$x$data) >= 1L)
    if (choice %in% c("volcano", "ma", "scatter")) {
      # downloads: repelled labels; interactive: no text layer, labels as annotations with leader lines
      stopifnot(any(vapply(q$layers, function(l) inherits(l$geom, "GeomTextRepel"), TRUE)), nrow(attr(q, "tvc_labels")) >= 5L)
      qi <- tvc_plot(choice, g, method, s, list(label_n = 5, interactive = TRUE)); w <- tvc_plotly(qi)
      stopifnot(!any(vapply(qi$layers, function(l) inherits(l$geom, "GeomTextRepel"), TRUE)), length(w$x$layout$annotations) >= 5L,
                isTRUE(w$x$layout$annotations[[length(w$x$layout$annotations)]]$showarrow))
      if (method != "gumbel") { dq <- tvc_plot(choice, g, method, s, list(direction = "depleted", label_n = 5, interactive = TRUE))
        stopifnot(!"enriched" %in% as.character(dq$data$call), "other direction" %in% as.character(dq$data$call)) }
    }
  }
  f <- tempfile(fileext = ".png"); ytab_save_ggplot(f, tvc_plot("volcano", g, method, s)); stopifnot(file.size(f) > 1000)
  cat(sprintf("%s: %s plots OK\n", method, paste(tvc_plot_choices(method), collapse = ", ")))
}
cat("PASS\n")
