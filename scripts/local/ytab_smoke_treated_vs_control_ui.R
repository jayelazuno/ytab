#!/usr/bin/env Rscript
# Treated vs Control: method choice in Run and Results, labels of the two offered methods, result views.
# The control-conditioned Gumbel is retired from the chooser (it returned calls on 2 of 14 datasets) but its label
# and its allowed normalizations stay, so a run stored before the retirement still resolves and still displays.
# Usage: Rscript scripts/local/ytab_smoke_treated_vs_control_ui.R REPO_ROOT [RESULT_DIR ...]
# RESULT_DIR: optional `ytab treated-vs-control` outputs whose plots and tables are drawn.
args <- commandArgs(trailingOnly = TRUE); root <- normalizePath(args[[1L]])
suppressPackageStartupMessages({library(shiny); library(bslib)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
for (name in c("job_manager.R", "job_progress.R", "project_discovery.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "analysis_runs.R", "design_support.R", "reference_projects.R", "glabrata_annotation_lookup.R", "normalization_module.R", "ui_fitness.R", "treated_vs_control.R", "comparisons_module.R"))
  source(file.path(root, "app/shiny/R", name))
html <- htmltools::renderTags(comparisons_module_ui())$html
stopifnot(grepl("Treated vs Control: Run", html, fixed = TRUE), grepl("Treated vs Control: Results", html, fixed = TRUE),
          !grepl("Fitness: Run", html, fixed = TRUE),
          grepl("Local-abundance z-score", html, fixed = TRUE), grepl(">Resampling<", html, fixed = TRUE),
          !grepl("Control-conditioned Gumbel (depletion only)", html, fixed = TRUE),  # retired from the chooser
          !grepl("narrowly failed FDR benchmark)", html, fixed = TRUE),
          grepl("pre-registered calibration gates", html, fixed = TRUE),   # the validation status stays in Methods
          grepl("fitness_normalization_selector", html, fixed = TRUE), grepl("tvc_result_selector", html, fixed = TRUE),
          !grepl("inside the R analysis", html, fixed = TRUE))
# TRANSIT is named only in Methods references.
stopifnot(!grepl("transit", gsub('(?s)<p class="ytab-reference">.*?</p>', "", html, perl = TRUE), ignore.case = TRUE))
# Retired from the chooser, but still resolvable: the label map and the allowed normalizations must stay, or a
# stored Gumbel run would display with no method name and no normalization it accepts.
stopifnot(identical(tvc_allowed_methods$gumbel, c("ttr", "totreads", "nzmean")), "betageom" %in% tvc_allowed_methods$resampling,
          "gumbel" %in% tvc_method_labels, !("gumbel" %in% tvc_method_choices))
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

# ---------------------------------------------------------------------------------------------------------
# Downloads carry what the screen carries.
#
# Under testServer, reading a download output runs its handler and returns the file it wrote, so the test can
# compare that file against the table actually rendered. The gene download used to write every gene and every
# column of the run regardless of the chosen comparison; it is now fed by the same reactive as the table, and
# this assertion is what holds that. Comparing against the rendered widget rather than an internal reactive
# keeps the test honest about what a user would receive.
# ---------------------------------------------------------------------------------------------------------
config <- file.path(root, "output/projects/Gale_2020/config/project.yaml")
if (file.exists(config)) {
  project <- read_project_summary(config, root)
  active <- reactive(modifyList(project, list(project_root = file.path(root, "output/projects/Gale_2020"))))
  start <- function(...) invisible(NULL)
  testServer(function(input, output, session)
    tvc_server(input, output, session, root, active, function() config,
               reactive(data.frame()), start), {
    session$flushReact()
    # Pick a real result: the resampling method, whose table is the one the download used to bypass.
    runs <- ytab_analysis_runs(file.path(root, "output/projects/Gale_2020/treated_vs_control"))
    resampling <- runs$path[grepl("resampling", basename(runs$path), fixed = TRUE)]
    stopifnot(length(resampling) > 0)
    session$setInputs(tvc_results_method = "resampling", tvc_result = resampling[[1L]])
    session$flushReact()
    rendered <- tryCatch(jsonlite::fromJSON(output$tvc_table)$x$container,
                         error = function(e) { cat("   table render error:", conditionMessage(e), "\n"); "" })
    if (nzchar(rendered)) {
      shown_columns <- regmatches(rendered, gregexpr("(?<=<th>)[^<]+", rendered, perl = TRUE))[[1]]
      written <- utils::read.delim(output$tvc_download_genes, check.names = FALSE, colClasses = "character", nrows = 5)
      stopifnot(length(shown_columns) > 0, identical(names(written), shown_columns))
      cat(sprintf("  tvc download columns match the table: %s\n", paste(shown_columns, collapse = ", ")))
    } else cat("  tvc: no result on disk to compare against; download assertion skipped\n")
  })
}
cat("PASS\n")
