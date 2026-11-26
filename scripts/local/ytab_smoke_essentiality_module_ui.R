#!/usr/bin/env Rscript
# Essentiality module: runs `ytab essentiality` on Gale 2020 (reusing the recorded HMM fit; adding the current RF and the
# hybrid), then checks the default libraries, summary, three interactive plots, gene table, downloads and run commands.
args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
suppressPackageStartupMessages({library(shiny); library(bslib)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
for (name in c("job_manager.R", "process_helpers.R", "project_discovery.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "analysis_runs.R", "design_support.R", "reference_projects.R", "glabrata_annotation_lookup.R", "qc_plot_utils.R", "essentiality_state.R", "essentiality_module.R"))
  source(file.path(root, "app/shiny/R", name))
python <- unname(Sys.which("python"))
config <- file.path(root, "output/projects/Gale_2020/config/project.yaml")
hmm <- file.path(root, "output/essentiality/hermes_hmm/benchmark_v2/gale_blinded/joint_run1")
# The feature table of whichever combined target is present, rather than a fixed one: the legacy MidLC targets
# (T<number>) are gone, and a target is named after its normalization method now.
features <- {
  targets <- Sys.glob(file.path(root, "output/projects/Gale_2020/summary_combined/*/combined_feature_table.*.txt"))
  if (length(targets)) sort(targets)[[1]] else ""
}
labels <- file.path(root, "resources/species/glabrata/essentiality_labels/labels.tsv")
stopifnot(nzchar(python), all(file.exists(c(config, file.path(hmm, "fit.json"), features, labels))))

project <- read_project_summary(config, root)
# The stacked state posterior keeps one coloured area per state (a per-point hover label must not split the groups).
fixture_windows <- data.frame(chromosome = "c1", start = c(1, 101, 201), P_ES = c(.9, .1, 0), P_GD = c(.1, .2, 0), P_NE = c(0, .7, 1), P_GA = 0, insertions = c(0, 3, 9), mappable = 1)
built <- ggplot2::ggplot_build(ess_ggplot("track", data.frame(), fixture_windows, "c1")$bottom)
stopifnot(setequal(unique(built$data[[1]]$fill), unname(essentiality_state_colours)))
stopifnot(identical(sort(essentiality_default_libraries(essentiality_parent_inventory(project))), c("SRR11561506", "SRR11561507", "SRR11561508")))

project_root <- file.path(tempfile("ess_project_"), "Gale_2020"); run_dir <- file.path(project_root, "essentiality", "essentiality_fixture")
dir.create(dirname(run_dir), recursive = TRUE)
status <- system2(python, c(file.path(root, "scripts/local/ytab_essentiality.py"), "--project-config", config, "--hmm-run", hmm,
                            "--features", features, "--labels", labels, "--output", run_dir), stdout = FALSE, stderr = FALSE)
stopifnot(status == 0L, file.exists(file.path(run_dir, "models/genes.tsv")))

essentiality_combine_ui <- essentiality_combined_summary_ui <- essentiality_classifier_ui <- function() tags$div("stage placeholder")
ui_html <- as.character(essentiality_module_ui())
stopifnot(!grepl("Hybrid Inputs", ui_html, fixed = TRUE), grepl("hybrid-inputs", ui_html, fixed = TRUE), grepl("ess_run", ui_html, fixed = TRUE), grepl("0.947", ui_html, fixed = TRUE),
          grepl("ess_download_genes", ui_html, fixed = TRUE), grepl("HapA", ui_html, fixed = TRUE),
          !grepl("Current random forest", ui_html, fixed = TRUE), !grepl("current random forest", ui_html, fixed = TRUE))

started <- NULL
fake_jobs <- list(job_is_running = function() FALSE, start_job = function(command, args, project_config, stage, ...) started <<- list(args = args, stage = stage))
active <- reactive(modifyList(project, list(project_root = project_root)))
testServer(function(input, output, session)
  ess <- essentiality_module_server(input, output, session, root, active, function() config, fake_jobs, function() python), {
  session$flushReact()
  m <- ess$manifest(); stopifnot(nrow(ess$runs()) == 1L, isTRUE(m$hybrid), setequal(unlist(m$models), c("hmm", "hybrid")))
  stopifnot(m$hmm$gene_calls$ES == 1382L, nrow(ess$genes()) == 5234L, nrow(ess$windows()) > 100000L)
  # The results open on the HMM even when a hybrid exists: the comparison is a view the user asks for, not the
  # default reading. Each view shows its own cards.
  stopifnot(grepl("Essential (ES)", output$ess_summary$html, fixed = TRUE),
            !grepl("Consensus essential", output$ess_summary$html, fixed = TRUE))
  session$setInputs(ess_view = "compare")
  stopifnot(grepl("Consensus essential", output$ess_summary$html, fixed = TRUE),
            grepl("HMM only", output$ess_summary$html, fixed = TRUE))
  session$setInputs(ess_view = "hybrid")
  stopifnot(grepl("genes (hybrid)", output$ess_summary$html, fixed = TRUE),
            !grepl("Growth defect", output$ess_summary$html, fixed = TRUE))
  session$setInputs(ess_view = "hmm")
  session$setInputs(ess_plot = "track", ess_chromosome = "CP048230.1"); stopifnot(grepl("plotly", output$ess_plotly, fixed = TRUE))
  for (choice in c("scores", "calls")) { session$setInputs(ess_plot = choice); stopifnot(grepl("plotly", output$ess_plotly, fixed = TRUE)) }
  # every view renders every plot it offers
  for (v in c("hmm", "hybrid", "compare")) { session$setInputs(ess_view = v)
    for (choice in c("scores", "calls")) { session$setInputs(ess_plot = choice)
      stopifnot(grepl("plotly", output$ess_plotly, fixed = TRUE)) } }
  # The gene table and the plot captions follow the view too: an Agreement column, or a caption describing the
  # HMM-vs-hybrid scatter, is the comparison reappearing where it was not asked for.
  # testServer hands back the rendered widget as JSON; the column headers are in the table's container HTML.
  has_col <- function(name) grepl(paste0(">", name, "<"), jsonlite::fromJSON(output$ess_genes_table)$x$container, fixed = TRUE)
  caption <- function() gsub("<[^>]*>", "", paste(as.character(output$ess_plot_caption), collapse = " "))
  session$setInputs(ess_view = "hmm", ess_plot = "scores")
  stopifnot(has_col("HMM call"), !has_col("Agreement"), !has_col("Hybrid score"),
            grepl("HMM score (P_ES_run)", caption(), fixed = TRUE))
  session$setInputs(ess_view = "hybrid")
  stopifnot(has_col("Hybrid call"), has_col("Hybrid score"), !has_col("HMM call"),
            grepl("hybrid score", caption(), fixed = TRUE))
  # Both models are read through the same plots: same list, same forms, only the model changes.
  plots <- function() { h <- as.character(output$ess_plot_selector$html %||% output$ess_plot_selector)
    vapply(c("track", "scores", "calls", "target", "roc"), function(v) grepl(paste0('value="', v, '"'), h, fixed = TRUE), logical(1)) }
  session$setInputs(ess_view = "hmm");     hmm_plots <- plots()
  session$setInputs(ess_view = "compare"); cmp_plots <- plots()
  session$setInputs(ess_view = "hybrid");  hyb_plots <- plots()
  stopifnot(identical(hmm_plots, cmp_plots), identical(hmm_plots, hyb_plots), all(hmm_plots))
  # The gene table carries coordinates, so the hybrid draws its own chromosome track.
  g <- ess$genes()
  stopifnot(all(c("start", "end") %in% names(g)), any(!is.na(g$start)), all(g$end[!is.na(g$start)] > g$start[!is.na(g$start)]))
  # The chromosome chooser decides the detail: none is the genome at a glance, one is the full track, several
  # put their insertion landscapes side by side. Every view does all three.
  for (v in c("hmm", "hybrid", "compare")) {
    session$setInputs(ess_view = v, ess_plot = "track", ess_chromosome = character(0))
    stopifnot(grepl("plotly", output$ess_plotly, fixed = TRUE))
    session$setInputs(ess_chromosome = "CP048230.1")
    stopifnot(grepl("plotly", output$ess_plotly, fixed = TRUE))
    session$setInputs(ess_chromosome = c("CP048230.1", "CP048231.1"))
    stopifnot(grepl("plotly", output$ess_plotly, fixed = TRUE))
  }
  session$setInputs(ess_view = "hmm", ess_chromosome = character(0))
  # A download is a copy of what is on screen. The plot follows the view and the chromosome selection; the
  # gene table carries the columns being shown, not every column the pipeline wrote.
  # The track has no server-side plot download any more -- it is captured from the figure -- so what is
  # checked per view is that the figure itself renders in each of the three selection modes.
  for (v in c("hmm", "hybrid", "compare")) {
    for (selection in list(character(0), "CP048230.1", c("CP048230.1", "CP048231.1"))) {
      session$setInputs(ess_view = v, ess_plot = "track", ess_chromosome = selection)
      stopifnot(grepl("plotly", output$ess_plotly, fixed = TRUE))
    }
  }
  session$setInputs(ess_view = "hmm", ess_chromosome = character(0))
  hmm_tsv <- readLines(output$ess_download_genes, n = 1)
  session$setInputs(ess_view = "hybrid")
  hyb_tsv <- readLines(output$ess_download_genes, n = 1)
  session$setInputs(ess_view = "compare")
  cmp_tsv <- readLines(output$ess_download_genes, n = 1)
  stopifnot(grepl("HMM call", hmm_tsv, fixed = TRUE), !grepl("Hybrid score", hmm_tsv, fixed = TRUE),
            grepl("Hybrid call", hyb_tsv, fixed = TRUE), !grepl("HMM call", hyb_tsv, fixed = TRUE),
            grepl("Agreement", cmp_tsv, fixed = TRUE))
  session$setInputs(ess_view = "hmm", ess_plot = "track")
  session$setInputs(ess_view = "compare", ess_plot = "scores")
  stopifnot(has_col("HMM call"), has_col("Hybrid call"), has_col("Agreement"),
            grepl("both models", caption(), fixed = TRUE))
  session$setInputs(ess_view = "hmm", ess_plot = "calls")
  stopifnot(grepl("per HMM call", caption(), fixed = TRUE), !grepl("hybrid", caption(), fixed = TRUE))
  # The chromosome track is saved by capturing the plotly figure on screen, so there is no second, server-side
  # drawing of it to drift: that is the whole point of the capture. What the test can check is that the track
  # is wired that way -- capture buttons present, no server-side rebuild left -- and that the plots still
  # written server-side keep working.
  ui <- as.character(htmltools::renderTags(essentiality_model_results_ui())$html)
  stopifnot(grepl("ess_capture_png", ui, fixed = TRUE), grepl("ess_capture_svg", ui, fixed = TRUE),
            !grepl("ess_download_png", ui, fixed = TRUE), !grepl("ess_download_pdf", ui, fixed = TRUE),
            !exists("ess_track_ggplot"))
  session$setInputs(ess_view = "hmm", ess_plot = "track", ess_chromosome = character(0))
  stopifnot(grepl("plotly", output$ess_plotly, fixed = TRUE))      # the figure the capture photographs
  for (choice in c("scores", "calls", "target")) {
    session$setInputs(ess_plot = choice)
    stopifnot(inherits(ess$plot_for_download(), "ggplot"))          # still rebuilt server-side, still works
  }
  session$setInputs(ess_plot = "track")
  session$setInputs(ess_view = "hmm", ess_plot = "track")
  # The hybrid's score is in the table under the views that show the hybrid, not under the HMM's own view.
  stopifnot(!grepl("Hybrid score", output$ess_genes_table, fixed = TRUE))
  session$setInputs(ess_view = "compare")
  stopifnot(grepl("Hybrid score", output$ess_genes_table, fixed = TRUE))
  session$setInputs(ess_view = "hmm")
  session$setInputs(ess_plot = "scores")
  for (id in c("ess_download_png", "ess_download_pdf", "ess_download_genes")) { path <- output[[id]]; stopifnot(file.size(path) > 100) }
  session$setInputs(ess_plot = "track")
  session$setInputs(ess_reuse_fit = FALSE, ess_libraries = c("SRR11561508", "SRR11561507"), ess_window = 100, ess_min_reads = 1, ess_model = "hmm", ess_run = 1)
  stopifnot(identical(started$stage, "essentiality_models"), all(c("SRR11561508", "SRR11561507") %in% started$args), !"--labels" %in% started$args)
  # The hybrid takes a normalization, not a classifier output: its gene features are built by the run (BGC default).
  session$setInputs(ess_reuse_fit = TRUE, ess_hmm_reuse = run_dir, ess_model = "hybrid", ess_labels = labels, ess_normalization = "ttr", ess_run = 2)
  stopifnot("--hmm-run" %in% started$args, "--labels" %in% started$args, !"--rf" %in% started$args, !"--libraries" %in% started$args,
            identical(started$args[[match("--normalization", started$args) + 1L]], "ttr"))
})
cat("PASS\n")
