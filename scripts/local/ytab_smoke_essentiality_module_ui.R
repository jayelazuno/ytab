#!/usr/bin/env Rscript
# Essentiality module: runs `ytab essentiality` on Gale 2020 (reusing the recorded HMM fit; adding the current RF and the
# hybrid), then checks the default libraries, summary, three interactive plots, gene table, downloads and run commands.
args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
suppressPackageStartupMessages({library(shiny); library(bslib)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
for (name in c("job_manager.R", "project_discovery.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "analysis_runs.R", "glabrata_annotation_lookup.R", "qc_plot_utils.R", "essentiality_state.R", "essentiality_module.R"))
  source(file.path(root, "app/shiny/R", name))
python <- unname(Sys.which("python"))
config <- file.path(root, "output/projects/Gale_2020/config/project.yaml")
hmm <- file.path(root, "output/essentiality/hermes_hmm/benchmark_v2/gale_blinded/joint_run1")
features <- file.path(root, "output/projects/Gale_2020/summary_combined/T013p239/combined_feature_table.T013p239.txt")
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
  stopifnot(grepl("Consensus essential", output$ess_summary$html, fixed = TRUE))
  session$setInputs(ess_plot = "track", ess_chromosome = "CP048230.1"); stopifnot(grepl("plotly", output$ess_plotly, fixed = TRUE))
  for (choice in c("scores", "calls")) { session$setInputs(ess_plot = choice); stopifnot(grepl("plotly", output$ess_plotly, fixed = TRUE)) }
  stopifnot(grepl("Hybrid score", output$ess_genes_table, fixed = TRUE))
  session$setInputs(ess_plot = "track")
  for (id in c("ess_download_png", "ess_download_pdf", "ess_download_genes")) { path <- output[[id]]; stopifnot(file.size(path) > 100) }
  session$setInputs(ess_hmm_mode = "fit", ess_libraries = c("SRR11561508", "SRR11561507"), ess_window = 100, ess_min_reads = 1, ess_use_labels = FALSE, ess_run = 1)
  stopifnot(identical(started$stage, "essentiality_models"), all(c("SRR11561508", "SRR11561507") %in% started$args), !"--labels" %in% started$args)
  # The hybrid takes a normalization, not a classifier output: its gene features are built by the run (BGC default).
  session$setInputs(ess_hmm_mode = "reuse", ess_hmm_reuse = run_dir, ess_use_labels = TRUE, ess_labels = labels, ess_normalization = "ttr", ess_run = 2)
  stopifnot("--hmm-run" %in% started$args, "--labels" %in% started$args, !"--rf" %in% started$args, !"--libraries" %in% started$args,
            identical(started$args[[match("--normalization", started$args) + 1L]], "ttr"))
})
cat("PASS\n")
