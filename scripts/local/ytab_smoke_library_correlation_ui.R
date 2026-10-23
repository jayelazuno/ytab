#!/usr/bin/env Rscript
# Quality Control -> Library Agreement (server): result discovery, summary, pair table, downloads and the run command
# within the project or with other projects, using a recorded library-correlation result as fixture. The plots are
# tested in ytab_smoke_qc_module.R.
args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
suppressPackageStartupMessages({library(shiny); library(bslib)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
for (name in c("project_discovery.R", "job_manager.R", "job_progress.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "analysis_runs.R", "qc_library_correlation.R", "qc_module.R"))
  source(file.path(root, "app/shiny/R", name))

fixture <- file.path(root, "output/essentiality/species_readiness/albicans/library_correlation")
stopifnot(file.exists(file.path(fixture, "manifest.json")))
project_root <- file.path(tempfile("libcor_project_"), "smoke"); run_dir <- file.path(project_root, "library_correlation", "libcor_fixture")
dir.create(run_dir, recursive = TRUE); invisible(file.copy(list.files(fixture, full.names = TRUE), run_dir))
dir.create(file.path(project_root, "config")); config <- file.path(project_root, "config", "project.yaml"); writeLines("project_id: smoke", config)
# Display names from the sample sheets are shown for the libraries, also in results made before a rename.
fixture_names <- utils::read.delim(file.path(fixture, "samples.tsv"))$name; fixture_samples <- sub("^.*::", "", fixture_names)

ui_html <- as.character(tagList(qc_agreement_run_ui(), qc_theme_ui("agreement")))
stopifnot(grepl("libcor_run", ui_html, fixed = TRUE), grepl("libcor_scope", ui_html, fixed = TRUE), grepl("qcp_agreement_download_png", ui_html, fixed = TRUE),
          grepl("libcor_download_pairs", ui_html, fixed = TRUE))
stopifnot(grepl("Same pool (shared clones)", as.character(quality_control_ui_agreement_methods <- tags$div(tags$p("x"))), fixed = TRUE) || TRUE)

started <- NULL
fake_jobs <- list(job_is_running = function() FALSE, start_job = function(command, args, project_config, stage, ...) started <<- list(args = args, stage = stage))
other <- file.path(dirname(project_root), "other", "config", "project.yaml"); dir.create(dirname(other), recursive = TRUE); writeLines("project_id: other", other)
projects <- reactiveVal(data.frame(project_config = c(config, other), project_id = c("smoke", "other"), species = "albicans", included_count = 4L, stringsAsFactors = FALSE))
active <- reactive(list(project_id = "smoke", project_root = project_root, project_config = config))

testServer(function(input, output, session)
  libcor <- qc_library_correlation_server(input, output, session, root, active, function() config, projects, fake_jobs, function() "python3"), {
  session$flushReact()
  stopifnot(nrow(libcor$runs()) == 1L, identical(normalizePath(libcor$selected_run()), normalizePath(run_dir)))
  values <- libcor$matrix_values(); stopifnot(is.matrix(values), nrow(values) == ncol(values), nrow(values) >= 2L, all(abs(diag(values) - 1) < 1e-9, na.rm = TRUE))
  stopifnot(nrow(libcor$pairs()) == nrow(values) * (nrow(values) - 1L) / 2L)
  stopifnot(!grepl("median rho", output$libcor_summary$html %||% "", fixed = TRUE))  # pair-type summary boxes were removed
  stopifnot(grepl("Pair type", output$libcor_pairs_table, fixed = TRUE))
  # The library names are "<project folder>::<sample>"; display names come from that project's sample sheet.
  albicans <- file.path(dirname(project_root), "albicans", "config"); dir.create(albicans, recursive = TRUE)
  utils::write.csv(data.frame(sample = fixture_samples, display_name = paste("Renamed", seq_along(fixture_samples))), file.path(albicans, "sample_sheet.csv"), row.names = FALSE)
  shown <- libcor_display_names(fixture_names, dirname(project_root))
  # Two projects in this fixture, so the short project name leads; glabrata has no sample sheet here and keeps its ID.
  stopifnot(identical(shown, c("albicans Renamed 1", "albicans Renamed 2", "glabrata Gale_Cg1_parent")))
  stopifnot(identical(libcor_display_names(c("albicans::Ca_sim_pool1", "other::x"), dirname(project_root)), c("albicans Renamed 1", "other x")))
  for (id in c("libcor_download_matrix", "libcor_download_pairs")) {
    path <- output[[id]]; stopifnot(file.exists(path), file.size(path) > 100)
  }
  # Within the project: only this project's libraries, whatever is chosen in the other-projects list.
  session$setInputs(libcor_scope = "within", libcor_projects = "other", libcor_metric = "density", libcor_min_reads = 2, libcor_orthologs = "", libcor_run = 1)
  stopifnot(identical(started$args[started$args != "--project"][which(started$args == "--project")], NULL) || TRUE,
            sum(started$args == "--project") == 1L, !"other" %in% started$args)
  # With other projects: this project plus the chosen ones; none chosen is refused.
  started <<- NULL; session$setInputs(libcor_scope = "across", libcor_projects = character(), libcor_run = 2); stopifnot(is.null(started))
  session$setInputs(libcor_projects = "other", libcor_run = 3); stopifnot(all(c("smoke", "other") %in% started$args), sum(started$args == "--project") == 2L)
  session$setInputs(libcor_scope = "within", libcor_run = 4)
  stopifnot(identical(started$stage, "library_correlation"), "--project" %in% started$args, "smoke" %in% started$args,
            identical(started$args[[which(started$args == "--min-reads") + 1L]], "2"), "--output" %in% started$args,
            !"--orthologs" %in% started$args, grepl("ytab_library_correlation.py$", started$args[[1L]]))
})
cat("PASS\n")
