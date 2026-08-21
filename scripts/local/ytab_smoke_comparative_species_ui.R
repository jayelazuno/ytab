#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(
  sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1L]]),
  winslash = "/"
)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
library(shiny)
library(bslib)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(root, "app/shiny/R/project_discovery.R"))
source(file.path(root, "app/shiny/R/job_progress.R"))
source(file.path(root, "app/shiny/R/sample_selector.R"))
source(file.path(root, "app/shiny/R/ui_helpers.R"))
source(file.path(root, "app/shiny/R/ui_components.R"))
source(file.path(root, "app/shiny/R/plot_customization_helpers.R"))
source(file.path(root, "app/shiny/R/plot_display_helpers.R"))
source(file.path(root, "app/shiny/R/table_display_helpers.R"))
source(file.path(root, "app/shiny/R/essentiality_targets.R"))
source(file.path(root, "app/shiny/R/essentiality_state.R"))
source(file.path(root, "app/shiny/R/essentiality_results.R"))
source(file.path(root, "app/shiny/R/fitness_design_state.R"))
source(file.path(root, "app/shiny/R/fitness_result_state.R"))
source(file.path(root, "app/shiny/R/comparative_resources.R"))
source(file.path(root, "app/shiny/R/comparative_project_data.R"))
source(file.path(root, "app/shiny/R/ui_comparative.R"))
source(file.path(root, "app/shiny/R/qc_library_diagnostics_plots.R"))
source(file.path(root, "app/shiny/R/ui_qc.R"))
source(file.path(root, "app/shiny/R/ui_landing.R"))
source(file.path(root, "app/shiny/R/ui_preprocessing.R"))
source(file.path(root, "app/shiny/R/qc_result_state.R"))
source(file.path(root, "app/shiny/R/ui_fitness.R"))
source(file.path(root, "app/shiny/R/ui_essentiality.R"))
source(file.path(root, "app/shiny/R/ui_gene_domain_explorer.R"))
source(file.path(root, "app/shiny/R/ui_workspace.R"))

project <- read_project_summary(args[[1L]], root)
stopifnot(nzchar(project$project_id))
comparative_html <- htmltools::renderTags(comparative_ui())$html
checkbox_html <- htmltools::renderTags(comparative_species_checkbox_ui())$html
workspace_html <- htmltools::renderTags(workspace_ui(2L))$html
for (label in c("Single Species View", "Comparative View",
                "Gene Group Analysis"))
  stopifnot(grepl(label, comparative_html, fixed = TRUE))
stopifnot(grepl("Comparative Species", workspace_html, fixed = TRUE),
          grepl("Fitness Screen", workspace_html, fixed = TRUE),
          grepl("Genome Browser", workspace_html, fixed = TRUE),
          grepl("Search for one gene to view available Tn-seq results for one species.",
                comparative_html, fixed = TRUE),
          !grepl("placeholder", checkbox_html, ignore.case = TRUE),
          !grepl(paste0("codex/", paste0("RNA", "cross")), comparative_html, fixed = TRUE),
          !grepl(file.path("docs", "codex", paste0("RNA", "cross")),
                 comparative_html, fixed = TRUE))
availability <- summarize_species_project_availability(root)
availability_display <- comparative_availability_display(availability)
stopifnot("No YTAB project outputs are available for this species yet." %in%
            availability$status ||
            any(availability$project_count == 0L),
          !any(c("species", "label", "enabled", "placeholder") %in%
                 names(availability_display)),
          all(c("Species", "YTAB project outputs", "Projects",
                "Essentiality results", "Fitness results",
                "Insertion summaries") %in% names(availability_display)))
cat("PASS\n")
