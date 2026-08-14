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
source(file.path(root, "app/shiny/R/essentiality_targets.R"))
source(file.path(root, "app/shiny/R/essentiality_state.R"))
source(file.path(root, "app/shiny/R/essentiality_results.R"))
source(file.path(root, "app/shiny/R/ui_helpers.R"))
source(file.path(root, "app/shiny/R/ui_essentiality.R"))

project <- read_project_summary(args[[1L]], root)
parents <- detect_essentiality_parent_samples(project)
recommendation <- essentiality_recommendation(project$project_root)
states <- essentiality_stage_state(
  project, parents, recommendation$target_tag, NULL, root
)
stopifnot("evaluate" %in% names(states))
stopifnot(identical(
  format_target_label(6400, "T6400", "user"),
  "Normalization target: 6400"
))
stopifnot(identical(
  format_target_label(6400, "T6400", "technical"),
  "Target value: 6400; technical tag: T6400"
))
stopifnot(identical(format_target_label(6400, "T6400", "badge"), "T6400"))

completed_panel <- essentiality_run_panel_ui(
  TRUE, "Rerun or reconfigure", "Configure stage",
  tagList(tags$details(tags$summary("Advanced settings"), tags$p("Advanced")),
          actionButton("hidden_run", "Run stage", class = "btn-primary"))
)
completed_html <- htmltools::renderTags(completed_panel)$html
stopifnot(grepl("^<details", completed_html),
          grepl("ytab-rerun-disclosure", completed_html, fixed = TRUE),
          !grepl("<details[^>]* open", completed_html),
          grepl("Advanced settings", completed_html, fixed = TRUE))

workspace_html <- htmltools::renderTags(essentiality_ui())$html
stopifnot(grepl("Normalize &amp; Choose Target", workspace_html, fixed = TRUE))
stopifnot(!grepl("Evaluate Targets", workspace_html, fixed = TRUE))
stopifnot(!grepl("ytab-workflow-tracker", workspace_html, fixed = TRUE))
stopifnot(!grepl("Current working target", workspace_html, fixed = TRUE))
stopifnot(grepl("summary_normalized_run_panel", workspace_html, fixed = TRUE))
stopifnot(grepl("The normalization target is a MidLC depth target", workspace_html, fixed = TRUE))
stopifnot(!grepl("Continue essentiality workflow", workspace_html, fixed = TRUE))
stopifnot(!grepl("essentiality_results_view", workspace_html, fixed = TRUE))
stopifnot(grepl("Search classifier predictions", workspace_html, fixed = TRUE))
stopifnot(grepl("Run workflow steps", workspace_html, fixed = TRUE))
stopifnot(grepl("Target provenance", workspace_html, fixed = TRUE))

overview <- essentiality_results_view_contract()
predictions <- essentiality_results_view_contract("predictions")
stopifnot(identical(overview$view, "overview"),
          !overview$renders_predictions_table,
          predictions$renders_predictions_table,
          !predictions$renders_visualizations)

server_text <- paste(
  readLines(file.path(root, "app/shiny/R/essentiality_server.R"), warn = FALSE),
  collapse = "\n"
)
stopifnot(!grepl("Parent-library normalization complete", server_text, fixed = TRUE),
          grepl("open = if (failed) NA else NULL", server_text, fixed = TRUE),
          grepl('expected_stage = "sample_normalization"', server_text, fixed = TRUE),
          grepl("Run feature-level target evaluation", server_text, fixed = TRUE),
          !grepl("Feature-level evaluation has not been run yet.", server_text, fixed = TRUE),
          grepl("selected_target_notice", server_text, fixed = TRUE),
          !grepl('expected_stage = "fitness', server_text, fixed = TRUE))
cat("PASS\n")
