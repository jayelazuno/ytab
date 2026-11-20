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
source(file.path(root, "app/shiny/R/essentiality_state.R"))
source(file.path(root, "app/shiny/R/ui_helpers.R"))
for (name in c("ui_components.R", "analysis_runs.R", "design_support.R", "reference_projects.R", "essentiality_module.R", "normalization_module.R")) source(file.path(root, "app/shiny/R", name))

project <- read_project_summary(args[[1L]], root)
parents <- essentiality_parent_inventory(project)
# Essentiality: Run / Results (hybrid inputs are built by Run); MidLC is a Site Normalization method (no subtab of its own).
workspace_html <- paste(htmltools::renderTags(essentiality_module_ui())$html, htmltools::renderTags(normalization_module_ui())$html)
stopifnot(!grepl("Hybrid Inputs", workspace_html, fixed = TRUE), grepl("Run essentiality", workspace_html, fixed = TRUE),
          !grepl('data-value="midlc"', workspace_html, fixed = TRUE), grepl('value="midlc"', workspace_html, fixed = TRUE), !grepl("MidLC (Classifier)", workspace_html, fixed = TRUE),
          !grepl("essentiality_target_selection", workspace_html, fixed = TRUE))
stopifnot(!grepl("Evaluate Targets", workspace_html, fixed = TRUE))
stopifnot(!grepl("ytab-workflow-tracker", workspace_html, fixed = TRUE))
stopifnot(!grepl("Current working target", workspace_html, fixed = TRUE))
stopifnot(!grepl("summary_normalized_run_panel", workspace_html, fixed = TRUE))  # retired with the MidLC subtab
stopifnot(!grepl("The normalization target is a MidLC depth target", workspace_html, fixed = TRUE))
stopifnot(!grepl("Continue essentiality workflow", workspace_html, fixed = TRUE))
stopifnot(!grepl("essentiality_results_view", workspace_html, fixed = TRUE))
stopifnot(!grepl("Search classifier predictions", workspace_html, fixed = TRUE))
stopifnot(!grepl("Run workflow steps", workspace_html, fixed = TRUE))  # legacy MidLC subtab retired

cat("PASS\n")
