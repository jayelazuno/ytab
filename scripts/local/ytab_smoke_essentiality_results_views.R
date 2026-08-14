#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(
  sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1L]]),
  winslash = "/"
)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(root, "app/shiny/R/project_discovery.R"))
source(file.path(root, "app/shiny/R/essentiality_targets.R"))
source(file.path(root, "app/shiny/R/essentiality_state.R"))
source(file.path(root, "app/shiny/R/essentiality_results.R"))

project <- read_project_summary(args[[1L]], root)
parents <- detect_essentiality_parent_samples(project)
results <- discover_classifier_results(project$project_root)
stopifnot(length(results) >= 1L)
selected_tag <- results[[1L]]$tag
selected_before <- choose_essentiality_result(results, selected = selected_tag)
source_state <- file.info(selected_before$path)[, c("size", "mtime"), drop = FALSE]

overview <- essentiality_results_view_contract("overview")
predictions <- essentiality_results_view_contract("predictions")
visualizations <- essentiality_results_view_contract("visualizations")
provenance_view <- essentiality_results_view_contract("provenance")
downloads_view <- essentiality_results_view_contract("downloads")
stopifnot(
  !overview$renders_predictions_table,
  "label_distribution" %in% overview$components,
  predictions$renders_predictions_table,
  all(c("search", "label_filter", "inclusion_filter") %in%
        predictions$components),
  !predictions$renders_visualizations,
  visualizations$renders_visualizations,
  !"predictions_table" %in% visualizations$components,
  provenance_view$renders_provenance,
  all(c("target_provenance", "technical_provenance") %in%
        provenance_view$components),
  downloads_view$renders_downloads
)

for (view in c("overview", "predictions", "visualizations", "provenance", "downloads")) {
  invisible(essentiality_results_view_contract(view))
  selected_after <- choose_essentiality_result(results, selected = selected_tag)
  stopifnot(identical(selected_after$tag, selected_before$tag))
}
provenance <- essentiality_result_provenance(
  project$project_root, selected_before, parents
)
stopifnot(identical(provenance$target_tag, selected_tag),
          is.logical(provenance$final))
downloads <- essentiality_download_paths(project$project_root, selected_before)
availability <- essentiality_download_availability(downloads)
stopifnot(isTRUE(availability[["predictions"]]),
          isTRUE(availability[["combined_feature"]]))
stopifnot(identical(
  file.info(selected_before$path)[, c("size", "mtime"), drop = FALSE],
  source_state
))
cat("PASS\n")
