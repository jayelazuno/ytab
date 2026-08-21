#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1]]), winslash = "/")
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(root, "app/shiny/R/project_discovery.R"))
source(file.path(root, "app/shiny/R/glabrata_annotation_lookup.R"))
source(file.path(root, "app/shiny/R/essentiality_targets.R"))
source(file.path(root, "app/shiny/R/essentiality_state.R"))
source(file.path(root, "app/shiny/R/essentiality_results.R"))
project <- read_project_summary(args[[1]], root)
parents <- detect_essentiality_parent_samples(project)
results <- discover_classifier_results(project$project_root)
stopifnot(length(results) >= 1L)
result <- choose_essentiality_result(results)
data <- read_essentiality_predictions(result$path)
stopifnot(nrow(data) == result$rows, nrow(data) > 0L)
summary <- essentiality_result_summary(result)
stopifnot(sum(summary$counts) == summary$total)
columns <- essentiality_prediction_columns(data)
label <- if (nzchar(columns$label)) as.character(data[[columns$label]][[1]]) else "All"
filtered <- filter_essentiality_results(data, label = label)
stopifnot(nrow(filtered) > 0L, nrow(filtered) <= nrow(data))
optional <- data.frame(feature_id = c("a", "b"), classifier_label = c("1", "2"),
                       stringsAsFactors = FALSE)
stopifnot(nrow(essentiality_visible_results(optional)) == 2L)
provenance <- essentiality_result_provenance(project$project_root, result, parents)
stopifnot(provenance$target_tag == result$tag)
stopifnot(is.logical(provenance$final), essentiality_smoke_project(project))
visible <- essentiality_visible_results(data)
stopifnot(!any(vapply(visible, function(column)
  any(grepl(project$project_root, as.character(column), fixed = TRUE)), logical(1))))
visible_mapped <- essentiality_visible_results(data, root)
stopifnot(nrow(visible_mapped) == nrow(visible),
          all(c("CAGL ID", "Gene name", "Cg-to-Sc relationship") %in% names(visible_mapped)),
          !"C. glabrata DESeq gene name" %in% names(visible_mapped))
downloads <- essentiality_download_paths(project$project_root, result)
stopifnot(file.exists(downloads$predictions), file.exists(downloads$combined_feature),
          file.exists(downloads$target_evaluation),
          file.exists(downloads$normalization_recommendation),
          file.exists(downloads$classifier_manifest))
server_text <- paste(readLines(file.path(root, "app/shiny/R/essentiality_server.R"),
                               warn = FALSE), collapse = "\n")
stopifnot(grepl("open = if \\(failed\\) NA else NULL", server_text))
cat("PASS\n")
