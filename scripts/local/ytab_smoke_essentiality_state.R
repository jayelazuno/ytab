#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1]]), winslash = "/")
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(root, "app/shiny/R/project_discovery.R"))
source(file.path(root, "app/shiny/R/essentiality_targets.R"))
source(file.path(root, "app/shiny/R/essentiality_state.R"))

project <- read_project_summary(args[[1]], root)
parents <- detect_essentiality_parent_samples(project)
stopifnot(length(parents) == 4L)
stopifnot(identical(parents, c("yH298-parent-pool1", "yH298-parent-pool2",
                              "yH299-parent-pool3", "yH299-parent-pool4")))
stopifnot(!any(grepl("treated", parents, ignore.case = TRUE)))
persisted_parents <- parents
persisted_target <- essentiality_recommendation(project$project_root)$tag
stopifnot(identical(persisted_parents, parents), nzchar(persisted_target))
available <- discover_essentiality_targets(project$project_root)
stopifnot(nrow(available) > 0L, persisted_target %in% available$Tag)
states <- essentiality_stage_state(project, parents, persisted_target, NULL, root)
stopifnot(states[["normalize"]] %in% c("complete", "cached"))
stopifnot(states[["evaluate"]] %in% c("complete", "cached"))
stopifnot(states[["combine"]] %in% c("complete", "cached"))
stopifnot(states[["combined_summary"]] %in% c("complete", "cached"))
missing <- essentiality_stage_state(project, parents, "T001", NULL, root)
stopifnot(missing[["combine"]] == "blocked",
          missing[["combined_summary"]] == "blocked",
          missing[["classifier"]] == "blocked")
files <- list.files(file.path(root, "app/shiny"), pattern = "\\.[Rr]$", recursive = TRUE,
                    full.names = TRUE)
text <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
stopifnot(!grepl("Existing project outputs are available under the project output directory.",
                 text, fixed = TRUE))
final <- essentiality_final_target(project$project_root)
stopifnot(!nzchar(final))
cat("PASS\n")
