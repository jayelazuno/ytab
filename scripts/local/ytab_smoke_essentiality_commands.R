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
source(file.path(root, "app/shiny/R/essentiality_commands.R"))
project <- read_project_summary(args[[1]], root)
parents <- detect_essentiality_parent_samples(project)
config <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
plans <- list(
  build_essentiality_command(root, config, "sample_normalization", "59.7", parents, "preview"),
  build_essentiality_command(root, config, "summary_normalized", "T059p7", parents, "run"),
  build_essentiality_command(root, config, "combined_hits", "T059p7", parents, "preview"),
  build_essentiality_command(root, config, "summary_combined", "T059p7", parents, "run", force = TRUE),
  build_essentiality_command(root, config, "classifier", "T059p7", parents, "run", seed = 0)
)
expected <- c("ytab_run_sample_normalization.py", "ytab_run_summary_normalized.py",
              "ytab_run_combine_hits.py", "ytab_run_summary_combined.py",
              "ytab_run_classifier.py")
for (index in seq_along(plans)) {
  plan <- plans[[index]]
  stopifnot(plan$project_config == config)
  stopifnot(basename(plan$full[[1]]) == expected[[index]])
  stopifnot(config %in% plan$full)
  stopifnot(all(parents %in% strsplit(plan$full[match("--samples", plan$full) + 1L],
                                      ",", fixed = TRUE)[[1]]) || index %in% c(4L, 5L))
  stopifnot(!any(grepl("treated", plan$parents, ignore.case = TRUE)))
}
stopifnot("--dry-run" %in% plans[[1]]$full, "--dry-run" %in% plans[[3]]$full)
stopifnot(!("--dry-run" %in% plans[[2]]$full), !("--dry-run" %in% plans[[4]]$full),
          !("--dry-run" %in% plans[[5]]$full))
stopifnot("--force" %in% plans[[4]]$full,
          !("--force" %in% plans[[1]]$full),
          !("--force" %in% plans[[5]]$full))
stopifnot("59.7" %in% plans[[1]]$full, "T059p7" %in% plans[[2]]$full,
          "T059p7" %in% plans[[5]]$full)
stopifnot(any(grepl("Repositories /ytab", plans[[1]]$full, fixed = TRUE)))
cat("PASS\n")
