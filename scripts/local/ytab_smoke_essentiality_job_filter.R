#!/usr/bin/env Rscript
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(
  sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1L]]),
  winslash = "/"
)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
library(shiny)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(root, "app/shiny/R/process_helpers.R"))
source(file.path(root, "app/shiny/R/job_manager.R"))
source(file.path(root, "app/shiny/R/job_progress.R"))

mapping <- c(
  normalize = "sample_normalization",
  evaluate = "summary_normalized",
  combine = "combined_hits",
  combined_summary = "summary_combined",
  classifier = "classifier"
)
rejected <- c(
  "treated_vs_parent", "fitness", "fitness_analysis", "mapfastq",
  "mapping", "create_hit_file", "summary_table", "summary",
  "library_diagnostics"
)
for (expected in unname(mapping)) {
  stopifnot(essentiality_job_matches(list(stage = expected), expected))
  for (stage in setdiff(unname(mapping), expected))
    stopifnot(!essentiality_job_matches(list(stage = stage), expected))
  for (stage in rejected)
    stopifnot(!essentiality_job_matches(list(stage = stage), expected))
}
stopifnot(!essentiality_job_matches(NULL, "classifier"),
          !essentiality_job_matches("classifier", "classifier"))

project_root <- tempfile("ytab-essentiality-job-filter-")
jobs_root <- file.path(project_root, "manifests", "jobs")
dir.create(jobs_root, recursive = TRUE)
on.exit(unlink(project_root, recursive = TRUE, force = TRUE), add = TRUE)
index <- 0L
for (stage in c(unname(mapping), rejected)) {
  index <- index + 1L
  atomic_json_write(
    file.path(jobs_root, sprintf("%03d.progress.json", index)),
    list(
      stage = stage, status = "success", target_tag = "T6400",
      execution_mode = "run", job_elapsed_seconds = index,
      updated_at = format(Sys.time() - index, tz = "UTC", usetz = TRUE)
    )
  )
}
for (expected in unname(mapping)) {
  history <- essentiality_stage_job_history(project_root, expected)
  stopifnot(nrow(history) == 1L,
            identical(history$Stage[[1L]], progress_stage_label(expected)),
            identical(history$Target[[1L]], "T6400"),
            !any(grepl("Fitness|Mapping|treated", history$Stage, ignore.case = TRUE)))
}
cat("PASS\n")
