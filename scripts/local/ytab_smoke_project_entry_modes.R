#!/usr/bin/env Rscript

`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
library(shiny)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]] %||% "scripts/local/ytab_smoke_project_entry_modes.R")
root <- normalizePath(file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)), "../.."), winslash = "/", mustWork = FALSE)
app_dir <- file.path(root, "app/shiny")
source(file.path(app_dir, "R", "ui_helpers.R"), local = TRUE)
source(file.path(app_dir, "R", "project_discovery.R"), local = TRUE)
source(file.path(app_dir, "R", "preprocessing_status.R"), local = TRUE)
source(file.path(app_dir, "R", "ui_landing.R"), local = TRUE)

html <- as.character(landing_ui(c("glabrata"), 2L, FALSE))
required <- c(
  "Start new project from FASTQ",
  "Run mapping and hit-file creation locally from FASTQ files.",
  "Import existing hit files",
  "Use CreateHitFile outputs that were already generated elsewhere, such as on an HPC.",
  "Open an existing YTAB project",
  "Continue analysis or view results for an existing YTAB project."
)
missing <- required[!vapply(required, grepl, logical(1), x = html, fixed = TRUE)]
if (length(missing)) stop("Landing entry mode text missing: ", paste(missing, collapse = "; "), call. = FALSE)

config <- file.path(root, "output/projects/H2O2_screen_v1/config/project.yaml")
if (file.exists(config)) {
  info <- read_project_summary(config, root)
  if (!identical(info$start_stage, "create_hit_file")) stop("H2O2 start_stage was not create_hit_file.", call. = FALSE)
  status <- build_sample_pipeline_status(info)
  if (!all(status$mapping == "imported_or_not_required")) stop("Imported project still requires Mapping.", call. = FALSE)
  if (!all(status$bam_index == "not_required")) stop("Imported project still requires BAM index files.", call. = FALSE)
  if (!all(status$hit_file == "imported_success")) stop("Imported hit-file status not preserved.", call. = FALSE)
}

cat("Project entry mode smoke passed.\n")
