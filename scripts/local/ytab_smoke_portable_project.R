#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[[1]], mustWork = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
source(file.path(root, "app/shiny/R/project_discovery.R"))
sandbox <- tempfile("ytab_portable_")
dir.create(sandbox)
sandbox <- normalizePath(sandbox, winslash = "/", mustWork = TRUE)
config_dir <- file.path(sandbox, "output/projects/demo/config")
dir.create(config_dir, recursive = TRUE)
reads <- file.path(sandbox, "external reads")
dir.create(reads)
yaml::write_yaml(list(project_id = "demo", repo_root = "../../../..",
  species = "glabrata", threads = 2, fastq_dir = reads,
  fastq_paths_relative_to = "fastq_dir",
  sample_sheet = "output/projects/demo/config/sample_sheet.csv"),
  file.path(config_dir, "project.yaml"))
write.csv(data.frame(sample = "parent", fastq_1 = "parent_1.fastq.gz",
  fastq_2 = "", condition = "parent", include = TRUE),
  file.path(config_dir, "sample_sheet.csv"), row.names = FALSE)
project <- read_project_summary(file.path(config_dir, "project.yaml"), sandbox)
stopifnot(identical(project$samples$fastq_1[[1]], file.path(reads, "parent_1.fastq.gz")),
          identical(project$samples$fastq_2[[1]], ""), project$parent_count == 1L,
          identical(project$export_root, file.path(sandbox, "output/exports/demo")))
cat("PASS: Shiny resolves FASTQ-relative filenames and standard downstream routes.\n")
