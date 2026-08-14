#!/usr/bin/env Rscript

usage <- paste(
  "Usage: Rscript scripts/ytab_fitness_condition_control_scatter.R",
  "--result-dir output/projects/PROJECT/treated_vs_parent/ANALYSIS",
  "[--output plots/condition_vs_control_loglog_scatter.png]"
)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) {
      cat(usage, "\n")
      quit(status = 0)
    }
    if (!startsWith(key, "--") || i == length(args)) stop("Invalid arguments.\n", usage)
    out[[gsub("-", "_", substring(key, 3L))]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo_root <- normalizePath(file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)), ".."),
                           winslash = "/", mustWork = TRUE)

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$result_dir) || !nzchar(args$result_dir)) stop("Missing --result-dir.\n", usage)

`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(repo_root, "app", "shiny", "R", "qc_plot_utils.R"))
source(file.path(repo_root, "app", "shiny", "R", "fitness_condition_control_plot.R"))

result_dir <- normalizePath(args$result_dir, winslash = "/", mustWork = TRUE)
result <- list(
  output_dir = result_dir,
  table = file.path(result_dir, "treated_vs_parent_results.csv")
)

output <- args$output %||% fitness_condition_control_plot_file(result)
if (!grepl("^(/|[A-Za-z]:)", output)) output <- file.path(result_dir, output)

path <- save_fitness_condition_control_scatter(result, output)
if (!nzchar(path) || !file.exists(path)) stop("Could not generate condition-versus-control scatter plot.")
cat(normalizePath(path, winslash = "/", mustWork = TRUE), "\n")
