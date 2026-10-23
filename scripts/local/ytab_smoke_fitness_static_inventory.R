#!/usr/bin/env Rscript
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]]))
root <- normalizePath(file.path(dirname(self), "../.."))
`%||%` <- function(x,y) if(is.null(x)||!length(x)) y else x
diagnostic_encode_relative_path <- function(path) utils::URLencode(path,reserved=TRUE)
source(file.path(root,"app/shiny/R/fitness_generated_plots.R"))
directory <- tempfile("ytab-static-plots-")
plots <- file.path(directory,"treated_vs_parent","arbitrary-analysis","plots")
dir.create(plots,recursive=TRUE)
names <- c("MA__14_labeled_top_hits.png","MA_unknown_background_labeled_top_hits.png",
  "MA_treated_vs_parent_combined.png","ma_future_pool.png","mapping_quality.png","other_diagnostic.png")
stopifnot(all(file.create(file.path(plots,names))))
inventory <- fitness_generated_plot_inventory(list(project_root=directory))
stopifnot(setequal(inventory$filename,c("mapping_quality.png","other_diagnostic.png")),
  all(file.exists(file.path(plots,names))))
cat("PASS: static MA images hidden across arbitrary backgrounds; files and unrelated diagnostics preserved\n")
