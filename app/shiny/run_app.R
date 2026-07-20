#!/usr/bin/env Rscript

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop(
    "The R package 'shiny' is required. Activate the ytab-local environment or install Shiny first.",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
launcher <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else file.path("app", "shiny", "run_app.R")
app_dir <- normalizePath(dirname(launcher), mustWork = TRUE)
shiny::runApp(appDir = app_dir, launch.browser = interactive())
