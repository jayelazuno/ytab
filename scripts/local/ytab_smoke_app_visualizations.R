#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
repo_root <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)

fail <- function(message) stop(message, call. = FALSE)
read_repo <- function(path) paste(readLines(file.path(repo_root, path), warn = FALSE), collapse = "\n")
must_contain <- function(path, pattern, label = pattern) {
  if (!grepl(pattern, read_repo(path), perl = TRUE)) fail(paste("Missing", label, "in", path))
}

invisible(parse(file = file.path(repo_root, "app/shiny/app.R")))

must_contain("app/shiny/R/ui_qc.R", "Library complexity", "Summary QC plot choices")
must_contain("app/shiny/R/ui_qc.R", "Combined features hit", "Combined features hit plot choice")
must_contain("app/shiny/R/ui_qc.R", "Library concordance", "Library concordance plot choice")
must_contain("app/shiny/R/ui_fitness.R", "ytab_plot_customization_controls\\(\"fitness\"", "Fitness display controls")
must_contain("app/shiny/app.R", "fitness_point_opacity", "Fitness point opacity binding")
must_contain("app/shiny/app.R", "fitness_point_size", "Fitness point size binding")
must_contain("app/shiny/R/ui_gene_domain_explorer.R", "gene_domain_display_height", "Gene Explorer image height control")
must_contain("app/shiny/R/gene_domain_explorer_server.R", "gene_domain_figure_ui", "Gene Explorer responsive image wrapper")
must_contain("app/shiny/www/ytab_release_ui.css", "\\.ytab-static-image img", "Static image sizing CSS")
must_contain("app/shiny/www/ytab_release_ui.css", "object-fit: contain", "Static image aspect policy")

cat("PASS\n")
