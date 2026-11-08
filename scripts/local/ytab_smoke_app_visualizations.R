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

# The QC plot menus moved from ui_qc.R into qc_module.R (qc_plot_menus) when Quality Control was rebuilt around the
# flagged QC table; these are the same three choices under their current names.
must_contain("app/shiny/R/qc_module.R", "Saturation curves", "Saturation (library complexity) plot choice")
must_contain("app/shiny/R/qc_module.R", "Genes hit, libraries of a role combined", "Combined features hit plot choice")
must_contain("app/shiny/R/qc_module.R", "Correlation between libraries", "Library concordance plot choice")
# Treated vs Control (all three methods), Multi-condition and GI share one display-options panel each.
must_contain("app/shiny/R/treated_vs_control.R", "ytab_plot_customization_controls\\(\"tvc\"", "Treated vs Control display controls")
must_contain("app/shiny/R/treated_vs_control.R", "ytab_plot_customization_values\\(input, \"tvc\"\\)", "Treated vs Control display binding")
must_contain("app/shiny/R/comparisons_module.R", "ytab_plot_customization_controls\\(\"mc\"", "Multi-condition display controls")
must_contain("app/shiny/R/comparisons_module.R", "ytab_plot_customization_controls\\(\"gi\"", "Genetic interactions display controls")
must_contain("app/shiny/R/treated_vs_control.R", "s[$]point_size", "Treated vs Control point size binding")
must_contain("app/shiny/R/ui_gene_domain_explorer.R", "gene_domain_display_height", "Gene Explorer image height control")
must_contain("app/shiny/R/ui_gene_domain_explorer.R", "imageOutput\\(\"gene_domain_figure\"", "Gene Explorer image element lives in the static UI, not a renderUI")
must_contain("app/shiny/R/gene_domain_explorer_server.R", "gene_domain_has_figure", "Gene Explorer figure visibility flag")
must_contain("app/shiny/www/ytab_release_ui.css", "\\.ytab-static-image img", "Static image sizing CSS")
must_contain("app/shiny/www/ytab_release_ui.css", "object-fit: contain", "Static image aspect policy")

cat("PASS\n")
