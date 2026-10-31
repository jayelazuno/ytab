#!/usr/bin/env Rscript
# Normalization module: runs `ytab normalize-project` on the Zn screen (BGC and TTR), then checks result discovery,
# one native before/after panel per ticked method (scaling statistic, geometric QQ + ranked counts), the replicate-agreement
# view, the per-sample statistics table, downloads, the run command and legacy navigation. Unticked methods disappear;
# ticked methods missing from a run are reported; runs made before the per-method summaries get them on first view.
args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
suppressPackageStartupMessages({library(shiny); library(bslib)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
for (name in c("job_manager.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "analysis_runs.R","reference_projects.R", "navigation.R", "normalization_module.R"))
  source(file.path(root, "app/shiny/R", name))
python <- unname(Sys.which("python"))
config <- file.path(root, "output/projects/Zn_toxicity_screen/config/project.yaml")
stopifnot(nzchar(python), file.exists(config))

project_root <- file.path(tempfile("norm_project_"), "zn"); dir.create(file.path(project_root, "config"), recursive = TRUE)
invisible(file.copy(file.path(dirname(config), "comparison_design.csv"), file.path(project_root, "config")))
run_dir <- file.path(project_root, "normalization", "norm_fixture")
status <- system2(python, c(file.path(root, "scripts/local/ytab_normalize_project.py"), "--project-config", config,
                            "--method", "betageom", "--method", "ttr", "--output", run_dir), stdout = FALSE, stderr = FALSE)
stopifnot(status == 0L, file.exists(file.path(run_dir, "manifest.json")), dir.exists(file.path(run_dir, "experiment/methods/betageom")))
# Emulate a run made before the per-method summaries existed.
stopifnot(file.exists(file.path(run_dir, "plots/native_stats.tsv"))); unlink(Sys.glob(file.path(run_dir, "plots/native_*.tsv")))

essentiality_normalize_ui <- function() tags$div("MidLC placeholder")
ui_html <- as.character(normalization_module_ui())
# The method chooser lists readable labels against method ids, with BGC chosen; unnormalized is never a choice.
stopifnot(grepl('value="betageom" selected', ui_html, fixed = TRUE), grepl("selectize", ui_html, fixed = TRUE),
          grepl("Beta-geometric correction (BGC)", ui_html, fixed = TRUE), !grepl('value="quantile" selected', ui_html, fixed = TRUE),
          !grepl('value="nonorm"', ui_html, fixed = TRUE), !grepl("baseline)", ui_html, fixed = TRUE))   # each panel already shows unnormalized counts
# The tab is "Normalization": not every method it offers works at site level (MidLC thins whole libraries).
stopifnot(grepl('data-value="site"', ui_html, fixed = TRUE), !grepl("Site Normalization", ui_html, fixed = TRUE),
          !grepl("Site normalization", ui_html, fixed = TRUE))
stopifnot(grepl("normalization_tabs", ui_html, fixed = TRUE), !grepl('data-value="midlc"', ui_html, fixed = TRUE), grepl("MidLC common saturation", ui_html, fixed = TRUE), !grepl("Classifier)", ui_html, fixed = TRUE),
          grepl("norm_download_pdf", ui_html, fixed = TRUE), grepl("12.4%", ui_html, fixed = TRUE))

started <- NULL
fake_jobs <- list(job_is_running = function() FALSE, start_job = function(command, args, project_config, stage, ...) started <<- list(args = args, stage = stage))
active <- reactive(list(project_id = "Zn_toxicity_screen", project_root = project_root, project_config = config))
testServer(function(input, output, session)
  norm <- normalization_server(input, output, session, root, active, function() config, fake_jobs, function() python), {
  session$flushReact()
  sel <- output$norm_result_selector$html
  stopifnot(grepl("Latest normalization:", sel, fixed = TRUE), grepl("8 samples", sel, fixed = TRUE), !grepl("Earlier runs", sel, fixed = TRUE),
            !grepl("norm_fixture", sel, fixed = TRUE))   # one run: no selector, no folder name
  stopifnot(nrow(norm$runs()) == 1L, identical(norm$manifest()$suggested, "betageom"), setequal(norm$manifest()$methods, c("nonorm", "betageom", "ttr")))
  session$setInputs(norm_methods = c("betageom", "ttr"), norm_plot = "native")
  d <- norm$data(); stopifnot(nrow(d$stats) == 3L * 8L, nrow(d$scaling) == 8L, all(d$scaling$method == "ttr"), nrow(d$qq) > 0, nrow(d$ranked) > 0)
  stopifnot(all(abs(d$scaling$after - 100) < 1e-6), length(unique(round(d$scaling$before, 3))) > 1)   # TTR equalizes density x trimmed mean at the target
  panels <- output$norm_native_panels$html
  stopifnot(grepl("Trimmed total reads (TTR)", panels, fixed = TRUE), grepl("Beta-geometric correction (BGC)", panels, fixed = TRUE),
            !grepl("norm_native_nonorm", panels, fixed = TRUE), !grepl("norm_native_quantile", panels, fixed = TRUE),
            grepl("empty sites included", panels, fixed = TRUE))
  for (m in c("ttr", "betageom")) stopifnot(grepl("plotly", output[[paste0("norm_native_", m)]], fixed = TRUE))
  stopifnot(grepl("Geometric quantile (log)", output$norm_native_betageom, fixed = TRUE), grepl("Rank (log)", output$norm_native_betageom, fixed = TRUE), grepl("Density × trimmed mean", output$norm_native_ttr, fixed = TRUE))
  session$setInputs(norm_methods = "ttr")    # unticking removes the panels and the table rows
  panels <- output$norm_native_panels$html; stopifnot(!grepl("norm_native_betageom", panels, fixed = TRUE), !grepl("BGC", output$norm_stats_table, fixed = TRUE))
  session$setInputs(norm_methods = c("ttr", "quantile"))   # not in this run: reported, no panel
  stopifnot(grepl("Not in this run: Quantile", output$norm_status$html, fixed = TRUE), !grepl("norm_native_quantile", output$norm_native_panels$html, fixed = TRUE))
  session$setInputs(norm_methods = c("betageom", "ttr"), norm_plot = "pairs")
  plot_json <- output$norm_plotly; stopifnot(grepl("BGC", plot_json, fixed = TRUE), grepl("Unnormalized", plot_json, fixed = TRUE))   # kept as the reference
  stopifnot(grepl("Nonzero median", output$norm_stats_table, fixed = TRUE), grepl("Skew", output$norm_stats_table, fixed = TRUE))
  for (id in c("norm_download_png", "norm_download_pdf", "norm_download_stats", "norm_download_pairs", "norm_native_png_betageom", "norm_native_pdf_ttr")) {
    path <- output[[id]]; stopifnot(file.size(path) > 100) }
  session$setInputs(norm_trim = .05, norm_target = 100, norm_run = 1)   # a run normalizes with every method
  stopifnot(identical(started$stage, "site_normalization"), !("--method" %in% started$args), grepl("ytab_normalize_project.py$", started$args[[1L]]))
})
legacy <- navigation_target_from_hash("#tab=essentiality&subtab=normalize")
old_midlc <- navigation_target_from_hash("#tab=normalization&subtab=midlc")
stopifnot(old_midlc$top == "normalization", old_midlc$nested_value == "site")
stopifnot(legacy$top == "normalization", legacy$nested_value == "site", navigation_state_defaults()$active_essentiality_subtab == "run")
cat("PASS\n")
