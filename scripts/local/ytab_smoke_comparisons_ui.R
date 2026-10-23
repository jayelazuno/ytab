#!/usr/bin/env Rscript
# Condition Comparisons: on the Zn screen runs the real chain normalize-project -> multicondition-project ->
# genetic-interactions, then checks both panels: inputs from the normalization run, condition and strain assignment,
# summaries (incl. the GI clone-sharing warning), interactive plots, tables, downloads and the run commands.
args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
suppressPackageStartupMessages({library(shiny); library(bslib)})
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
for (name in c("job_manager.R", "project_discovery.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R",
               "table_display_helpers.R", "analysis_runs.R", "normalization_module.R", "ui_fitness.R", "treated_vs_control.R", "comparisons_module.R"))
  source(file.path(root, "app/shiny/R", name))
python <- unname(Sys.which("python")); config <- file.path(root, "output/projects/Zn_toxicity_screen/config/project.yaml")
cli <- function(script, ...) { status <- system2(python, c(file.path(root, "scripts/local", script), ...), stdout = FALSE, stderr = FALSE); stopifnot(status == 0L) }
project <- read_project_summary(config, root)
project_root <- file.path(tempfile("cmp_project_"), "Zn"); dir.create(project_root, recursive = TRUE)
norm <- file.path(project_root, "normalization", "norm_fixture")
cli("ytab_normalize_project.py", "--project-config", config, "--method", "betageom", "--output", norm)
obs <- utils::read.delim(file.path(norm, "experiment/matrix/observations.tsv"), stringsAsFactors = FALSE)
condition <- ifelse(grepl("mock$", obs$observation_id), "mock", "Zn"); strain <- sub("::.*", "", sub("^[^:]*::", "", obs$library_id))
cli("ytab_multicondition_project.py", "--project-config", config, "--normalized", file.path(norm, "experiment/methods"), "--normalization", "betageom",
    as.vector(rbind("--sample", paste0(obs$observation_id, "=", condition))), "--reference-condition", "mock", "--stratify-libraries",
    "--output", file.path(project_root, "multicondition", "mc_fixture"))
gff <- file.path(root, "resources/species/glabrata/reference_genome/GCA_014217725.1_ASM1421772v1_genomic.gff"); fasta <- sub("\\.gff$", ".fna", gff)
cli("ytab_genetic_interactions.py", "--normalized", file.path(norm, "experiment/methods"), "--normalization", "betageom", "--gff", gff, "--reference", fasta,
    "--output", file.path(project_root, "genetic_interactions", "gi_fixture"), "--reference-strain", "yH298", "--query-strain", "yH299",
    as.vector(rbind("--sample", paste0(obs$observation_id, "=", strain, ":", obs$role))))

fitness_design_ui <- fitness_run_ui <- fitness_results_ui <- function() tags$div("fitness placeholder")
ui_html <- as.character(comparisons_module_ui())
stopifnot(grepl("Genetic Interactions", ui_html, fixed = TRUE), grepl("Multi-condition", ui_html, fixed = TRUE), grepl("fitness_tabs", ui_html, fixed = TRUE),
          grepl("not validated", ui_html, fixed = TRUE), grepl("gi_download_genes", ui_html, fixed = TRUE))

started <- NULL
fake_jobs <- list(job_is_running = function() FALSE, start_job = function(command, args, project_config, stage, ...) started <<- list(args = args, stage = stage))
active <- reactive(modifyList(project, list(project_root = project_root)))
testServer(function(input, output, session)
  cmp <- comparisons_server(input, output, session, root, active, function() config, fake_jobs, function() python), {
  session$flushReact()
  stopifnot(nrow(cmp$norm_runs()) == 1L, nrow(cmp$mc_runs()) == 1L, nrow(cmp$gi_runs()) == 1L)
  # Multi-condition: inputs, assignment from the sample sheet, run command
  session$setInputs(mc_norm_run = norm, mc_method = "betageom", mc_condition_column = "treatment")
  stopifnot(grepl("Normalization method", output$mc_method_selector$html, fixed = TRUE), grepl("mc_cond_8", output$mc_condition_table$html, fixed = TRUE))
  values <- ifelse(grepl("mock$", obs$observation_id), "mock", "1_5mM-Zn")
  do.call(session$setInputs, stats::setNames(as.list(values), paste0("mc_cond_", seq_along(values))))
  session$setInputs(mc_reference = "mock", mc_stratify = TRUE, mc_test = "both", mc_fdr = .05, mc_run = 1)
  stopifnot(setequal(cmp$mc_assignment()$condition, c("mock", "1_5mM-Zn")), identical(started$stage, "multicondition"), "--stratify-libraries" %in% started$args,
            sum(started$args == "--sample") == 8L, grepl("ytab_multicondition_project.py$", started$args[[1L]]))
  # Multi-condition results
  stopifnot(grepl("ANOVA: genes that vary", output$mc_summary$html, fixed = TRUE))
  session$setInputs(mc_plot = "volcano", mc_pair = "Zn vs mock"); stopifnot(grepl("plotly", output$mc_plotly, fixed = TRUE), grepl("q_all_pairs", output$mc_table, fixed = TRUE))
  session$setInputs(mc_plot = "correlation"); stopifnot(grepl("two non-reference conditions", output$mc_plot_frame$html, fixed = TRUE))
  session$setInputs(mc_plot = "volcano")
  for (id in c("mc_download_png", "mc_download_pdf", "mc_download_anova", "mc_download_pairs")) stopifnot(file.size(output[[id]]) > 100)
  # Genetic interactions: strains from the sample sheet, run command
  session$setInputs(gi_norm_run = norm, gi_method = "betageom", gi_strain_column = "background", gi_reference = "yH298", gi_query = "yH299",
                    gi_fdr = .05, gi_rope = .5, gi_correction = "bh", gi_confidence = .95, gi_run = 1)
  stopifnot(identical(started$stage, "genetic_interactions"), any(grepl("=yH298:control$", started$args)), any(grepl("=yH299:treated$", started$args)),
            identical(started$args[[which(started$args == "--rope") + 1L]], "0.5"))
  # GI results: Zn control and treated cultures share no clones, so the low-power warning must show
  stopifnot(grepl("Genes tested", output$gi_summary$html, fixed = TRUE), grepl("Low power", output$gi_summary$html, fixed = TRUE))
  for (choice in c("scatter", "forest", "volcano")) { session$setInputs(gi_plot = choice); stopifnot(grepl("plotly", output$gi_plotly, fixed = TRUE)) }
  stopifnot(grepl("delta_lower", output$gi_table, fixed = TRUE))
  for (id in c("gi_download_png", "gi_download_pdf", "gi_download_genes")) stopifnot(file.size(output[[id]]) > 100)
})
cat("PASS\n")
