#!/usr/bin/env Rscript
# Run the R smoke tests with the arguments each one expects.
#
# They do not share a convention: some take REPO_ROOT, some a project.yaml, a few take two paths, and several
# shell out to `python`. Run by hand with the wrong argument -- or with a `python` that is not the project
# environment -- a healthy test fails, which is how a third of them came to look broken. This runner holds the
# right invocation for each, so a failure here is a real one.
#
#   Rscript scripts/local/ytab_run_smoke_tests.R [PROJECT_YAML]
#
# PROJECT_YAML defaults to the local smoke project. Tests needing a fixture this workspace does not have are
# skipped by name, with the reason, rather than counted as failures.

args <- commandArgs(trailingOnly = TRUE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]]), winslash = "/")
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
project <- if (length(args) >= 1L) normalizePath(args[[1]], winslash = "/") else
  file.path(root, "output/projects/local_glabrata_smoke_v1/config/project.yaml")

# The tests that shell out must find the project's interpreter, not whatever `python` the shell resolves first
# (usually a base environment without the dependencies).
env_bin <- file.path(path.expand("~"), "miniforge3", "envs", "ytab-local", "bin")
if (dir.exists(env_bin)) Sys.setenv(PATH = paste(env_bin, Sys.getenv("PATH"), sep = .Platform$path.sep))

takes_repo_root <- c(
  "ytab_smoke_analysis_empty_states.R",
  "ytab_smoke_app_visualizations.R",
  "ytab_smoke_comparative_resources.R",
  "ytab_smoke_comparisons_ui.R",
  "ytab_smoke_downstream_status.R",
  "ytab_smoke_essentiality_job_filter.R",
  "ytab_smoke_essentiality_module_ui.R",
  "ytab_smoke_fitness_execution_mode.R",
  "ytab_smoke_fitness_job_status.R",
  "ytab_smoke_fitness_plot_edge_cases.R",
  "ytab_smoke_gene_domain_explorer_ui.R",
  "ytab_smoke_gene_explorer_ui_state.R",
  "ytab_smoke_glabrata_display_lookup.R",
  "ytab_smoke_go_analysis_ui.R",
  "ytab_smoke_help_center.R",
  "ytab_smoke_hermes_mapping.R",
  "ytab_smoke_job_progress.R",
  "ytab_smoke_job_progress_close.R",
  "ytab_smoke_job_reconciliation.R",
  "ytab_smoke_json_serialization.R",
  "ytab_smoke_landing_ui.R",
  "ytab_smoke_library_correlation_ui.R",
  "ytab_smoke_navigation_state.R",
  "ytab_smoke_normalization_ui.R",
  "ytab_smoke_original_poisson_planner.R",
  "ytab_smoke_plot_number_format.R",
  "ytab_smoke_portable_project.R",
  "ytab_smoke_process_paths.R",
  "ytab_smoke_project_discovery.R",
  "ytab_smoke_qc_execution_mode.R",
  "ytab_smoke_qc_output_presentation.R",
  "ytab_smoke_run_all_semantics.R",
  "ytab_smoke_run_mapping.R",
  "ytab_smoke_sample_names.R",
  "ytab_smoke_selectize_requests.R",
  "ytab_smoke_summary_qc_group_labels.R",
  "ytab_smoke_treated_vs_control_ui.R",
  NULL)
takes_project <- c(
  "ytab_smoke_app_initialization.R",
  "ytab_smoke_app_shell.R",
  "ytab_smoke_comparative_species_ui.R",
  "ytab_smoke_diagnostic_plots.R",
  "ytab_smoke_essentiality_layout.R",
  "ytab_smoke_fitness_design.R",
  "ytab_smoke_fitness_result_state.R",
  "ytab_smoke_fitness_selected_comparisons.R",
  "ytab_smoke_library_diagnostics_ui.R",
  "ytab_smoke_preprocessing_selection_state.R",
  "ytab_smoke_qc_result_state.R",
  "ytab_smoke_qc_sample_selector.R",
  "ytab_smoke_qc_tables.R",
  "ytab_smoke_reactive_context.R",
  "ytab_smoke_sample_selector.R",
  "ytab_smoke_workspace_state.R",
  NULL)
# name -> why it cannot run here. These are fixture gaps, not failures.
needs_fixture <- c(
  "ytab_smoke_essentiality_state.R" = "needs a project with four parent libraries and a completed evaluate stage",
  "ytab_smoke_essentiality_commands.R" = "asserts a repository path containing a space",
  "ytab_smoke_fitness_h2o2_test.R" = "runs the fitness CLI end to end; needs its H2O2 fixture project",
  "ytab_smoke_new_project_readiness.R" = "needs a FRESH project whose next stage is Mapping",
  "ytab_smoke_preprocessing_commands.R" = "needs a fresh project",
  "ytab_smoke_project_entry_modes.R" = "needs its own fixture projects",
  "ytab_smoke_navigation.R" = "needs a fresh project (opens one through the server)",
  "ytab_smoke_open_existing_project.R" = "needs a fresh project",
  "ytab_smoke_qc_module.R" = "needs QC summary rows for a single library",
  "ytab_smoke_tvc_agreement.R" = "needs a treated-vs-control run with calls in both directions",
  "ytab_smoke_rnacross_import.R" = "needs codex/RNAcross, which is not in this workspace",
  "ytab_smoke_job_manager_integration.R" = "starts a background worker; flaky outside a job host",
  "ytab_smoke_mutant_identification.R" = "needs REPO_ROOT plus a fixture directory",
  "ytab_smoke_screen_coverage_project_switch.R" = "needs REPO_ROOT plus a second project",
  NULL)

run_one <- function(name) {
  script <- file.path(root, "scripts/local", name)
  if (!file.exists(script)) return("MISSING")
  arg <- if (name %in% takes_repo_root) root else if (name %in% takes_project) project else character()
  out <- suppressWarnings(system2("Rscript", c(script, arg), stdout = TRUE, stderr = TRUE))
  if (length(out) && grepl("^PASS", utils::tail(out, 1L))) "PASS" else "FAIL"
}

all_tests <- sort(basename(list.files(file.path(root, "scripts/local"), "ytab_smoke_.*[.]R$")))
pass <- fail <- skip <- character()
for (name in all_tests) {
  if (name %in% names(needs_fixture)) { skip <- c(skip, name); next }
  result <- run_one(name)
  if (identical(result, "PASS")) pass <- c(pass, name) else fail <- c(fail, name)
  cat(sprintf("%-6s %s\n", result, name))
}
cat(sprintf("\n%d passed, %d failed, %d skipped (fixture not in this workspace)\n",
            length(pass), length(fail), length(skip)))
for (name in skip) cat(sprintf("  skipped  %-44s %s\n", name, needs_fixture[[name]]))
if (length(fail)) { cat("\nfailed:\n"); for (name in fail) cat("  ", name, "\n"); quit(status = 1L) }
quit(status = 0L)
