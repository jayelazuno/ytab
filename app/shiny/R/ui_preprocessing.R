keep_going_supported_stage <- function(stage) stage %in% c(
  "mapfastq", "create_hit_file", "summary", "sample_normalization",
  "summary_normalized", "fitness_analysis"
)

resolve_keep_going <- function(stage, item_count, run_all = FALSE,
                               selected_value = TRUE) {
  keep_going_supported_stage(stage) && as.integer(item_count) > 1L &&
    (isTRUE(run_all) || isTRUE(selected_value))
}

keep_going_argument <- function(stage, item_count, run_all = FALSE,
                                selected_value = TRUE) {
  if (resolve_keep_going(stage, item_count, run_all, selected_value))
    "--keep-going" else character()
}

run_scope_message <- function(run_all = FALSE) {
  if (isTRUE(run_all))
    "All eligible items will be processed. Failures are recorded and remaining items continue."
  else
    "Selected items will be processed. For multi-item selections, remaining items continue after a recorded failure."
}

stage_run_controls <- function(prefix,all_label) tagList(
  uiOutput(paste0(prefix,"_selected_names")),
  tags$div(class="ytab-actions",actionButton(paste0("change_",prefix,"_selection"),"Change sample selection")),
  uiOutput(paste0(prefix,"_run_scope_help")),
  tags$div(class="ytab-actions",actionButton(paste0("preview_",prefix),"Preview selected"),actionButton(paste0("run_",prefix),"Run selected samples",class="btn-primary"),actionButton(paste0("run_all_",prefix),all_label),actionButton(paste0("refresh_",prefix,"_status"),"Refresh status")))

stage_status_details <- function(label, output_id) {
  tags$details(
    tags$summary(label),
    DT::DTOutput(output_id)
  )
}

preprocessing_ui <- function() navset_tab(id="preprocessing_tabs",
  nav_panel("Samples",value="samples",panel_card("Preprocessing samples",uiOutput("project_created_card"),sample_selector_ui("preprocessing_samples"))),
  nav_panel("Reference",value="reference",panel_card("Reference readiness",uiOutput("reference_readiness_ui"),uiOutput("reference_actions_ui"))),
  nav_panel("Mapping",value="mapping",panel_card("Mapping",job_progress_ui("mapping_job"),numericInput("preprocessing_threads","CPU threads",2,min=2,step=1),tags$p(class="text-muted","Use two threads on systems with approximately 4–8 GB RAM."),stage_run_controls("mapfastq","Run all included samples"),stage_status_details("Review mapping status by sample","mapping_status_table"))),
  nav_panel("Hit Files",value="hit_files",panel_card("Create hit files",stage_run_controls("create_hit_file","Run all mapped samples"),stage_status_details("Review hit-file status by sample","hits_status_table"))),
  nav_panel("Summary Tables",value="summary_tables",panel_card("Raw SummaryTable",stage_run_controls("summary_table","Run all samples with hit files"),uiOutput("analysis_ready_action"),stage_status_details("Review SummaryTable status by sample","summary_status_table"))),
  nav_panel("Status",value="status",panel_card("Core preprocessing status",uiOutput("preprocessing_progress_ui"),DT::DTOutput("pipeline_matrix"),uiOutput("next_ready_action_ui"),tags$details(class="ytab-technical-details",tags$summary("Technical details"),verbatimTextOutput("pipeline_job_text_preprocessing"),verbatimTextOutput("command_output")))))
