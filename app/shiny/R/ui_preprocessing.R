keep_going_supported_stage <- function(stage) stage %in% c(
  "preprocessing", "mapfastq", "create_hit_file", "summary", "sample_normalization",
  "summary_normalized", "fitness_analysis"
)

mapping_read_arguments <- function(mate="project", agreement="project") {
  stopifnot(mate %in% c("project","R1","R2"), agreement %in% c("project","yes","no"))
  c(if(mate %in% c("R1","R2")) c("--junction-mate",mate),
    if(agreement=="yes") "--require-mate-alignment",
    if(agreement=="no") "--no-require-mate-alignment")
}

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

run_profiles_ui <- function(detected_cpu) panel_card("Run profiles",
  tags$p(class="text-muted","Run several stages in order with one command."),
  selectInput("pipeline_profile","Profile",c("core","classifier","fitness","all"),"all"), sliderInput("pipeline_threads","Threads",1,max(2,detected_cpu),2),
  textInput("pipeline_target","Classifier target","recommended"), textInput("pipeline_analysis_id","Fitness analysis ID","treated_vs_parent_raw_cpm"),
  actionButton("preview_pipeline","Preview plan"), actionButton("run_pipeline","Run selected profile",class="btn-primary"), actionButton("cancel_pipeline","Cancel current job"), tags$details(class="ytab-technical-details",tags$summary("Technical details"),verbatimTextOutput("pipeline_job_text")))

run_module_ui <- function(detected_cpu = 2L) navset_tab(id="run_tabs",
  nav_panel("Samples",value="samples",panel_card("Preprocessing samples",uiOutput("project_created_card"),sample_selector_ui("preprocessing_samples"))),
  nav_panel("Reference",value="reference",panel_card("Reference readiness",uiOutput("reference_readiness_ui"),uiOutput("reference_actions_ui"))),
  # Mapping, hit files and SummaryTable run as one step (ytab_run_preprocessing.py); stages already done are skipped.
  nav_panel("Mapping",value="mapping",panel_card("Mapping",
    tags$p(class="text-muted","Maps reads, creates hit files and builds SummaryTables for the selected samples in one run. ",
      "Stages already done are skipped, including ones run elsewhere: a hit file or SummaryTable made on the HPC counts as done."),
    uiOutput("preprocessing_progress_ui"),job_progress_ui("mapping_job"),
    numericInput("preprocessing_threads","CPU threads",2,min=2,step=1),tags$p(class="text-muted","Use two threads on systems with approximately 4–8 GB RAM."),
    selectInput("mapping_junction_mate","Hermes/junction read",c("Use project settings (default R1)"="project","Read 1"="R1","Read 2"="R2")),
    selectInput("mapping_mate_agreement","Require both mates to agree",c("Use project settings (default no)"="project","No — junction read only"="no","Yes — concordant pairs"="yes")),
    tags$p(class="text-muted","Detected Hermes and adaptor sequence is trimmed automatically. Reads already starting in genomic DNA are retained without requiring a Hermes match."),
    stage_run_controls("preprocess","Run all included samples"),uiOutput("analysis_ready_action"),
    stage_status_details("Review status by sample","mapping_status_table"),
    tags$details(class="ytab-technical-details",tags$summary("Technical details"),verbatimTextOutput("pipeline_job_text_preprocessing"),verbatimTextOutput("command_output")))),
  nav_panel("Run Profiles",value="profiles",run_profiles_ui(detected_cpu)))
