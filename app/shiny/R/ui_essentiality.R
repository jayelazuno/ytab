essentiality_stage_labels <- c(
  normalize = "Normalize & Choose Target",
  evaluate = "Evaluate",
  combine = "Combine",
  combined_summary = "Combined Summary",
  classifier = "Classifier",
  results = "Results"
)

essentiality_stage_backend <- c(
  normalize = "sample_normalization",
  evaluate = "summary_normalized",
  combine = "combined_hits",
  combined_summary = "summary_combined",
  classifier = "classifier"
)

essentiality_stage_icons <- c(
  blocked = "\u26d4", ready = "\u25cb", running = "\u25d4",
  complete = "\u2713", cached = "\u21ba", failed = "\u26a0"
)

essentiality_execution_controls <- function(prefix) {
  tagList(
    radioButtons(
      paste0(prefix, "_execution_mode"), "Execution mode",
      choices = c("Preview command" = "preview", "Run stage" = "run"),
      selected = "preview", inline = TRUE
    ),
    uiOutput(paste0(prefix, "_action"))
  )
}

essentiality_technical_details <- function(id) uiOutput(id)

essentiality_stage_job_slot <- function(prefix) {
  uiOutput(paste0(prefix, "_job_disclosure"))
}

essentiality_review_details <- function(label, ...) {
  tags$details(class = "ytab-more-options", tags$summary(label), ...)
}

essentiality_run_panel_ui <- function(completed, rerun_label, ready_title, content) {
  if (isTRUE(completed))
    tags$details(
      class = "ytab-rerun-disclosure",
      tags$summary(rerun_label),
      content
    )
  else tags$section(
    class = "ytab-run-card",
    tags$h4(ready_title),
    content
  )
}

essentiality_normalization_configuration_ui <- function() {
  tagList(
    radioButtons(
      "normalization_target_mode", "Target mode",
      choices = c("Auto scan and recommend" = "auto",
                  "Manual target list" = "manual"),
      selected = "auto", inline = TRUE
    ),
    conditionalPanel(
      condition = "input.normalization_target_mode == 'manual'",
      textInput("normalization_manual_targets", "Targets", "20,20.5,59.7,100"),
      helpText("Enter one or more positive targets, such as 20, 20.5, 59.7, 100."),
      DT::DTOutput("normalization_target_preview")
    ),
    conditionalPanel(
      condition = "input.normalization_target_mode == 'auto'",
      tags$details(
        tags$summary("Advanced scan settings"),
        tags$div(
          class = "ytab-config-grid",
          numericInput("normalization_auto_min", "Minimum target", value = NA, min = 0),
          numericInput("normalization_auto_max", "Maximum target", value = NA, min = 0),
          numericInput("normalization_auto_step", "Target step", value = 5, min = .000001),
          numericInput("normalization_min_retention", "Minimum site retention",
                       value = .95, min = .01, max = 1, step = .01)
      )
    )
  ),
    essentiality_review_details("Run details", uiOutput("normalization_prerun_summary")),
    essentiality_execution_controls("normalization")
  )
}

essentiality_evaluation_configuration_ui <- function() {
  tagList(
    radioButtons(
      "summary_normalized_target_source", "Target source",
      choices = c("All normalized targets" = "all",
                  "Recommended target" = "recommended",
                  "Selected existing targets" = "existing",
                  "Manual target list" = "manual"),
      selected = "all"
    ),
    conditionalPanel(
      condition = "input.summary_normalized_target_source == 'existing'",
      selectizeInput("summary_normalized_existing_targets", "Existing targets",
                     choices = character(), multiple = TRUE)
    ),
    conditionalPanel(
      condition = "input.summary_normalized_target_source == 'manual'",
      tags$details(
        tags$summary("Advanced settings"),
        textInput("summary_normalized_manual_targets", "Manual targets", "59.7")
      )
    ),
    numericInput("summary_normalized_min_retention", "Minimum feature retention",
                 value = .95, min = .01, max = 1, step = .01),
    essentiality_review_details("Run details", uiOutput("summary_normalized_prerun_summary")),
    essentiality_execution_controls("summary_normalized")
  )
}

essentiality_classifier_configuration_ui <- function() {
  tagList(
    tags$details(
      tags$summary("Advanced settings"),
      textInput("classifier_seed", "Random seed", "0"),
      helpText("The seed controls reproducibility and does not change the input target.")
    ),
    essentiality_review_details("Run details", uiOutput("classifier_prerun_summary")),
    essentiality_execution_controls("classifier")
  )
}

essentiality_normalize_ui <- function() {
  tagList(
    panel_card(
      "Normalize & Choose Target",
      tags$p(class = "ytab-stage-purpose",
             title = "The normalization target is a MidLC depth target, not a percent. YTAB down-samples libraries only when needed and never upsamples.",
             "Choose the recommended MidLC target, then continue."),
      tags$section(
        class = "ytab-review-section",
        uiOutput("normalization_recommendation_card"),
        uiOutput("essentiality_target_selection")
      ),
      tags$details(class="ytab-more-options",tags$summary("Run workflow steps"),
        uiOutput("normalization_run_panel"),
        uiOutput("summary_normalized_run_panel")),
      tags$details(class="ytab-more-options",tags$summary("Visualizations"),
        uiOutput("target_comparison_ui")),
      tags$details(class="ytab-more-options",tags$summary("Review details"),
        tags$h4("Parent libraries"),
        DT::DTOutput("essentiality_parent_table"),
        uiOutput("essentiality_parent_override"),
        uiOutput("normalization_rerun_panel"),
        uiOutput("summary_normalized_rerun_panel"),
        essentiality_stage_job_slot("normalization"),
        essentiality_stage_job_slot("summary_normalized")),
      tags$details(class="ytab-more-options",tags$summary("Technical details"),
        essentiality_technical_details("normalization_technical"),
        essentiality_technical_details("summary_normalized_technical"))
    )
  )
}

essentiality_combine_ui <- function() {
  tagList(
    panel_card(
      "Combine normalized parent libraries",
      tags$p(class = "ytab-stage-purpose",
             "Merge the normalized parent hit files for the selected target into one parent library."),
      uiOutput("combined_hits_status"),
      uiOutput("combine_prerequisite_summary"),
      uiOutput("combine_run_panel"),
      tags$details(class="ytab-more-options",tags$summary("Review details"),
        uiOutput("combined_hits_result"),
        tags$h4("Parent inputs"),DT::DTOutput("combine_parent_inputs"),
        uiOutput("combine_rerun_panel"),essentiality_stage_job_slot("combine"),
        essentiality_technical_details("combine_technical"))
    )
  )
}

essentiality_combined_summary_ui <- function() {
  tagList(
    panel_card(
      "Build combined parent feature table",
      tags$p(class = "ytab-stage-purpose",
             "Build the classifier-ready feature table from the combined normalized parent library."),
      uiOutput("combined_summary_status"),
      uiOutput("combined_summary_prerequisite"),
      uiOutput("combined_summary_run_panel"),
      tags$details(class="ytab-more-options",tags$summary("Review details"),
        uiOutput("combined_summary_result"),
        uiOutput("combined_summary_rerun_panel"),essentiality_stage_job_slot("combined_summary"),
        essentiality_technical_details("combined_summary_technical"))
    )
  )
}

essentiality_classifier_ui <- function() {
  tagList(
    panel_card(
      "Essentiality classifier",
      tags$p(class = "ytab-stage-purpose",
             "Apply the established classifier to the combined parent feature table for the selected normalization target."),
      uiOutput("classifier_species_support"),
      uiOutput("classifier_result_summary"),
      uiOutput("classifier_run_panel"),
      tags$details(class="ytab-more-options",tags$summary("Review details"),
        uiOutput("classifier_readiness_summary"),
        uiOutput("classifier_recommended_target"),
        tags$h4("Classifier resources"),DT::DTOutput("classifier_resource_table"),
        uiOutput("classifier_rerun_panel"),essentiality_stage_job_slot("classifier"),
        essentiality_technical_details("classifier_technical"))
    )
  )
}

essentiality_results_overview_ui <- function(has_label = TRUE) {
  tagList(
    uiOutput("classifier_summary_cards"),
    if (isTRUE(has_label))
      tags$div(class = "ytab-plot-card ytab-overview-plot",
               tags$h4("Prediction-label distribution"),
               plotOutput("classifier_label_plot")),
    uiOutput("classifier_final_target_action")
  )
}

essentiality_results_predictions_ui <- function() {
  tagList(
    tags$div(
      class = "ytab-filter-row",
      textInput("classifier_gene_search", "Search gene or feature", ""),
      selectInput("classifier_label_filter", "Classifier label",
                  choices = c("All labels" = "All")),
      selectInput("classifier_inclusion_filter", "Features",
                  choices = c("All" = "All", "Included" = "Included",
                              "Excluded" = "Excluded"))
    ),
    uiOutput("classifier_filtered_count"),
    DT::DTOutput("classifier_predictions_table")
  )
}

essentiality_results_visualizations_ui <- function(columns) {
  tagList(uiOutput("essentiality_visualization_selector"),
          uiOutput("essentiality_selected_visualization"))
}

essentiality_results_ui <- function() {
  tagList(
    panel_card(
      "Classifier results",
      tags$p(class = "ytab-stage-purpose",
             "Search classifier predictions for the selected normalization target."),
      uiOutput("essentiality_smoke_warning"),
      uiOutput("classifier_results_empty"),
      ytab_two_column_layout(
        controls = ytab_control_panel(
          "Classifier display",
          uiOutput("classifier_result_selector"),
          uiOutput("classifier_result_state"),
          uiOutput("essentiality_visualization_selector"),
          ytab_plot_customization_controls("essentiality", include_bars = TRUE,
                                           include_heatmap = TRUE,
                                           default_height = "medium"),
          tags$details(class="ytab-more-options",tags$summary("Result summary"),
            uiOutput("classifier_summary_cards")),
          tags$details(class="ytab-more-options",tags$summary("Target provenance"),
            uiOutput("classifier_provenance"),
            essentiality_technical_details("classifier_results_technical"))
        ),
        main = tagList(
          ytab_plot_card("Classifier visualization", uiOutput("essentiality_selected_visualization"),
                         description = "App-rendered plots use classifier result tables; generated classifier images are shown as static provenance images."),
          tags$details(class="ytab-more-options",tags$summary("Tables / downloads"),
            essentiality_results_predictions_ui(),
            uiOutput("essentiality_download_cards"))
        )
      )
    )
  )
}

essentiality_ui <- function() {
  navset_tab(
    id = "essentiality_tabs",
    nav_panel("Normalize & Choose Target", value = "normalize",
              essentiality_normalize_ui()),
    nav_panel("Combine Parents", value = "combine", essentiality_combine_ui()),
    nav_panel("Combined Summary", value = "combined_summary",
              essentiality_combined_summary_ui()),
    nav_panel("Classifier", value = "classifier", essentiality_classifier_ui()),
    nav_panel("Results", value = "results", essentiality_results_ui())
  )
}
