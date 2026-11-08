fitness_design_ui <- function() tagList(
  panel_card(
    "Comparison design",
    uiOutput("fitness_design_state"),
    uiOutput("fitness_design_cards"),
    tags$div(
      class = "ytab-actions",
      uiOutput("fitness_design_generate_action"),
      actionButton("fitness_refresh_design", "Refresh design", class = "btn-secondary"),
      actionButton("fitness_review_issues", "Review design issues", class = "btn-secondary"),
      actionButton("fitness_continue_analysis", "Continue to analysis", class = "btn-primary")
    ),
    uiOutput("fitness_design_presets"),
    DT::DTOutput("fitness_design_table"),
    tags$details(id = "fitness_design_issues_details", tags$summary("Design issues"),
                 DT::DTOutput("fitness_design_issues_table")),
    tags$details(tags$summary("Show input details"), DT::DTOutput("fitness_design_input_details"))
  )
)

fitness_run_ui <- function() tagList(
  panel_card(
    "Local-abundance z-score",
    uiOutput("fitness_normalization_selector"),
    helpText("BGC from the latest normalization run is the default; gene values are the normalized site counts summed per gene, per million. Total reads (CPM) needs no normalization run."),
    selectizeInput("fitness_comparisons", "Comparisons to analyze", choices = character(),
                   multiple = TRUE, options = list(placeholder = "Select valid comparisons…")),
    uiOutput("fitness_selection_summary"),
    tags$div(
      class = "ytab-actions",
      actionButton("fitness_select_all", "Select all valid", class = "btn-secondary"),
      actionButton("fitness_clear_selection", "Clear selection", class = "btn-secondary"),
      actionButton("fitness_change_design", "Change design", class = "btn-secondary")
    ),
    tags$details(tags$summary("Review selected comparisons"), DT::DTOutput("fitness_selected_preview")),
    radioButtons("fitness_execution_mode", "Execution mode",
                 choices = list("Preview analysis" = "preview", "Run analysis" = "run"),
                 selected = "preview", inline = TRUE),
    checkboxInput("fitness_annotate_classifier", "Annotate with classifier", FALSE),
    uiOutput("fitness_classifier_controls"),
    tags$details(tags$summary("Advanced settings"),
                 textInput("fitness_analysis_id", "Analysis ID", "treated_vs_parent_raw_cpm")),
    uiOutput("fitness_run_summary"),
    uiOutput("fitness_run_action"),
    job_progress_ui("fitness_job"),
    tags$details(class = "ytab-technical-details", tags$summary("Technical details"),
                 tags$div(class = "ytab-technical-console", verbatimTextOutput("fitness_technical")))
  )
)


