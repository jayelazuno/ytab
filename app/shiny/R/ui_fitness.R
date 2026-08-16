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
    "Configure fitness analysis",
    tags$div(
      class = "alert alert-info",
      tags$b("Fitness analysis uses raw per-sample SummaryTable outputs and performs CPM normalization inside the R analysis."),
      tags$p("MidLC-normalized inputs are not used for this workflow.")
    ),
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
                 choices = list("Preview analysis" = "preview", "Run fitness analysis" = "run"),
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

fitness_results_ui <- function() tagList(
  panel_card(
    "Fitness results",
    uiOutput("fitness_results_state"),
    ytab_two_column_layout(
      controls = ytab_control_panel(
        "Fitness display",
        uiOutput("fitness_result_selector"),
        uiOutput("fitness_visualization_selector"),
        ytab_plot_customization_controls("fitness", include_points = TRUE,
                                         include_labels = TRUE,
                                         default_height = "large"),
        tags$details(class = "ytab-more-options", tags$summary("Result summary"),
                     uiOutput("fitness_smoke_warning"), uiOutput("fitness_summary_cards")),
        tags$details(class = "ytab-technical-details", tags$summary("Technical details"),
                     uiOutput("fitness_result_technical"))
      ),
      main = tagList(
        ytab_plot_card("Fitness visualization", uiOutput("fitness_selected_visualization"),
                       description = "App-rendered plots are live views from existing result tables. Generated PNGs are shown as static image cards."),
        tags$details(
          class = "ytab-more-options",
          tags$summary("Tables / downloads"),
          tags$div(
            class = "ytab-filter-row",
            selectInput("fitness_call_filter", "Fitness call", choices = "All"),
            textInput("fitness_gene_search", "Search gene or feature ID")
          ),
          uiOutput("fitness_filtered_count"),
          DT::DTOutput("fitness_results_table"),
          tags$details(tags$summary("Comparison-level results"),
                       DT::DTOutput("fitness_comparison_results")),
          tags$div(
            class = "ytab-actions",
            downloadButton("download_fitness_full", "Full fitness result"),
            downloadButton("download_fitness_filtered", "Filtered fitness result"),
            downloadButton("download_fitness_design", "Comparison design"),
            downloadButton("download_fitness_comparison", "Comparison-level summary"),
            downloadButton("download_fitness_manifest", "Run manifest")
          )
        )
      )
    )
  )
)

fitness_ui <- function() navset_tab(
  id = "fitness_tabs",
  nav_panel("Design", value = "fitness_design", fitness_design_ui()),
  nav_panel("Run Analysis", value = "fitness_run", fitness_run_ui()),
  nav_panel("Results", value = "fitness_results", fitness_results_ui())
)
