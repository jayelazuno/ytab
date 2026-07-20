library(shiny)

`%||%` <- function(left, right) if (is.null(left)) right else left
app_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL) %||% file.path("app", "shiny", "app.R")
repo_root <- normalizePath(file.path(dirname(app_file), "../.."), mustWork = TRUE)

detected_cpu <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cpu)) detected_cpu <- 2L
detected_cpu <- as.integer(detected_cpu)
default_thread_count <- max(2L, min(4L, detected_cpu))
species_dir <- file.path(repo_root, "resources", "species")
species_choices <- if (dir.exists(species_dir)) {
  candidates <- list.dirs(species_dir, recursive = FALSE, full.names = TRUE)
  selectable <- vapply(candidates, function(path) {
    reference_dir <- file.path(path, "reference_genome")
    dir.exists(reference_dir) && length(list.files(
      reference_dir,
      pattern = "(\\.(fna|fa|fasta|gff|gff3|gtf)$|feature_table.*\\.txt$|chromosomal_feature.*\\.tab$)",
      ignore.case = TRUE
    )) > 0L
  }, logical(1))
  sort(basename(candidates[selectable]))
} else character()

reference_status <- function(species) {
  reference_dir <- file.path(species_dir, species, "reference_genome")
  files <- if (dir.exists(reference_dir)) list.files(reference_dir) else character()
  fasta <- any(grepl("\\.(fna|fa|fasta)$", files, ignore.case = TRUE))
  annotation <- any(grepl(
    "(\\.(gff|gff3|gtf)$|feature_table.*\\.txt$|chromosomal_feature.*\\.tab$)",
    files, ignore.case = TRUE
  ))
  prefixes <- sub("\\.1\\.(bt2|bt2l)$", "", files[grepl("\\.1\\.(bt2|bt2l)$", files)])
  complete <- any(vapply(prefixes, function(prefix) {
    ext <- if (all(file.exists(file.path(reference_dir, paste0(prefix, c(
      ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
    )))))) "bt2" else ""
    if (nzchar(ext)) return(TRUE)
    all(file.exists(file.path(reference_dir, paste0(prefix, c(
      ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"
    )))))
  }, logical(1)))
  list(runnable = fasta && annotation && complete, index = complete)
}

ui <- fluidPage(
  titlePanel("YTAB: Yeast Transposon Analysis Browser"),
  sidebarLayout(
    sidebarPanel(
      h3("Project Setup"),
      textInput("project_id", "Project ID"),
      textInput("fastq_dir", "FASTQ directory"),
      selectInput("species", "Species/reference", choices = species_choices),
      sliderInput("threads", "CPU threads", min = 2L, max = max(2L, detected_cpu),
                  value = min(default_thread_count, max(2L, detected_cpu)), step = 1L),
      helpText(sprintf("Detected CPU count: %d", detected_cpu)),
      conditionalPanel("input.threads > 4",
                       tags$div(class = "text-warning",
                                "Using more than 4 threads may increase memory use on low-memory machines.")),
      actionButton("initialize", "Initialize local project", class = "btn-primary")
    ),
    mainPanel(
      h3("Initialization output"),
      verbatimTextOutput("command_output"),
      hr(),
      h3("Reference Preparation"),
      textOutput("reference_status"),
      actionButton("prepare_reference", "Prepare reference / build index"),
      hr(),
      h3("Run Mapping"),
      textInput("mapping_project_config", "Project config path",
                value = "output/projects/local_glabrata_smoke_v1/config/project.yaml"),
      textInput("mapping_samples", "Samples (comma-separated; blank means all)", value = ""),
      checkboxInput("mapping_dry_run", "Dry run", value = TRUE),
      actionButton("run_mapfastq", "Run MapFastq"),
      hr(),
      h3("Create Hit Files"),
      textInput("hits_project_config", "Project config path",
                value = "output/projects/local_glabrata_smoke_v1/config/project.yaml"),
      textInput("hits_samples", "Samples (comma-separated; blank means all)", value = ""),
      checkboxInput("hits_dry_run", "Dry run", value = TRUE),
      actionButton("run_create_hit_file", "Run CreateHitFile"),
      hr(),
      h3("Build Summary Tables"),
      textInput("summary_project_config", "Project config path",
                value = "output/projects/local_glabrata_smoke_v1/config/project.yaml"),
      textInput("summary_samples", "Samples (comma-separated; blank means all)", value = ""),
      checkboxInput("summary_dry_run", "Dry run", value = TRUE),
      actionButton("run_summary_table", "Run SummaryTable"),
      hr(),
      h3("Library Diagnostics"),
      textInput("diagnostics_project_config", "Project config path",
                value = "output/projects/local_glabrata_smoke_v1/config/project.yaml"),
      textInput("diagnostics_samples", "Samples (comma-separated; blank means all)", value = ""),
      checkboxInput("diagnostics_dry_run", "Dry run", value = TRUE),
      actionButton("run_library_diagnostics", "Run LibraryDiagnostics"),
      hr(),
      h3("Normalization Explorer"),
      textInput("normalization_project_config", "Project config path",
                value = "output/projects/local_glabrata_smoke_v1/config/project.yaml"),
      selectInput("normalization_target_mode", "Target mode",
                  choices = c("Auto recommend target" = "auto", "Manual targets" = "manual")),
      conditionalPanel(
        "input.normalization_target_mode == 'auto'",
        numericInput("normalization_min_retention", "Minimum site retention", value = 0.95,
                     min = 0.01, max = 1, step = 0.01),
        textInput("normalization_auto_min", "Auto minimum target (optional)", value = ""),
        textInput("normalization_auto_max", "Auto maximum target (optional)", value = ""),
        numericInput("normalization_auto_step", "Auto target step", value = 5, min = 0.01)
      ),
      conditionalPanel(
        "input.normalization_target_mode == 'manual'",
        textInput("normalization_manual_targets", "Manual targets (comma-separated)",
                  value = "20,20.5,59.7,100")
      ),
      selectInput("normalization_sample_mode", "Sample mode",
                  choices = c("parents", "all", "treated"), selected = "parents"),
      checkboxInput("normalization_dry_run", "Dry run", value = TRUE),
      checkboxInput("normalization_force", "Force rerun", value = FALSE),
      actionButton("run_sample_normalization", "Run normalization"),
      helpText("Auto mode recommends a conservative target based on insertion-site retention. Feature-level retention will be confirmed after SummaryTable is run on normalized hits."),
      helpText("Final classifier target selection will be added after normalized SummaryTable and combined parent hits are implemented."),
      h4("Normalization recommendation"),
      uiOutput("normalization_recommendation_ui"),
      h4("Normalization comparison"),
      uiOutput("normalization_comparison_ui"),
      hr(),
      h3("Normalized Summary / Target Evaluation"),
      textInput("summary_normalized_project_config", "Project config path",
                value = "output/projects/local_glabrata_smoke_v1/config/project.yaml"),
      selectInput("summary_normalized_target_mode", "Targets",
                  choices = c("recommended", "all", "manual"), selected = "recommended"),
      conditionalPanel(
        "input.summary_normalized_target_mode == 'manual'",
        textInput("summary_normalized_manual_targets", "Manual targets or tags",
                  value = "T059p7")
      ),
      selectInput("summary_normalized_sample_mode", "Sample mode",
                  choices = c("parents", "all", "treated"), selected = "parents"),
      numericInput("summary_normalized_min_retention", "Minimum feature retention",
                   value = 0.95, min = 0.01, max = 1, step = 0.01),
      checkboxInput("summary_normalized_dry_run", "Dry run", value = TRUE),
      checkboxInput("summary_normalized_force", "Force rerun", value = FALSE),
      actionButton("run_summary_normalized", "Run normalized SummaryTable"),
      helpText("This step evaluates feature-level retention after normalization. Parent-pool combining and classifier input selection come next."),
      h4("Target evaluation"),
      uiOutput("normalization_target_evaluation_ui"),
      h4("Target summary"),
      uiOutput("normalization_target_summary_ui"),
      h4("Feature recommendation"),
      uiOutput("normalization_feature_recommendation_ui"),
      hr(),
      h3("Combine Parent Hits"),
      textInput("combine_hits_project_config", "Project config path",
                value = "output/projects/local_glabrata_smoke_v1/config/project.yaml"),
      selectInput("combine_hits_target_mode", "Target",
                  choices = c("recommended", "site-recommended", "feature-recommended", "manual"),
                  selected = "recommended"),
      conditionalPanel(
        "input.combine_hits_target_mode == 'manual'",
        textInput("combine_hits_manual_target", "Manual target or tag", value = "T059p7")
      ),
      textInput("combine_hits_samples", "Parent samples override (comma-separated; blank means all parents)", value = ""),
      checkboxInput("combine_hits_dry_run", "Dry run", value = TRUE),
      checkboxInput("combine_hits_force", "Force rerun", value = FALSE),
      actionButton("run_combine_hits", "Combine normalized parent hits"),
      helpText("This step combines normalized parent libraries for the selected target. SummaryTable on the combined parent library and classifier input generation come next."),
      h4("Combined hits status"),
      uiOutput("combined_hits_status_ui"),
      hr(),
      h3("Combined Parent Summary"),
      textInput("summary_combined_project_config", "Project config path",
                value = "output/projects/local_glabrata_smoke_v1/config/project.yaml"),
      selectInput("summary_combined_target_mode", "Target",
                  choices = c("recommended", "site-recommended", "feature-recommended", "manual"),
                  selected = "recommended"),
      conditionalPanel(
        "input.summary_combined_target_mode == 'manual'",
        textInput("summary_combined_manual_target", "Manual target or tag", value = "T059p7")
      ),
      checkboxInput("summary_combined_dry_run", "Dry run", value = TRUE),
      checkboxInput("summary_combined_force", "Force rerun", value = FALSE),
      actionButton("run_summary_combined", "Run combined parent SummaryTable"),
      helpText("This step creates the combined parent feature table. Classifier input generation and treated-vs-parent analysis remain separate later steps."),
      h4("Combined summary status"),
      uiOutput("summary_combined_status_ui"),
      hr(),
      h3("Essentiality Classifier"),
      textInput("classifier_project_config", "Project config path",
                value = "output/projects/local_glabrata_smoke_v1/config/project.yaml"),
      selectInput("classifier_target_mode", "Target",
                  choices = c("recommended", "feature-recommended", "site-recommended", "available", "manual"),
                  selected = "recommended"),
      conditionalPanel(
        "input.classifier_target_mode == 'available'",
        selectInput("classifier_available_target", "Available combined-summary target", choices = character(0))
      ),
      conditionalPanel(
        "input.classifier_target_mode == 'manual'",
        textInput("classifier_manual_target", "Manual target or tag", value = "T059p7")
      ),
      textInput("classifier_seed", "Random seed (blank uses classifier default)", value = ""),
      checkboxInput("classifier_dry_run", "Dry run", value = TRUE),
      checkboxInput("classifier_force", "Force rerun", value = FALSE),
      checkboxInput("classifier_save_final_target", "Save this as the final classifier target", value = FALSE),
      actionButton("run_classifier", "Run essentiality classifier"),
      helpText("Results from shallow FASTQ subsets validate the workflow but should not be interpreted as final biological essentiality calls."),
      helpText("Treated-versus-parent fitness analysis is a separate downstream step and is not run by the essentiality classifier."),
      h4("Classifier status"),
      uiOutput("classifier_status_ui"),
      h4("Classifier summary"),
      uiOutput("classifier_summary_ui"),
      h4("Essentiality predictions"),
      uiOutput("classifier_predictions_ui"),
      hr(),
      h3("Treated vs Parent Fitness"),
      textInput("tvp_project_config", "Project config path",
                value = "output/projects/local_glabrata_smoke_v1/config/project.yaml"),
      actionButton("tvp_design", "Generate / refresh comparison design"),
      h4("Comparison design"), uiOutput("tvp_design_ui"),
      h4("Validation issues"), uiOutput("tvp_issues_ui"),
      textInput("tvp_comparisons", "Included comparison IDs (comma-separated; blank means all)", value = ""),
      textInput("tvp_analysis_id", "Analysis ID", value = "treated_vs_parent_raw_cpm"),
      textInput("tvp_input_mode", "Scientific input mode", value = "Raw per-sample SummaryTable"),
      helpText("CPM normalization is performed inside the treated-versus-parent R analysis."),
      checkboxInput("tvp_annotate", "Annotate results with classifier predictions", FALSE),
      conditionalPanel("input.tvp_annotate", textInput("tvp_classifier_target", "Classifier target", value = "recommended")),
      checkboxInput("tvp_dry_run", "Dry run", TRUE), checkboxInput("tvp_force", "Force rerun", FALSE),
      checkboxInput("tvp_keep_going", "Keep going", FALSE), actionButton("tvp_run", "Run treated-versus-parent analysis"),
      helpText("Treated-versus-parent analysis uses raw per-sample SummaryTable outputs and performs CPM normalization inside the R analysis."),
      helpText("MidLC normalization targets used by the essentiality classifier are not used for treated-versus-parent analysis."),
      helpText("Results from the 10,000-read smoke subset validate software behavior and should not be interpreted as a final H2O2 fitness screen."),
      h4("Analysis status"), uiOutput("tvp_status_ui"), h4("Comparison summary"), uiOutput("tvp_summary_ui"),
      h4("Fitness results"), uiOutput("tvp_results_ui")
    )
  )
)

server <- function(input, output, session) {
  log_text <- reactiveVal("Enter project settings, then initialize. Alignment is not run in Step 1.")
  output$command_output <- renderText(log_text())
  output$reference_status <- renderText({
    req(nzchar(input$species))
    status <- reference_status(input$species)
    sprintf("Selected species: %s\nReference runnable: %s\nBowtie2 index exists: %s",
            input$species, if (status$runnable) "yes" else "no",
            if (status$index) "yes" else "no")
  })

  normalization_export_dir <- reactive({
    config_path <- input$normalization_project_config
    project_id <- sub("/config/project\\.yaml$", "", gsub("\\\\", "/", config_path))
    project_id <- basename(project_id)
    file.path(repo_root, "output", "exports", project_id, "qc", "sample_normalization")
  })

  summary_normalized_export_dir <- reactive({
    config_path <- input$summary_normalized_project_config
    project_id <- sub("/config/project\\.yaml$", "", gsub("\\\\", "/", config_path))
    project_id <- basename(project_id)
    file.path(repo_root, "output", "exports", project_id, "qc", "sample_normalization")
  })

  combined_hits_status_path <- reactive({
    config_path <- input$combine_hits_project_config
    project_id <- sub("/config/project\\.yaml$", "", gsub("\\\\", "/", config_path))
    project_id <- basename(project_id)
    file.path(repo_root, "output", "projects", project_id, "manifests", "combined_hits", "combined_hits_status.csv")
  })

  summary_combined_status_path <- reactive({
    config_path <- input$summary_combined_project_config
    project_id <- sub("/config/project\\.yaml$", "", gsub("\\\\", "/", config_path))
    project_id <- basename(project_id)
    file.path(repo_root, "output", "projects", project_id, "manifests", "summary_combined", "summary_combined_status.csv")
  })

  classifier_project_dir <- reactive({
    config_path <- input$classifier_project_config
    project_id <- sub("/config/project\\.yaml$", "", gsub("\\\\", "/", config_path))
    file.path(repo_root, "output", "projects", basename(project_id))
  })

  classifier_status_path <- reactive({
    file.path(classifier_project_dir(), "manifests", "classifier", "classifier_status.csv")
  })

  tvp_project_dir <- reactive({
    project_id <- basename(sub("/config/project\\.yaml$", "", gsub("\\\\", "/", input$tvp_project_config)))
    file.path(repo_root, "output", "projects", project_id)
  })
  tvp_analysis_dir <- reactive(file.path(tvp_project_dir(), "treated_vs_parent", trimws(input$tvp_analysis_id)))
  tvp_table_ui <- function(path, id) {
    if (!file.exists(path)) return(helpText("No output available yet."))
    if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput(paste0(id, "_dt")) else tableOutput(paste0(id, "_table"))
  }
  output$tvp_design_ui <- renderUI(tvp_table_ui(file.path(tvp_project_dir(), "config", "comparison_design.csv"), "tvp_design"))
  output$tvp_issues_ui <- renderUI(tvp_table_ui(file.path(tvp_project_dir(), "config", "comparison_design_issues.csv"), "tvp_issues"))
  output$tvp_status_ui <- renderUI(tvp_table_ui(file.path(tvp_project_dir(), "manifests", "treated_vs_parent", "treated_vs_parent_status.csv"), "tvp_status"))
  output$tvp_summary_ui <- renderUI(tvp_table_ui(file.path(tvp_analysis_dir(), "treated_vs_parent_comparison_summary.csv"), "tvp_summary"))
  output$tvp_results_ui <- renderUI(tvp_table_ui(file.path(tvp_analysis_dir(), "treated_vs_parent_results.csv"), "tvp_results"))
  tvp_paths <- list(
    tvp_design=function() file.path(tvp_project_dir(),"config","comparison_design.csv"),
    tvp_issues=function() file.path(tvp_project_dir(),"config","comparison_design_issues.csv"),
    tvp_status=function() file.path(tvp_project_dir(),"manifests","treated_vs_parent","treated_vs_parent_status.csv"),
    tvp_summary=function() file.path(tvp_analysis_dir(),"treated_vs_parent_comparison_summary.csv"),
    tvp_results=function() file.path(tvp_analysis_dir(),"treated_vs_parent_results.csv"))
  for (id in names(tvp_paths)) local({
    current <- id; path_fun <- tvp_paths[[id]]
    output[[paste0(current,"_table")]] <- renderTable({ req(file.exists(path_fun())); read.csv(path_fun(),check.names=FALSE) })
    if (requireNamespace("DT",quietly=TRUE)) output[[paste0(current,"_dt")]] <- DT::renderDT({
      req(file.exists(path_fun())); DT::datatable(read.csv(path_fun(),check.names=FALSE),filter="top",options=list(pageLength=15,scrollX=TRUE)) })
  })

  observe({
    base <- file.path(classifier_project_dir(), "summary_combined")
    targets <- if (dir.exists(base)) basename(list.dirs(base, recursive = FALSE, full.names = TRUE)) else character(0)
    selected <- if (length(targets)) sort(targets)[[1]] else character(0)
    updateSelectInput(session, "classifier_available_target", choices = sort(targets), selected = selected)
  })

  classifier_display_tag <- reactive({
    mode <- input$classifier_target_mode
    if (identical(mode, "available")) return(input$classifier_available_target)
    if (identical(mode, "manual")) return(trimws(input$classifier_manual_target))
    recommendation <- if (identical(mode, "site-recommended")) "normalization_recommendation.json" else "normalization_feature_recommendation.json"
    path <- file.path(classifier_project_dir(), "sample_normalization", recommendation)
    if (file.exists(path) && requireNamespace("jsonlite", quietly = TRUE)) {
      return(as.character(jsonlite::fromJSON(path)$recommended_target_tag))
    }
    input$classifier_available_target
  })

  output$normalization_recommendation_ui <- renderUI({
    path <- file.path(normalization_export_dir(), "normalization_recommendation.csv")
    if (!file.exists(path)) return(helpText("No recommendation file yet."))
    if (requireNamespace("DT", quietly = TRUE)) {
      DT::DTOutput("normalization_recommendation_dt")
    } else {
      tableOutput("normalization_recommendation_table")
    }
  })
  output$normalization_recommendation_table <- renderTable({
    path <- file.path(normalization_export_dir(), "normalization_recommendation.csv")
    req(file.exists(path))
    read.csv(path, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)
  if (requireNamespace("DT", quietly = TRUE)) {
    output$normalization_recommendation_dt <- DT::renderDT({
      path <- file.path(normalization_export_dir(), "normalization_recommendation.csv")
      req(file.exists(path))
      DT::datatable(read.csv(path, check.names = FALSE), options = list(pageLength = 5))
    })
  }

  output$summary_combined_status_ui <- renderUI({
    path <- summary_combined_status_path()
    if (!file.exists(path)) return(helpText("No combined-summary status file yet."))
    if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("summary_combined_status_dt") else tableOutput("summary_combined_status_table")
  })
  output$summary_combined_status_table <- renderTable({
    path <- summary_combined_status_path(); req(file.exists(path)); read.csv(path, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)
  if (requireNamespace("DT", quietly = TRUE)) {
    output$summary_combined_status_dt <- DT::renderDT({
      path <- summary_combined_status_path(); req(file.exists(path))
      DT::datatable(read.csv(path, check.names = FALSE), options = list(pageLength = 10, scrollX = TRUE))
    })
  }

  output$classifier_status_ui <- renderUI({
    path <- classifier_status_path()
    if (!file.exists(path)) return(helpText("No classifier status file yet."))
    if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("classifier_status_dt") else tableOutput("classifier_status_table")
  })
  output$classifier_status_table <- renderTable({
    path <- classifier_status_path(); req(file.exists(path)); read.csv(path, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)
  if (requireNamespace("DT", quietly = TRUE)) {
    output$classifier_status_dt <- DT::renderDT({
      path <- classifier_status_path(); req(file.exists(path))
      DT::datatable(read.csv(path, check.names = FALSE), filter = "top", options = list(pageLength = 10, scrollX = TRUE))
    })
  }

  output$classifier_summary_ui <- renderUI({
    tag <- classifier_display_tag(); req(nzchar(tag))
    path <- file.path(classifier_project_dir(), "classifier", tag, sprintf("classifier_summary.%s.csv", tag))
    if (!file.exists(path)) return(helpText("No classifier summary for the selected target yet."))
    if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("classifier_summary_dt") else tableOutput("classifier_summary_table")
  })
  output$classifier_summary_table <- renderTable({
    tag <- classifier_display_tag(); path <- file.path(classifier_project_dir(), "classifier", tag, sprintf("classifier_summary.%s.csv", tag))
    req(file.exists(path)); read.csv(path, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)
  if (requireNamespace("DT", quietly = TRUE)) {
    output$classifier_summary_dt <- DT::renderDT({
      tag <- classifier_display_tag(); path <- file.path(classifier_project_dir(), "classifier", tag, sprintf("classifier_summary.%s.csv", tag))
      req(file.exists(path)); DT::datatable(read.csv(path, check.names = FALSE), options = list(scrollX = TRUE))
    })
  }

  output$classifier_predictions_ui <- renderUI({
    tag <- classifier_display_tag(); req(nzchar(tag))
    path <- file.path(classifier_project_dir(), "classifier", tag, sprintf("essentiality_predictions.%s.csv", tag))
    if (!file.exists(path)) return(helpText("No essentiality predictions for the selected target yet."))
    if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("classifier_predictions_dt") else tableOutput("classifier_predictions_table")
  })
  output$classifier_predictions_table <- renderTable({
    tag <- classifier_display_tag(); path <- file.path(classifier_project_dir(), "classifier", tag, sprintf("essentiality_predictions.%s.csv", tag))
    req(file.exists(path)); head(read.csv(path, check.names = FALSE), 25)
  }, striped = TRUE, bordered = TRUE)
  if (requireNamespace("DT", quietly = TRUE)) {
    output$classifier_predictions_dt <- DT::renderDT({
      tag <- classifier_display_tag(); path <- file.path(classifier_project_dir(), "classifier", tag, sprintf("essentiality_predictions.%s.csv", tag))
      req(file.exists(path)); DT::datatable(read.csv(path, check.names = FALSE), filter = "top", options = list(pageLength = 25, scrollX = TRUE))
    })
  }

  output$combined_hits_status_ui <- renderUI({
    path <- combined_hits_status_path()
    if (!file.exists(path)) return(helpText("No combined-hits status file yet."))
    if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("combined_hits_status_dt") else tableOutput("combined_hits_status_table")
  })
  output$combined_hits_status_table <- renderTable({
    path <- combined_hits_status_path(); req(file.exists(path)); read.csv(path, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)
  if (requireNamespace("DT", quietly = TRUE)) {
    output$combined_hits_status_dt <- DT::renderDT({
      path <- combined_hits_status_path(); req(file.exists(path))
      DT::datatable(read.csv(path, check.names = FALSE), options = list(pageLength = 10, scrollX = TRUE))
    })
  }

  output$normalization_target_evaluation_ui <- renderUI({
    path <- file.path(summary_normalized_export_dir(), "normalization_target_evaluation.csv")
    if (!file.exists(path)) return(helpText("No target evaluation file yet."))
    if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("normalization_target_evaluation_dt") else tableOutput("normalization_target_evaluation_table")
  })
  output$normalization_target_evaluation_table <- renderTable({
    path <- file.path(summary_normalized_export_dir(), "normalization_target_evaluation.csv")
    req(file.exists(path)); read.csv(path, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)
  if (requireNamespace("DT", quietly = TRUE)) {
    output$normalization_target_evaluation_dt <- DT::renderDT({
      path <- file.path(summary_normalized_export_dir(), "normalization_target_evaluation.csv")
      req(file.exists(path)); DT::datatable(read.csv(path, check.names = FALSE), options = list(pageLength = 10, scrollX = TRUE))
    })
  }

  output$normalization_target_summary_ui <- renderUI({
    path <- file.path(summary_normalized_export_dir(), "normalization_target_summary.csv")
    if (!file.exists(path)) return(helpText("No target summary file yet."))
    if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("normalization_target_summary_dt") else tableOutput("normalization_target_summary_table")
  })
  output$normalization_target_summary_table <- renderTable({
    path <- file.path(summary_normalized_export_dir(), "normalization_target_summary.csv")
    req(file.exists(path)); read.csv(path, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)
  if (requireNamespace("DT", quietly = TRUE)) {
    output$normalization_target_summary_dt <- DT::renderDT({
      path <- file.path(summary_normalized_export_dir(), "normalization_target_summary.csv")
      req(file.exists(path)); DT::datatable(read.csv(path, check.names = FALSE), options = list(pageLength = 10, scrollX = TRUE))
    })
  }

  output$normalization_feature_recommendation_ui <- renderUI({
    path <- file.path(summary_normalized_export_dir(), "normalization_feature_recommendation.csv")
    if (!file.exists(path)) return(helpText("No feature recommendation file yet."))
    if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput("normalization_feature_recommendation_dt") else tableOutput("normalization_feature_recommendation_table")
  })
  output$normalization_feature_recommendation_table <- renderTable({
    path <- file.path(summary_normalized_export_dir(), "normalization_feature_recommendation.csv")
    req(file.exists(path)); read.csv(path, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)
  if (requireNamespace("DT", quietly = TRUE)) {
    output$normalization_feature_recommendation_dt <- DT::renderDT({
      path <- file.path(summary_normalized_export_dir(), "normalization_feature_recommendation.csv")
      req(file.exists(path)); DT::datatable(read.csv(path, check.names = FALSE), options = list(pageLength = 5))
    })
  }

  output$normalization_comparison_ui <- renderUI({
    path <- file.path(normalization_export_dir(), "normalization_comparison.csv")
    if (!file.exists(path)) return(helpText("No comparison file yet."))
    if (requireNamespace("DT", quietly = TRUE)) {
      DT::DTOutput("normalization_comparison_dt")
    } else {
      tableOutput("normalization_comparison_table")
    }
  })
  output$normalization_comparison_table <- renderTable({
    path <- file.path(normalization_export_dir(), "normalization_comparison.csv")
    req(file.exists(path))
    read.csv(path, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)
  if (requireNamespace("DT", quietly = TRUE)) {
    output$normalization_comparison_dt <- DT::renderDT({
      path <- file.path(normalization_export_dir(), "normalization_comparison.csv")
      req(file.exists(path))
      DT::datatable(read.csv(path, check.names = FALSE), options = list(pageLength = 10, scrollX = TRUE))
    })
  }

  observeEvent(input$prepare_reference, {
    req(nzchar(input$species))
    python <- Sys.which("python")
    if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) {
      log_text("Python was not found on PATH.")
      return()
    }
    script <- file.path(repo_root, "scripts", "local", "ytab_prepare_reference.py")
    args <- c(shQuote(script), "--species", shQuote(input$species),
              "--threads", as.character(input$threads), "--repo-root", shQuote(repo_root))
    result <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) conditionMessage(w),
      error = function(e) paste("Reference preparation failed:", conditionMessage(e))
    )
    status <- attr(result, "status")
    if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$initialize, {
    req(nzchar(input$project_id), nzchar(input$fastq_dir), nzchar(input$species))
    python <- Sys.which("python")
    if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) {
      log_text("Python was not found on PATH.")
      return()
    }
    script <- file.path(repo_root, "scripts", "local", "ytab_init_local_project.py")
    args <- c(shQuote(script), "--project-id", shQuote(input$project_id),
              "--fastq-dir", shQuote(input$fastq_dir), "--species", shQuote(input$species),
              "--threads", as.character(input$threads), "--repo-root", shQuote(repo_root))
    result <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) conditionMessage(w),
      error = function(e) paste("Initialization failed:", conditionMessage(e))
    )
    status <- attr(result, "status")
    if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$run_mapfastq, {
    req(nzchar(input$mapping_project_config))
    python <- Sys.which("python")
    if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) {
      log_text("Python was not found on PATH.")
      return()
    }
    config_path <- input$mapping_project_config
    if (!grepl("^(/|[A-Za-z]:)", config_path)) config_path <- file.path(repo_root, config_path)
    script <- file.path(repo_root, "scripts", "local", "ytab_run_mapfastq.py")
    args <- c(shQuote(script), "--project-config", shQuote(config_path),
              "--threads", as.character(input$threads))
    if (nzchar(trimws(input$mapping_samples))) {
      args <- c(args, "--samples", shQuote(trimws(input$mapping_samples)))
    }
    if (isTRUE(input$mapping_dry_run)) args <- c(args, "--dry-run")
    result <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) conditionMessage(w),
      error = function(e) paste("MapFastq failed:", conditionMessage(e))
    )
    status <- attr(result, "status")
    if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$run_create_hit_file, {
    req(nzchar(input$hits_project_config))
    python <- Sys.which("python")
    if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) {
      log_text("Python was not found on PATH.")
      return()
    }
    config_path <- input$hits_project_config
    if (!grepl("^(/|[A-Za-z]:)", config_path)) config_path <- file.path(repo_root, config_path)
    script <- file.path(repo_root, "scripts", "local", "ytab_run_create_hit_file.py")
    args <- c(shQuote(script), "--project-config", shQuote(config_path),
              "--threads", as.character(input$threads))
    if (nzchar(trimws(input$hits_samples))) {
      args <- c(args, "--samples", shQuote(trimws(input$hits_samples)))
    }
    if (isTRUE(input$hits_dry_run)) args <- c(args, "--dry-run")
    result <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) conditionMessage(w),
      error = function(e) paste("CreateHitFile failed:", conditionMessage(e))
    )
    status <- attr(result, "status")
    if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$run_summary_table, {
    req(nzchar(input$summary_project_config))
    python <- Sys.which("python")
    if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) {
      log_text("Python was not found on PATH.")
      return()
    }
    config_path <- input$summary_project_config
    if (!grepl("^(/|[A-Za-z]:)", config_path)) config_path <- file.path(repo_root, config_path)
    script <- file.path(repo_root, "scripts", "local", "ytab_run_summary_table.py")
    args <- c(shQuote(script), "--project-config", shQuote(config_path),
              "--threads", as.character(input$threads))
    if (nzchar(trimws(input$summary_samples))) {
      args <- c(args, "--samples", shQuote(trimws(input$summary_samples)))
    }
    if (isTRUE(input$summary_dry_run)) args <- c(args, "--dry-run")
    result <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) conditionMessage(w),
      error = function(e) paste("SummaryTable failed:", conditionMessage(e))
    )
    status <- attr(result, "status")
    if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$run_library_diagnostics, {
    req(nzchar(input$diagnostics_project_config))
    python <- Sys.which("python")
    if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) {
      log_text("Python was not found on PATH.")
      return()
    }
    config_path <- input$diagnostics_project_config
    if (!grepl("^(/|[A-Za-z]:)", config_path)) config_path <- file.path(repo_root, config_path)
    script <- file.path(repo_root, "scripts", "local", "ytab_run_library_diagnostics.py")
    args <- c(shQuote(script), "--project-config", shQuote(config_path))
    if (nzchar(trimws(input$diagnostics_samples))) {
      args <- c(args, "--samples", shQuote(trimws(input$diagnostics_samples)))
    }
    if (isTRUE(input$diagnostics_dry_run)) args <- c(args, "--dry-run")
    result <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) conditionMessage(w),
      error = function(e) paste("LibraryDiagnostics failed:", conditionMessage(e))
    )
    status <- attr(result, "status")
    if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$run_sample_normalization, {
    req(nzchar(input$normalization_project_config))
    python <- Sys.which("python")
    if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) {
      log_text("Python was not found on PATH.")
      return()
    }
    config_path <- input$normalization_project_config
    if (!grepl("^(/|[A-Za-z]:)", config_path)) config_path <- file.path(repo_root, config_path)
    script <- file.path(repo_root, "scripts", "local", "ytab_run_sample_normalization.py")
    targets <- if (identical(input$normalization_target_mode, "auto")) "auto" else trimws(input$normalization_manual_targets)
    args <- c(shQuote(script), "--project-config", shQuote(config_path),
              "--targets", shQuote(targets), "--sample-mode", input$normalization_sample_mode,
              "--threads", as.character(input$threads))
    if (identical(input$normalization_target_mode, "auto")) {
      args <- c(args, "--min-site-retention", as.character(input$normalization_min_retention),
                "--auto-step", as.character(input$normalization_auto_step))
      if (nzchar(trimws(input$normalization_auto_min))) {
        args <- c(args, "--auto-min-target", shQuote(trimws(input$normalization_auto_min)))
      }
      if (nzchar(trimws(input$normalization_auto_max))) {
        args <- c(args, "--auto-max-target", shQuote(trimws(input$normalization_auto_max)))
      }
    }
    if (isTRUE(input$normalization_dry_run)) args <- c(args, "--dry-run")
    if (isTRUE(input$normalization_force)) args <- c(args, "--force")
    result <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) conditionMessage(w),
      error = function(e) paste("Sample normalization failed:", conditionMessage(e))
    )
    status <- attr(result, "status")
    if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$run_summary_normalized, {
    req(nzchar(input$summary_normalized_project_config))
    python <- Sys.which("python")
    if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) {
      log_text("Python was not found on PATH.")
      return()
    }
    config_path <- input$summary_normalized_project_config
    if (!grepl("^(/|[A-Za-z]:)", config_path)) config_path <- file.path(repo_root, config_path)
    targets <- if (identical(input$summary_normalized_target_mode, "manual")) {
      trimws(input$summary_normalized_manual_targets)
    } else input$summary_normalized_target_mode
    script <- file.path(repo_root, "scripts", "local", "ytab_run_summary_normalized.py")
    args <- c(shQuote(script), "--project-config", shQuote(config_path),
              "--targets", shQuote(targets), "--sample-mode", input$summary_normalized_sample_mode,
              "--min-feature-retention", as.character(input$summary_normalized_min_retention),
              "--threads", as.character(input$threads))
    if (isTRUE(input$summary_normalized_dry_run)) args <- c(args, "--dry-run")
    if (isTRUE(input$summary_normalized_force)) args <- c(args, "--force")
    result <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) conditionMessage(w),
      error = function(e) paste("Normalized SummaryTable failed:", conditionMessage(e))
    )
    status <- attr(result, "status")
    if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$run_combine_hits, {
    req(nzchar(input$combine_hits_project_config))
    python <- Sys.which("python")
    if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) {
      log_text("Python was not found on PATH.")
      return()
    }
    config_path <- input$combine_hits_project_config
    if (!grepl("^(/|[A-Za-z]:)", config_path)) config_path <- file.path(repo_root, config_path)
    target <- if (identical(input$combine_hits_target_mode, "manual")) {
      trimws(input$combine_hits_manual_target)
    } else input$combine_hits_target_mode
    script <- file.path(repo_root, "scripts", "local", "ytab_run_combine_hits.py")
    args <- c(shQuote(script), "--project-config", shQuote(config_path), "--target", shQuote(target))
    if (nzchar(trimws(input$combine_hits_samples))) {
      args <- c(args, "--samples", shQuote(trimws(input$combine_hits_samples)))
    }
    if (isTRUE(input$combine_hits_dry_run)) args <- c(args, "--dry-run")
    if (isTRUE(input$combine_hits_force)) args <- c(args, "--force")
    result <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) conditionMessage(w),
      error = function(e) paste("Combine normalized parent hits failed:", conditionMessage(e))
    )
    status <- attr(result, "status")
    if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$run_summary_combined, {
    req(nzchar(input$summary_combined_project_config))
    python <- Sys.which("python")
    if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) {
      log_text("Python was not found on PATH.")
      return()
    }
    config_path <- input$summary_combined_project_config
    if (!grepl("^(/|[A-Za-z]:)", config_path)) config_path <- file.path(repo_root, config_path)
    target <- if (identical(input$summary_combined_target_mode, "manual")) {
      trimws(input$summary_combined_manual_target)
    } else input$summary_combined_target_mode
    script <- file.path(repo_root, "scripts", "local", "ytab_run_summary_combined.py")
    args <- c(shQuote(script), "--project-config", shQuote(config_path), "--target", shQuote(target),
              "--threads", as.character(input$threads))
    if (isTRUE(input$summary_combined_dry_run)) args <- c(args, "--dry-run")
    if (isTRUE(input$summary_combined_force)) args <- c(args, "--force")
    result <- tryCatch(
      system2(python, args = args, stdout = TRUE, stderr = TRUE),
      warning = function(w) conditionMessage(w),
      error = function(e) paste("Combined parent SummaryTable failed:", conditionMessage(e))
    )
    status <- attr(result, "status")
    if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$run_classifier, {
    req(nzchar(input$classifier_project_config))
    python <- Sys.which("python"); if (!nzchar(python)) python <- Sys.which("python3")
    if (!nzchar(python)) { log_text("Python was not found on PATH."); return() }
    config_path <- input$classifier_project_config
    if (!grepl("^(/|[A-Za-z]:)", config_path)) config_path <- file.path(repo_root, config_path)
    target <- if (identical(input$classifier_target_mode, "manual")) trimws(input$classifier_manual_target) else if (identical(input$classifier_target_mode, "available")) input$classifier_available_target else input$classifier_target_mode
    script <- file.path(repo_root, "scripts", "local", "ytab_run_classifier.py")
    args <- c(shQuote(script), "--project-config", shQuote(config_path), "--target", shQuote(target))
    if (nzchar(trimws(input$classifier_seed))) args <- c(args, "--seed", shQuote(trimws(input$classifier_seed)))
    if (isTRUE(input$classifier_dry_run)) args <- c(args, "--dry-run")
    if (isTRUE(input$classifier_force)) args <- c(args, "--force")
    if (isTRUE(input$classifier_save_final_target)) args <- c(args, "--save-final-target")
    result <- tryCatch(system2(python, args = args, stdout = TRUE, stderr = TRUE),
                       warning = function(w) conditionMessage(w),
                       error = function(e) paste("Essentiality classifier failed:", conditionMessage(e)))
    status <- attr(result, "status"); if (!is.null(status)) result <- c(result, sprintf("Exit status: %d", status))
    log_text(paste(result, collapse = "\n"))
  })

  observeEvent(input$tvp_design, {
    python <- Sys.which("python"); if (!nzchar(python)) python <- Sys.which("python3")
    config_path <- input$tvp_project_config; if (!grepl("^(/|[A-Za-z]:)",config_path)) config_path <- file.path(repo_root,config_path)
    args <- c(shQuote(file.path(repo_root,"scripts","local","ytab_init_comparison_design.py")),"--project-config",shQuote(config_path),"--overwrite","--print-design")
    result <- tryCatch(system2(python,args=args,stdout=TRUE,stderr=TRUE),error=function(e) conditionMessage(e)); log_text(paste(result,collapse="\n"))
  })
  observeEvent(input$tvp_run, {
    python <- Sys.which("python"); if (!nzchar(python)) python <- Sys.which("python3")
    config_path <- input$tvp_project_config; if (!grepl("^(/|[A-Za-z]:)",config_path)) config_path <- file.path(repo_root,config_path)
    args <- c(shQuote(file.path(repo_root,"scripts","local","ytab_run_treated_vs_parent.py")),"--project-config",shQuote(config_path),"--analysis-id",shQuote(input$tvp_analysis_id),"--input-mode","raw-summary")
    if (nzchar(trimws(input$tvp_comparisons))) args <- c(args,"--comparisons",shQuote(trimws(input$tvp_comparisons)))
    if (isTRUE(input$tvp_annotate)) args <- c(args,"--annotate-classifier","--classifier-target",shQuote(input$tvp_classifier_target))
    if (isTRUE(input$tvp_dry_run)) args <- c(args,"--dry-run"); if (isTRUE(input$tvp_force)) args <- c(args,"--force"); if (isTRUE(input$tvp_keep_going)) args <- c(args,"--keep-going")
    result <- tryCatch(system2(python,args=args,stdout=TRUE,stderr=TRUE),error=function(e) conditionMessage(e)); log_text(paste(result,collapse="\n"))
  })
}

shinyApp(ui, server)
