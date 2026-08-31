library(shiny)
library(bslib)

`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
app_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL) %||% file.path(getwd(), "app.R")
app_dir <- normalizePath(dirname(app_file), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(app_dir, "../.."), winslash = "/", mustWork = TRUE)
for (helper in c("project_discovery.R", "project_state.R", "process_helpers.R", "preprocessing_status.R", "sample_selector.R", "job_manager.R", "job_progress.R", "navigation.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R", "table_display_helpers.R", "glabrata_annotation_lookup.R", "ui_landing.R", "ui_preprocessing.R", "qc_result_state.R", "qc_plot_utils.R", "qc_mapping_stats_plot.R", "qc_summary_library_plots.R", "qc_library_diagnostics_plots.R", "ui_qc.R", "fitness_design_state.R", "fitness_result_state.R", "fitness_generated_plots.R", "fitness_condition_control_plot.R", "fitness_plot_utils.R", "ui_fitness.R", "essentiality_targets.R", "essentiality_state.R", "essentiality_results.R", "essentiality_commands.R", "essentiality_generated_plots.R", "ui_essentiality.R", "essentiality_server.R", "gene_domain_explorer_state.R", "ui_gene_domain_explorer.R", "gene_domain_explorer_server.R", "genome_browser_state.R", "ui_genome_browser.R", "comparative_resources.R", "comparative_project_data.R", "ui_comparative.R", "ui_workspace.R"))
  source(file.path(app_dir, "R", helper), local = TRUE)

detected_cpu <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cpu) || detected_cpu < 1L) detected_cpu <- 2L
detected_cpu <- as.integer(detected_cpu)
species_dir <- file.path(repo_root, "resources", "species")
species_choices <- if (dir.exists(species_dir)) sort(basename(list.dirs(species_dir, recursive = FALSE))) else character()
shinyfiles_available <- requireNamespace("shinyFiles", quietly = TRUE)
launcher_config <- read_launcher_config(app_dir)
demo_project_id <- as.character(launcher_config$demo_project_id %||% "")

ui <- tagList(
  tags$head(tags$title("YTAB — Yeast Transposon Analysis Browser"), tags$link(rel="stylesheet", href="ytab.css"),
    tags$link(rel="stylesheet", href="ytab-landing.css"), tags$link(rel="stylesheet", href="ytab-qc.css"), tags$link(rel="stylesheet",href="job-progress.css"), tags$link(rel="stylesheet",href="ytab_release_ui.css?v=20260828_2"),
    tags$script(src="https://cdn.jsdelivr.net/npm/igv@3/dist/igv.min.js"),
    tags$script(src="ytab_igv_browser.js?v=20260828_5")),
  uiOutput("app_shell")
)

server <- function(input, output, session) {
  state <- new_project_state()
  projects <- reactiveVal(discover_ytab_projects(repo_root, include_tests = FALSE))
  new_project_open <- reactiveVal(FALSE)
  import_project_open <- reactiveVal(FALSE)
  fastq_directory <- reactiveVal("")
  fastq_files <- reactiveVal(data.frame())
  fastq_scan_state <- reactiveVal("no_directory")
  fastq_scan_warning <- reactiveVal("")
  initialization_running <- reactiveVal(FALSE)
  log_text <- reactiveVal("No command has run in this session.")
  exports_log <- reactiveVal("")
  jobs <- new_job_manager()
  status_tick <- reactiveVal(0L)
  project_created <- reactiveVal(FALSE)
  completed_job_seen <- reactiveVal("")
  last_project_path <- reactiveVal("")
  navigation_location <- reactiveVal(navigation_state_defaults("overview"))
  go_to <- function(top,nested_id=NULL,nested_value=NULL,results_view=NULL){
    navigation_location(navigate_workspace(
      session, navigation_location(), top, nested_id, nested_value, results_view
    ))
  }
  for (navigation_input in c("workspace_tabs", names(navigation_nested_fields))) local({
    key <- navigation_input
    observeEvent(input[[key]], {
      navigation_location(record_navigation_selection(
        navigation_location(), key, input[[key]]
      ))
    }, ignoreInit = TRUE)
  })

  initial_raw <- Sys.getenv("YTAB_PROJECT_CONFIG", "")
  if (nzchar(initial_raw)) {
    checked <- validate_project_config_path(initial_raw, repo_root)
    if (checked$valid) {
      session$onFlushed(function() updateSelectInput(session, "project_selector", selected = checked$path), once = TRUE)
    } else state$values$landing_warning <- checked$message
  }

  output$app_shell <- renderUI({
    if (identical(state$values$view, "landing")) landing_ui(species_choices, detected_cpu, shinyfiles_available) else workspace_ui(detected_cpu)
  })
  output$new_project_visible <- reactive(new_project_open()); outputOptions(output, "new_project_visible", suspendWhenHidden = FALSE)
  output$import_project_visible <- reactive(import_project_open()); outputOptions(output, "import_project_visible", suspendWhenHidden = FALSE)
  output$landing_warning <- renderUI(if (nzchar(state$values$landing_warning)) tags$div(class="alert alert-warning", state$values$landing_warning))

  project_choices <- reactive({
    p <- projects(); if (!nrow(p)) return(list())
    labels <- vapply(seq_len(nrow(p)), function(i) project_display_label(read_project_summary(p$project_config[[i]], repo_root)), "")
    as.list(setNames(unname(p$project_config), unname(labels)))
  })
  observe({
    choices <- project_choices(); selected <- isolate(input$project_selector %||% "")
    initial_path <- if (nzchar(initial_raw)) ytab_resolve_path(initial_raw, repo_root) else ""
    values<-unname(as.character(unlist(choices,use.names=FALSE)));if (nzchar(initial_path) && initial_path %in% values) selected <- initial_path
    if (!selected %in% values) selected <- character()
    updateSelectizeInput(session, "project_selector", choices = c(list("Select a project…" = ""), choices), selected = selected, server = TRUE)
  })
  refresh_projects <- function() projects(discover_ytab_projects(repo_root, include_tests = isTRUE(input$show_test_projects)))
  observeEvent(input$refresh_projects, { refresh_projects(); showNotification("Project list refreshed.") })
  observeEvent(input$show_test_projects, refresh_projects(), ignoreInit = TRUE)
  observeEvent(input$show_new_project, new_project_open(!new_project_open()))
  observeEvent(input$show_import_project, import_project_open(!import_project_open()))

  selected_info <- reactive({ req(input$project_selector); tryCatch(read_project_summary(input$project_selector, repo_root), error = function(e) NULL) })
  output$open_project_action <- renderUI(actionButton("open_project", "Open project", class="btn-primary ytab-button",
    disabled = if (is.null(selected_info())) "disabled" else NULL))
  output$selected_project_preview <- renderUI({
    info <- selected_info(); req(info)
    tags$div(class="ytab-project-preview",
      tags$span(class="ytab-preview-badge", info$species),
      tags$span(class="ytab-preview-badge", sprintf("%d samples", info$included_count)),
      tags$span(class="ytab-preview-badge", if(info$analysis_ready) "analysis ready" else "preprocessing incomplete"),
      if (info$project_type %in% c("test", "temporary")) tags$span(class="ytab-preview-badge ytab-test-badge", "Test"))
  })
  demo_project <- reactive(select_demo_project(projects(), demo_project_id))
  output$demo_project_card <- renderUI({
    demo <- demo_project(); if (is.null(demo)) return(NULL)
    panel_card("Workflow demonstration project", tags$p(project_display_label(as.list(demo[1, ]))),
      tags$p(class="ytab-warning", "Workflow demonstration only; not final biological results."),
      actionButton("open_demo", "Open demonstration project", class="btn-primary ytab-button"))
  })
  open_path <- function(path) {
    checked <- validate_project_config_path(path, repo_root)
    if (!checked$valid) { showNotification(checked$message, type="error", duration=NULL); return(FALSE) }
    state$load(checked$path, repo_root); state$enter(); status_tick(status_tick()+1L);last_project_path(checked$path);jobs$recover_project(checked$path)
    workspace_tab <- isolate(state$values$workspace_tab)
    location <- navigation_state_defaults(workspace_tab)
    navigation_location(location)
    session$onFlushed(function() restore_navigation_state(session, location), once=TRUE)
    TRUE
  }
  observeEvent(input$open_project, { info <- selected_info(); if (!is.null(info)) open_path(info$project_config) })
  observeEvent(input$open_demo, { demo <- demo_project(); if (!is.null(demo)) open_path(demo$project_config[[1]]) })
  switch_to_landing<-function(){last<-last_project_path();state$clear();jobs$clear_display_state();navigation_location(navigation_state_defaults("overview"));projects(discover_ytab_projects(repo_root,include_tests=FALSE));if(nzchar(last))last_project_path(last);session$onFlushed(function()updateSelectizeInput(session,"project_selector",selected=character()),once=TRUE)}
  observeEvent(input$switch_project,{
    if(jobs$job_is_running()){
      showModal(modalDialog(title="Pipeline job running","A pipeline job is running. Switching projects will not stop the job. You may continue monitoring it from its project later.",footer=tagList(modalButton("Stay here"),actionButton("confirm_switch_project","Switch project",class="btn-danger"))))
    } else switch_to_landing()
  },ignoreInit=TRUE)
  observeEvent(input$confirm_switch_project,{removeModal();switch_to_landing()},ignoreInit=TRUE)
  observeEvent(input$return_from_new_project,new_project_open(FALSE),ignoreInit=TRUE)
  observeEvent(input$return_from_import_project,import_project_open(FALSE),ignoreInit=TRUE)

  output$documentation_guidance <- renderUI({
    req((input$show_launch_command %||% 0) + (input$view_setup_guide %||% 0) > 0)
    if ((input$view_setup_guide %||% 0) > (input$show_launch_command %||% 0))
      tags$p("See ", tags$code("docs/setup_local.md"), " for local setup and preprocessing guidance.")
    else tags$pre("mamba activate ytab-local\n./scripts/local/ytab_launch_app.sh")
  })

  if (shinyfiles_available) {
    roots <- c(Home = normalizePath("~", mustWork=TRUE), `Repository parent` = dirname(repo_root))
    if (.Platform$OS.type == "windows") roots <- shinyFiles::getVolumes() else {
      if (dir.exists("/Volumes")) roots <- c(roots, Volumes="/Volumes")
      if (dir.exists("/mnt")) roots <- c(roots, Mounts="/mnt")
    }
    shinyFiles::shinyDirChoose(input, "fastq_dir_choose", roots=roots, session=session, restrictions=system.file(package="base"))
    observeEvent(input$fastq_dir_choose, {
      parsed <- shinyFiles::parseDirPath(roots, input$fastq_dir_choose)
      if (length(parsed) && nzchar(parsed[[1]]) && dir.exists(parsed[[1]])) {
        path <- normalize_local_directory(parsed[[1]])
        fastq_directory(path); fastq_files(data.frame()); fastq_scan_state("directory_selected"); fastq_scan_warning("")
        updateTextInput(session, "fastq_dir_path", value=path)
      }
    }, ignoreInit=TRUE)
  }
  observeEvent(input$use_fastq_path, {
    path <- tryCatch(normalize_local_directory(input$fastq_dir_path), error=function(e)e)
    if(inherits(path,"error")){showNotification(conditionMessage(path),type="error");return()}
    fastq_directory(path); fastq_files(data.frame()); fastq_scan_state("directory_selected"); fastq_scan_warning("")
  }, ignoreInit=TRUE)
  observeEvent(input$scan_fastqs, {
    if(!nzchar(fastq_directory()) || !dir.exists(fastq_directory())){showNotification("Select a valid FASTQ directory first.",type="error");return()}
    fastq_scan_state("scanning")
    files <- withProgress(message="Scanning top-level FASTQ files", value=0.2, {
      result <- scan_fastq_directory(fastq_directory(), limit=5000L)
      incProgress(0.8); result
    })
    truncated <- isTRUE(attr(files,"truncated")); fastq_files(files); fastq_scan_state("scan_complete")
    fastq_scan_warning(if(truncated)"More than 5,000 FASTQ files were found; the preview is limited to the first 5,000 files." else "")
  }, ignoreInit=TRUE)
  output$fastq_dir_text <- renderText(if(nzchar(fastq_directory())) fastq_directory() else "No directory selected")
  output$fastq_preview <- DT::renderDT(DT::datatable(fastq_files(), rownames=FALSE, options=list(pageLength=8,scrollX=TRUE)))
  output$fastq_message <- renderUI({
    if(nzchar(fastq_scan_warning())) return(tags$div(class="alert alert-warning",fastq_scan_warning()))
    if(identical(fastq_scan_state(),"scan_complete") && !nrow(fastq_files())) tags$div(class="alert alert-warning","No FASTQ files were found in the selected directory.")
    else tags$p(class="text-muted",switch(fastq_scan_state(),no_directory="Select a directory, then scan FASTQs.",directory_selected="Directory selected. Click Scan FASTQs to review files.",scanning="Scanning…",scan_complete=sprintf("Scan complete: %d FASTQ files.",nrow(fastq_files())),initializing="Initializing project…",failed="Initialization failed.",""))
  })
  validate_project_id <- function(id) {
    id<-trimws(id %||% ""); if(!nzchar(id)) return("Project ID is required.")
    if(!grepl("^[A-Za-z0-9_-]+$",id) || grepl("\\.\\.|[/\\\\]",id)) return("Project ID may contain only letters, numbers, hyphens, and underscores.")
    if(dir.exists(file.path(repo_root,"output","projects",id))) return("This project already exists. Open it from the existing-project list.")
    ""
  }
  validate_import_project_id <- function(id) {
    id<-trimws(id %||% ""); if(!nzchar(id)) return("Project ID is required.")
    if(!grepl("^[A-Za-z0-9_-]+$",id) || grepl("\\.\\.|[/\\\\]",id)) return("Project ID may contain only letters, numbers, hyphens, and underscores.")
    ""
  }
  init_args <- reactive(c(file.path(repo_root,"scripts/local/ytab_init_local_project.py"),"--project-id",trimws(input$project_id %||% ""),"--fastq-dir",fastq_directory(),"--species",input$species %||% "","--threads",as.character(input$threads %||% 2L),"--repo-root",repo_root))
  import_args <- reactive({
    args <- c(file.path(repo_root,"scripts/local/ytab_import_existing_hit_project.py"),
      "--project-id",trimws(input$import_project_id %||% ""),"--species",input$import_species %||% "",
      "--hit-dir",input$import_hit_dir_path %||% "","--threads",as.character(input$import_threads %||% 2L),
      "--repo-root",repo_root)
    metadata <- trimws(input$import_metadata_path %||% "")
    if(nzchar(metadata)) args <- c(args,"--sample-metadata",metadata)
    args
  })
  output$init_command <- renderText(tryCatch(format_command_for_display(locate_python_executable(),init_args()),error=function(e)conditionMessage(e)))
  initialize_ready <- reactive(nzchar(fastq_directory()) && identical(fastq_scan_state(),"scan_complete") && nrow(fastq_files())>0L && !nzchar(validate_project_id(input$project_id)) && nzchar(input$species %||% "") && as.integer(input$threads %||% 0L)>=2L && !initialization_running())
  output$initialize_button <- renderUI(actionButton("initialize", "Initialize project", class="btn-primary", disabled=if(initialize_ready())NULL else "disabled"))
  observeEvent(input$initialize, {
    if(initialization_running())return()
    error <- validate_project_id(input$project_id); if(nzchar(error)){showNotification(error,type="error");return()}
    if(!identical(fastq_scan_state(),"scan_complete") || !nrow(fastq_files())){showNotification("Scan the FASTQ directory and confirm at least one file before initialization.",type="error");return()}
    if(!nzchar(input$species %||% "")){showNotification("Choose a species/reference.",type="error");return()}
    if(as.integer(input$threads %||% 0L)<2L){showNotification("Use at least two threads.",type="error");return()}
    script<-init_args()[[1]];if(!file.exists(script)||dir.exists(script)){showNotification("Initialization script is unavailable.",type="error");return()}
    initialization_running(TRUE);fastq_scan_state("initializing");on.exit(initialization_running(FALSE),add=TRUE)
    result<-tryCatch(run_process_sync(locate_python_executable(),init_args(),wd=repo_root),error=function(e)list(status=1L,stdout="",stderr=conditionMessage(e)))
    log_text(paste(c(format_command_for_display(locate_python_executable(),init_args()),result$stdout,result$stderr,paste("Exit status:",result$status)),collapse="\n")); status<-result$status
    generated<-file.path(repo_root,"output","projects",trimws(input$project_id),"config","project.yaml")
    if(identical(as.integer(status),0L) && validate_project_config_path(generated,repo_root)$valid){
      status_args<-c(file.path(repo_root,"scripts/local/ytab_project_status.py"),"--project-config",generated,"--show-next")
      status_result<-run_process_sync(locate_python_executable(),status_args,wd=repo_root)
      log_text(paste(c(log_text(),status_result$stdout,status_result$stderr),collapse="\n"))
      projects(discover_ytab_projects(repo_root));project_created(TRUE);open_path(generated)
      go_to("preprocessing","preprocessing_tabs","samples")
      showNotification("Project created. Review samples and run mapping; mapping has not started.",type="message")
    } else {fastq_scan_state("failed");showNotification("Project initialization failed; see command output.",type="error",duration=NULL)}
  }, ignoreInit=TRUE)
  observeEvent(input$import_existing_hit_project, {
    error <- validate_import_project_id(input$import_project_id); if(nzchar(error)){showNotification(error,type="error");return()}
    if(!nzchar(trimws(input$import_hit_dir_path %||% ""))){showNotification("Provide a CreateHitFile output directory.",type="error");return()}
    if(!nzchar(input$import_species %||% "")){showNotification("Choose a species/reference.",type="error");return()}
    script<-import_args()[[1]];if(!file.exists(script)||dir.exists(script)){showNotification("Hit-file import script is unavailable.",type="error");return()}
    result<-tryCatch(run_process_sync(locate_python_executable(),import_args(),wd=repo_root),error=function(e)list(status=1L,stdout="",stderr=conditionMessage(e)))
    log_text(paste(c(format_command_for_display(locate_python_executable(),import_args()),result$stdout,result$stderr,paste("Exit status:",result$status)),collapse="\n"))
    generated<-file.path(repo_root,"output","projects",trimws(input$import_project_id),"config","project.yaml")
    if(identical(as.integer(result$status),0L) && validate_project_config_path(generated,repo_root)$valid){
      status_args<-c(file.path(repo_root,"scripts/local/ytab_project_status.py"),"--project-config",generated,"--show-next")
      status_result<-run_process_sync(locate_python_executable(),status_args,wd=repo_root)
      log_text(paste(c(log_text(),status_result$stdout,status_result$stderr),collapse="\n"))
      projects(discover_ytab_projects(repo_root));project_created(TRUE);open_path(generated)
      go_to("preprocessing","preprocessing_tabs","summary_tables")
      showNotification("Hit-file project imported. SummaryTable is the next downstream step.",type="message")
    } else showNotification("Hit-file import failed; see command output.",type="error",duration=NULL)
  }, ignoreInit=TRUE)

  active_project_path <- reactive({ req(state$values$path); checked<-validate_project_config_path(state$values$path,repo_root); validate(need(checked$valid,checked$message)); checked$path })
  active <- reactive({ status_tick(); req(state$values$project); state$values$project })
  observeEvent(status_tick(), {
    if (identical(state$values$view, "workspace")) {
      location <- navigation_location()
      session$onFlushed(function() restore_navigation_state(session, location),
                        once = TRUE)
    }
  }, ignoreInit = TRUE)
  observeEvent(active_project_path(),{
    for(route in c("ytab-diagnostics","ytab-diagnostics-project","ytab-diagnostics-export","ytab-project-output","ytab-reference","ytab-igv-aliases","ytab-igv-assets"))
      suppressWarnings(try(removeResourcePath(route),silent=TRUE))
    roots<-diagnostic_resource_roots(active())
    for(route in names(roots))addResourcePath(route,roots[[route]])
    addResourcePath("ytab-project-output",active()$project_root)
    reference_dir<-genome_browser_reference_dir(active())
    if(nzchar(reference_dir)&&dir.exists(reference_dir))addResourcePath("ytab-reference",reference_dir)
    alias_dir<-genome_browser_alias_dir(active())
    if(nzchar(alias_dir)&&dir.exists(alias_dir))addResourcePath("ytab-igv-aliases",alias_dir)
    asset_dir<-genome_browser_asset_dir(active())
    if(nzchar(asset_dir)&&dir.exists(asset_dir))addResourcePath("ytab-igv-assets",asset_dir)
  },ignoreInit=FALSE)
  sample_pipeline_status <- reactive({status_tick();build_sample_pipeline_status(active())})
  preprocessing_selected_samples <- reactiveVal(character())
  selection_project <- reactiveVal("")
  qc_selected_samples <- reactiveVal(character())
  qc_selection_project <- reactiveVal("")
  diagnostic_inventory_state <- reactiveVal(data.frame())
  diagnostic_selected_run <- reactiveVal("All")
  genome_browser_custom_tracks <- reactiveVal(character())
  fitness_selected_comparisons <- reactiveVal(character())
  fitness_selection_project <- reactiveVal("")
  fitness_design_tick <- reactiveVal(0L)
  fitness_selected_result <- reactiveVal("")
  observeEvent(active_project_path(),{
    path<-active_project_path();if(!identical(path,selection_project())){d<-active()$samples;inc<-if("include"%in%names(d))tolower(as.character(d$include))%in%c("true","1","yes")else rep(TRUE,nrow(d));preprocessing_selected_samples(as.character(d$sample[inc]));selection_project(path)}
  },ignoreInit=FALSE)
  observeEvent(active_project_path(),{
    path<-active_project_path();if(!identical(path,qc_selection_project())){eligible<-qc_sample_eligibility(active(),sample_pipeline_status());qc_selected_samples(eligible$Sample[eligible$Eligible=="Yes"]);qc_selection_project(path);diagnostic_inventory_state(build_diagnostic_file_inventory(active()));diagnostic_selected_run("All")}
  },ignoreInit=FALSE)
  sample_selector_server("preprocessing_samples",sample_data=reactive(active()$samples),status_data=sample_pipeline_status,selected_state=preprocessing_selected_samples)
  job_progress_server("global_job",jobs,active,compact=TRUE)
  job_progress_server("mapping_job",jobs,active,compact=FALSE)
  job_progress_server("library_diagnostics_job",jobs,active,compact=FALSE)
  job_progress_server("fitness_job",jobs,active,compact=FALSE)
  included_samples <- reactive({d<-active()$samples;inc<-if("include"%in%names(d))tolower(as.character(d$include))%in%c("true","1","yes")else rep(TRUE,nrow(d));as.character(d$sample[inc])})
  output$active_project_badge <- renderUI({p<-active();tagList(tags$b(p$project_id),tags$span(p$species),tags$span(if(p$analysis_ready)"analysis ready" else "preprocessing"))})
  overview_stage_rows <- function(project, status) {
    next_stage <- next_preprocessing_stage(project, repo_root, status)
    core_complete <- identical(next_stage, "Complete")
    has_files <- function(path, pattern = NULL) {
      if (!dir.exists(path)) return(FALSE)
      files <- list.files(path, pattern = pattern, recursive = TRUE, full.names = TRUE)
      any(file.exists(files) & !dir.exists(files) & file.info(files)$size > 0)
    }
    library_ready <- tryCatch(nrow(build_diagnostic_file_inventory(project)) > 0L,
                              error = function(e) has_files(file.path(project$project_root, "library_diagnostics")))
    classifier_ready <- tryCatch(length(discover_classifier_results(project$project_root)) > 0L,
                                 error = function(e) has_files(file.path(project$project_root, "classifier"),
                                                               "essentiality_predictions.*\\.(csv|tsv)$"))
    fitness_ready <- tryCatch(length(discover_fitness_results(project)) > 0L,
                              error = function(e) has_files(file.path(project$project_root, "treated_vs_parent"),
                                                            "^treated_vs_parent_results\\.csv$"))
    exports_ready <- has_files(project$export_root)
    data.frame(
      Stage = c("Core preprocessing", "Library diagnostics", "Essentiality classifier",
                "Fitness analysis", "Reports / exports"),
      Status = c(if (core_complete) "Complete" else paste("Next:", next_stage),
                 if (library_ready) "Available" else "Not run",
                 if (classifier_ready) "Available" else "Not run",
                 if (fitness_ready) "Available" else "Not run",
                 if (exports_ready) "Available" else "Not run"),
      Complete = c(core_complete, library_ready, classifier_ready, fitness_ready, exports_ready),
      stringsAsFactors = FALSE
    )
  }
  output$overview_ui <- renderUI({
    p<-active();s<-sample_pipeline_status();ref<-reference_readiness(p,repo_root)
    next_stage<-next_preprocessing_stage(p,repo_root,s)
    preprocessing_complete<-identical(next_stage,"Complete")
    stages<-overview_stage_rows(p,s)
    progress<-round(100*sum(stages$Complete)/nrow(stages))
    analysis_label<-if(preprocessing_complete||isTRUE(p$analysis_ready))"Ready"else"Not ready"
    next_label<-if(preprocessing_complete)"Complete"else next_stage
    stage_cards<-lapply(seq_len(nrow(stages)),function(i)tags$div(tags$b(stages$Status[[i]]),stages$Stage[[i]]))
    tagList(
      tags$div(class="ytab-stat-grid",
        tags$div(tags$b(p$project_id)," Project"),
        tags$div(tags$b(p$included_count)," Included samples"),
        tags$div(tags$b(ref$label)," Reference"),
        tags$div(tags$b(analysis_label)," Analysis")),
      tags$div(class="progress",tags$div(class="progress-bar",style=sprintf("width:%d%%",progress),paste0(progress,"%"))),
      tags$p(tags$b(if(preprocessing_complete)"Status: "else"Next action: "),next_label),
      tags$div(class="ytab-stat-grid",stage_cards)
    )
  })
  output$overview_actions <- renderUI({
    stage<-next_preprocessing_stage(active(),repo_root,sample_pipeline_status())
    if(identical(stage,"Complete")) {
      tags$div(class="ytab-actions",
        actionButton("overview_enter_analysis","Enter analysis browser",class="btn-primary"),
        actionButton("refresh_project_status","Refresh status"))
    } else {
      tags$div(class="ytab-actions",
        actionButton("overview_continue","Continue preprocessing",class="btn-primary"),
        actionButton("refresh_project_status","Refresh status"))
    }
  })
  output$project_metadata_ui <- renderUI({p<-active();tags$dl(class="ytab-meta",tags$dt("Project ID"),tags$dd(p$project_id),tags$dt("Project config"),tags$dd(tags$code(p$project_config)),tags$dt("FASTQ directory"),tags$dd(tags$code(p$fastq_directory)),tags$dt("Sample sheet"),tags$dd(tags$code(p$sample_sheet)),tags$dt("Threads"),tags$dd(p$threads),tags$dt("Reference"),tags$dd(p$reference$fasta %||% "Not resolved"))})
  output$sample_metadata_summary<-renderUI({p<-active();tags$p(sprintf("%d included samples; %d parent; %d treated. Temporary stage selection does not edit project inclusion.",p$included_count,p$parent_count,p$treated_count))})
  output$reference_readiness_ui_project<-renderUI({r<-reference_readiness(active(),repo_root);tags$dl(class="ytab-meta",tags$dt("Status"),tags$dd(r$label),tags$dt("FASTA"),tags$dd(r$fasta),tags$dt("Annotations"),tags$dd(paste(r$annotations,collapse="; ")),tags$dt("Bowtie2 prefix"),tags$dd(r$index_prefix),tags$dt("Index complete"),tags$dd(if(r$index_complete)"yes"else"no"))})
  output$project_files_ui<-renderUI({p<-active();tagList(tags$p("Project root: ",tags$code(p$project_root)),tags$p("Mapping: ",tags$code(file.path(p$project_root,"mapfastq"))),tags$p("Hit files: ",tags$code(file.path(p$project_root,"create_hit_file"))),tags$p("Raw summaries: ",tags$code(file.path(p$project_root,"summary"))),tags$p("Exports: ",tags$code(p$export_root)))})
  output$sample_sheet_table <- DT::renderDT({
    DT::datatable(compact_sample_table(active()$samples, sample_pipeline_status()),
                  rownames=FALSE,filter="top",selection="none",
                  options=list(scrollX=TRUE,pageLength=12))
  })
  output$terminal_command_ui <- renderUI(tags$pre(paste("python scripts/local/ytab_run_pipeline.py \\\n  --project-config",shQuote(active_project_path()),"\\\n  --profile core \\\n  --threads 2 \\\n  --keep-going")))
  output$help_ui <- renderUI({p<-active();tagList(tags$p("Core preprocessing ends when raw SummaryTable succeeds. MidLC is used only by the classifier branch; CPM is used inside treated-versus-parent analysis."),tags$p("Project outputs: ",tags$code(p$project_root)),tags$p("Exports: ",tags$code(p$export_root)),tags$p(class="ytab-warning","Smoke-project outputs demonstrate software behavior, not final biological conclusions."))})
  output$provenance_ui <- renderUI({p<-active();tagList(tags$p("Project status: ",tags$code(file.path(p$project_root,"manifests","project_status.json"))),tags$p("Logs: ",tags$code(file.path(p$project_root,"logs"))),tags$p("Repository commit: ",git_commit_or_unavailable(repo_root)))})
  output$last_job_summary <- renderUI({
    progress <- jobs$last_progress()
    if (is.null(progress)) return(tags$p(class="text-muted","No persisted job is available for this project."))
    tags$div(class="ytab-job-summary",
      tags$h4("Last job"),
      tags$p(tags$b(progress_status_label(progress$status %||% "unknown")),
        " · ", progress$stage %||% "pipeline",
        " · ", sprintf("%s/%s items", progress$processed_items %||% 0, progress$total_items %||% 0)),
      tags$p(class="text-muted", progress$message %||% ""),
      tags$p("Progress record: ",tags$code(jobs$current_job()$progress_file %||% "unavailable")))
  })
  output$project_created_card<-renderUI(if(project_created())tags$div(class="alert alert-success","Project is ready. Continue with the next available analysis step."))
  output$reference_readiness_ui<-renderUI({r<-reference_readiness(active(),repo_root);tags$dl(class="ytab-meta",tags$dt("Reference status"),tags$dd(r$label),tags$dt("FASTA"),tags$dd(if(r$fasta_found)paste("found —",r$fasta)else"missing"),tags$dt("Annotation"),tags$dd(if(r$annotation_found)"found"else"missing"),tags$dt("Bowtie2 index"),tags$dd(if(r$index_complete)paste("complete —",r$index_prefix)else"incomplete"),tags$dt("Mapping"),tags$dd(if(r$runnable)"ready"else"blocked"))})
  output$reference_actions_ui<-renderUI({r<-reference_readiness(active(),repo_root);tags$div(class="ytab-actions",actionButton("refresh_reference_status","Refresh reference status"),if(r$can_prepare&&!r$runnable)actionButton("prepare_reference","Prepare reference",class="btn-primary"),if(!r$fasta_found)tags$span(class="alert alert-warning","Reference files required"))})
  output$mapping_readiness <- renderUI(if(!reference_readiness(active(),repo_root)$runnable) blocked_card("Mapping blocked","Reference is not runnable.","Refresh or prepare the reference",active()$reference$reference_dir %||% "") else NULL)
  output$hits_readiness <- renderUI({s<-sample_pipeline_status();if(!any(s$included&s$hit_file%in%c("ready","success")))blocked_card("Hit-file creation blocked","Mapped BAM files are unavailable.","Run MapFastq",file.path(active()$project_root,"mapfastq"))})
  output$summary_readiness <- renderUI({s<-sample_pipeline_status();if(!any(s$included&s$summary%in%c("ready","success")))blocked_card("SummaryTable blocked","Hit files are unavailable.","Create hit files",file.path(active()$project_root,"create_hit_file"))})
  selected_names_ui<-function(){x<-preprocessing_selected_samples();tags$p(class="text-muted",sprintf("Selected samples: %d",length(x)))}
  for(prefix in c("mapfastq","create_hit_file","summary_table"))local({key<-prefix;output[[paste0(key,"_selected_names")]]<-renderUI(selected_names_ui())})
  for(prefix in c("mapfastq","create_hit_file","summary_table"))local({
    key<-prefix
    output[[paste0(key,"_run_scope_help")]]<-renderUI(tags$p(
      class="text-muted",
      "Run all always continues through every eligible sample. ",
      if(length(preprocessing_selected_samples())>1L)
        "Selected multi-sample runs continue by default after a recorded failure."
      else "A single-sample run stops when that sample finishes or fails."
    ))
  })
  counts_output<-function(prefix,column,stage){for(metric in c("selected","success","skipped","remaining","failed"))local({m<-metric;id<-paste0(prefix,"_",m);output[[id]]<-renderText({if(jobs$job_is_running())invalidateLater(1500,session);p<-jobs$current_progress();if(!is.null(p)&&identical(p$stage,stage)){values<-list(selected=p$total_items%||%0,success=p$successful_items%||%0,skipped=p$skipped_items%||%0,remaining=p$remaining_items%||%0,failed=p$failed_items%||%0);return(values[[m]])};rows<-sample_pipeline_status();values<-rows[[column]];if(m=="selected")length(preprocessing_selected_samples())else if(m=="success")sum(values%in%c("success","imported_success"))else if(m=="skipped")sum(values%in%c("skipped","imported_or_not_required"))else if(m=="failed")sum(values=="failed")else sum(!values%in%complete_stage_statuses)})})}
  counts_output("mapping","mapping","mapfastq");counts_output("hits","hit_file","create_hit_file");counts_output("summary","summary","summary")
  output$mapping_status_table<-DT::renderDT({if(jobs$job_is_running())invalidateLater(1500,session);s<-sample_pipeline_status();p<-jobs$current_progress();s$run_status<-s$mapping;s$current_phase<-"";s$current_elapsed<-"";s$last_update<-"";if(!is.null(p)&&identical(p$stage,"mapfastq")&&!is.null(p$current_item)){i<-match(p$current_item,s$sample);if(!is.na(i)){s$run_status[[i]]<-"running";s$current_phase[[i]]<-p$current_phase%||%"running";s$current_elapsed[[i]]<-format_duration(p$current_item_elapsed_seconds);s$last_update[[i]]<-p$updated_at%||%""}};d<-active()$samples;imported<-project_starts_from_hit_files(active());s$FASTQ<-if(imported)rep("not required",nrow(s))else ifelse(file.exists(as.character(d$fastq_1[match(s$sample,d$sample)])),"available","missing");s$BAM<-if(imported)rep("not required",nrow(s))else ifelse(file.exists(s$mapped_bam),"available","pending");s$Stats<-if(imported)rep("not required",nrow(s))else ifelse(file.exists(file.path(active()$project_root,"mapfastq",s$sample,paste0(s$sample,".mapping_stats.csv"))),"available","pending");DT::datatable(s[,c("sample","run_status","current_phase","current_elapsed","last_update","FASTQ","BAM","Stats","mapping_message")],rownames=FALSE,filter="top",selection="none",options=list(scrollX=TRUE,pageLength=12))})
  output$hits_status_table<-DT::renderDT({s<-sample_pipeline_status();imported<-project_starts_from_hit_files(active());visible<-data.frame(Sample=s$sample,BAM=if(imported)rep("not required",nrow(s))else ifelse(file.exists(s$mapped_bam),"available","missing"),`Hit file`=s$hit_file,Output=ifelse(file.exists(s$hits_path),"available","pending"),Elapsed=s$hit_elapsed,Message=s$hit_message,check.names=FALSE);DT::datatable(visible,rownames=FALSE,filter="top",selection="none",options=list(scrollX=TRUE,pageLength=12))})
  output$summary_status_table<-DT::renderDT({s<-sample_pipeline_status();visible<-data.frame(Sample=s$sample,`Hit file`=s$hit_file,`Summary table`=s$summary,`Feature table`=ifelse(file.exists(s$feature_table),"available","pending"),Elapsed=s$summary_elapsed,Message=s$summary_message,check.names=FALSE);DT::datatable(visible,rownames=FALSE,filter="top",selection="none",options=list(scrollX=TRUE,pageLength=12))})
  output$pipeline_matrix<-DT::renderDT({s<-sample_pipeline_status();DT::datatable(s[s$included,c("sample","fastq","mapping","bam_index","hit_file","summary")],rownames=FALSE,options=list(dom="t",paging=FALSE))})
  output$preprocessing_progress_ui<-renderUI({s<-sample_pipeline_status();make<-function(col,label){x<-stage_counts(s,col);pct<-if(x$total)round(100*x$complete/x$total)else 0;tags$div(tags$p(sprintf("%s: %d / %d",label,x$complete,x$total)),tags$div(class="progress",tags$div(class="progress-bar",style=sprintf("width:%d%%",pct))))};tagList(make("mapping","Mapping"),make("hit_file","Hit files"),make("summary","Summary tables"))})
  output$next_ready_action_ui<-renderUI({stage<-next_preprocessing_stage(active(),repo_root,sample_pipeline_status());tagList(tags$p(tags$b("Next ready action: "),stage),if(stage!="Complete")actionButton("run_next_ready_stage",paste("Open",stage),class="btn-primary")else actionButton("enter_analysis_browser","Enter analysis browser",class="btn-primary"),actionButton("refresh_preprocessing_status","Refresh status"))})
  output$analysis_ready_action<-renderUI({s<-sample_pipeline_status();x<-stage_counts(s,"summary");if(x$total>0&&x$complete==x$total)tags$div(class="alert alert-success","Preprocessing complete. Raw SummaryTable is available for all included samples.",actionButton("enter_analysis_browser_summary","Enter analysis browser",class="btn-primary"))})
  observeEvent(input$run_next_ready_stage,{destination<-resolve_continue_destination(active(),repo_root,sample_pipeline_status());go_to(destination$top,"preprocessing_tabs",destination$nested)},ignoreInit=TRUE)
  observeEvent(input$enter_analysis_browser,go_to("quality_control","qc_tabs","mapping_qc"),ignoreInit=TRUE);observeEvent(input$enter_analysis_browser_summary,go_to("quality_control","qc_tabs","mapping_qc"),ignoreInit=TRUE);observeEvent(input$overview_enter_analysis,go_to("quality_control","qc_tabs","mapping_qc"),ignoreInit=TRUE)
  observeEvent(input$overview_continue,{destination<-resolve_continue_destination(active(),repo_root,sample_pipeline_status());go_to(destination$top,"preprocessing_tabs",destination$nested)},ignoreInit=TRUE)
  output$qc_selected_samples<-renderUI(selected_names_ui())
  output$mapping_qc_readiness<-renderUI({if(!nrow(mapping_qc_data(active())))tags$p(class="text-muted","Mapping QC becomes available after MapFastq completes.")})
  output$summary_qc_readiness<-renderUI({if(!nrow(summary_qc_data(active())))tags$p(class="text-muted","Summary QC becomes available after raw SummaryTable completes.")})
  mapping_data<-reactive(mapping_qc_data(active()));summary_data<-reactive(summary_qc_data(active()))
  qc_download_slug<-function(...) {
    value<-paste(unlist(list(...)),collapse=".")
    value<-gsub("[^A-Za-z0-9_]+","_",value)
    value<-gsub("_+","_",value)
    value<-gsub("^_|_$","",value)
    if(nzchar(value)) value else "ytab"
  }
  qc_download_dimensions<-function(prefix,default_height="medium") {
    width_choice<-input[[paste0(prefix,"_plot_width")]]%||%"standard"
    height_choice<-input[[paste0(prefix,"_plot_height")]]%||%default_height
    css_height<-suppressWarnings(as.integer(sub("px$","",ytab_plot_height_px(height_choice,default=default_height))))
    if(is.na(css_height))css_height<-520L
    width<-switch(width_choice,wide=2000L,full=2400L,standard=1600L,1600L)
    list(width=width,height=max(700L,as.integer(css_height*2L)))
  }
  qc_download_plot_png<-function(file,prefix,expr,default_height="medium") {
    dim<-qc_download_dimensions(prefix,default_height)
    grDevices::png(file,width=dim$width,height=dim$height,res=150)
    on.exit(grDevices::dev.off(),add=TRUE)
    ytab_with_plot_display_options(input,prefix,expr)
  }
  qc_mapping_plot_data_for_download<-function(project) qc_mapping_stats_plot_data(project)
  qc_summary_plot_data_for_download<-function(project,choice,combined_group) {
    switch(choice,
           features=qc_summary_stats_plot_data(project),
           combined_features={
             data<-qc_summary_combined_features_hit_data(project)
             if(identical(combined_group,"both"))data[data$group%in%c("control","treated"),,drop=FALSE]
             else data[as.character(data$group)==combined_group,,drop=FALSE]
           },
           reads_per_hit=qc_summary_stats_plot_data(project),
           feature_intergenic=qc_summary_stats_plot_data(project),
           genome_bins=qc_summary_binned_plot_data(project),
           pairwise=qc_summary_pairwise_plot_data(project),
           qc_summary_stats_plot_data(project))
  }
  qc_render_current_summary_plot<-function(project,choice,combined_group) {
    switch(choice,
           features=plot_qc_summary_metric(project,"features"),
           combined_features=plot_qc_summary_combined_features_hit(project,combined_group),
           reads_per_hit=plot_qc_summary_metric(project,"reads_per_hit"),
           feature_intergenic=plot_qc_summary_metric(project,"feature_intergenic"),
           genome_bins=plot_qc_summary_genome_bins(project),
           pairwise=plot_qc_summary_pairwise_correlations(project),
           plot_qc_summary_metric(project,"complexity"))
  }
  qc_library_plot_data_for_download<-function(project,choice,group,panel) {
    data<-switch(choice,
                 jackpot=qc_library_summary_plot_data(project),
                 centromere=qc_library_centromere_plot_data(project),
                 metaplots=qc_library_metaplot_plot_data(project),
                 sequence_bias=qc_library_sequence_bias_plot_data(project),
                 qc_library_midlc_plot_data(project))
    if(!is.data.frame(data)||!nrow(data))return(data.frame())
    if(identical(choice,"metaplots")&&"feature"%in%names(data)&&!identical(panel,"all")) {
      data<-data[as.character(data$feature)==panel,,drop=FALSE]
    }
    qc_library_filter_data(project,data,group)
  }
  qc_render_current_library_plot<-function(project,choice,group,color_by,panel) {
    switch(choice,
           jackpot=plot_qc_library_jackpot_depth(project,group,color_by),
           centromere=plot_qc_library_centromere_bias(project,group,color_by),
           metaplots={
             panels<-if(identical(panel,"all"))intersect(c("tss","tts","trna"),unique(as.character(qc_library_metaplot_plot_data(project)$feature)))else panel
             panel_one<-panels[[1]]%||%"tss"
             plot_qc_library_metaplot(project,panel_one,group,color_by)
           },
           sequence_bias=plot_qc_library_sequence_bias(project,group),
           plot_qc_library_midlc(project,group,color_by))
  }
  output$mapping_qc_plot_controls<-renderUI({
    ytab_plot_customization_controls("mapping_qc", include_bars = TRUE,
                                     include_value_labels = TRUE,
                                     default_height = "medium",
                                     default_width = "wide",
                                     default_show_value_labels = FALSE)
  })
  output$mapping_qc_selected_plot<-renderUI(tagList(ytab_plot_frame(plotOutput("mapping_qc_stats_plot",width="100%",height=ytab_plot_height_px(input$mapping_qc_plot_height%||%"medium")),input$mapping_qc_plot_width%||%"standard","app-rendered"),uiOutput("mapping_qc_plot_key")))
  output$mapping_qc_plot_key<-renderUI({if(!identical(input$mapping_qc_plot_choice%||%"read_counts","read_counts"))return(NULL);tags$div(class="ytab-inline-legend",tags$span(class="ytab-inline-legend-item",tags$span(class="ytab-inline-swatch",style="background:#d8e4ef;border-color:#8ca8bf;"),"Total reads"),tags$span(class="ytab-inline-legend-item",tags$span(class="ytab-inline-swatch",style="background:#2f6f9f;border-color:#1f4f73;"),"Mapped reads"))})
  output$mapping_qc_stats_plot<-renderPlot(ytab_with_plot_display_options(input,"mapping_qc",plot_qc_mapping_stats(active(),input$mapping_qc_plot_choice%||%"read_counts")))
  output$mapping_qc_table<-DT::renderDT({data<-mapping_data();if(!nrow(data))return(NULL);compact_qc_table(data)})
  output$mapping_qc_details<-DT::renderDT({details<-attr(mapping_data(),"details");if(is.null(details)||!nrow(details))return(NULL);compact_qc_table(details)})
  output$download_mapping_qc_plot<-downloadHandler(filename=function()paste0(qc_download_slug(active()$project_id,"mapping_qc",input$mapping_qc_plot_choice%||%"read_counts"),".png"),content=function(file)qc_download_plot_png(file,"mapping_qc",plot_qc_mapping_stats(active(),input$mapping_qc_plot_choice%||%"read_counts"),"medium"))
  output$download_mapping_qc_plotted_data<-downloadHandler(filename=function()paste0(qc_download_slug(active()$project_id,"mapping_qc","mapping_summary"),".csv"),content=function(file)write.csv(qc_mapping_plot_data_for_download(active()),file,row.names=FALSE))
  output$download_mapping_qc_table<-downloadHandler(filename=function()"mapping_qc_summary.csv",content=function(file)write.csv(mapping_data(),file,row.names=FALSE))
  output$download_mapping_qc_details<-downloadHandler(filename=function()"mapping_qc_file_details.csv",content=function(file){details<-attr(mapping_data(),"details");if(is.null(details))details<-data.frame();write.csv(details,file,row.names=FALSE)})
  output$summary_qc_cards<-renderUI(summary_qc_cards(active()))
  output$summary_qc_combined_group_selector<-renderUI({
    if(!identical(input$summary_qc_plot_choice%||%"complexity","combined_features"))return(NULL)
    choices<-qc_summary_combined_feature_group_choices(active())
    if(!length(choices))return(tags$p(class="text-muted","No combined control/treated groups are available for this project."))
    selected<-input$summary_qc_combined_group%||%qc_summary_combined_feature_default_group(active())
    values<-unname(unlist(choices,use.names=FALSE))
    if(!selected%in%values)selected<-qc_summary_combined_feature_default_group(active())
    selectInput("summary_qc_combined_group","Combined group",choices=choices,selected=selected)
  })
  output$summary_qc_plot_controls<-renderUI({
    choice<-input$summary_qc_plot_choice%||%"complexity"
    bar_choices<-c("complexity","features","combined_features","reads_per_hit","feature_intergenic")
    ytab_plot_customization_controls("summary_qc",
                                     include_bars = choice %in% bar_choices,
                                     include_value_labels = choice %in% c(bar_choices,"pairwise"),
                                     default_height = if(choice%in%c("pairwise","genome_bins"))"large"else"medium",
                                     default_show_value_labels = FALSE)
  })
  output$summary_qc_selected_plot<-renderUI({
    choice<-input$summary_qc_plot_choice%||%"complexity"
    output_id<-switch(choice,features="summary_qc_features_plot",combined_features="summary_qc_combined_features_plot",reads_per_hit="summary_qc_reads_per_hit_plot",feature_intergenic="summary_qc_feature_intergenic_plot",genome_bins="summary_qc_genome_bins_plot",pairwise="summary_qc_pairwise_plot","summary_qc_complexity_plot")
    default_height<-if(choice%in%c("pairwise","genome_bins"))"large"else"medium"
    ytab_plot_frame(plotOutput(output_id,width="100%",height=ytab_plot_height_px(input$summary_qc_plot_height%||%default_height,default=default_height)),input$summary_qc_plot_width%||%"standard","app-rendered")
  })
  output$summary_qc_complexity_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_metric(active(),"complexity")))
  output$summary_qc_features_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_metric(active(),"features")))
  output$summary_qc_combined_features_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_combined_features_hit(active(),input$summary_qc_combined_group%||%qc_summary_combined_feature_default_group(active()))))
  output$summary_qc_reads_per_hit_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_metric(active(),"reads_per_hit")))
  output$summary_qc_feature_intergenic_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_metric(active(),"feature_intergenic")))
  output$summary_qc_genome_bins_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_genome_bins(active())))
  output$summary_qc_pairwise_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_pairwise_correlations(active())))
  output$summary_qc_table<-DT::renderDT({data<-summary_data();if(!nrow(data))return(NULL);compact_qc_table(data)})
  output$summary_qc_details<-DT::renderDT({details<-attr(summary_data(),"details");if(is.null(details)||!nrow(details))return(NULL);compact_qc_table(details)})
  output$download_summary_qc_plot<-downloadHandler(filename=function(){choice<-input$summary_qc_plot_choice%||%"complexity";group<-if(identical(choice,"combined_features"))input$summary_qc_combined_group%||%qc_summary_combined_feature_default_group(active())else"all";paste0(qc_download_slug(active()$project_id,"summary_qc",choice,group),".png")},content=function(file){choice<-input$summary_qc_plot_choice%||%"complexity";group<-input$summary_qc_combined_group%||%qc_summary_combined_feature_default_group(active());qc_download_plot_png(file,"summary_qc",qc_render_current_summary_plot(active(),choice,group),if(choice%in%c("pairwise","genome_bins"))"large"else"medium")})
  output$download_summary_qc_plotted_data<-downloadHandler(filename=function(){choice<-input$summary_qc_plot_choice%||%"complexity";group<-if(identical(choice,"combined_features"))input$summary_qc_combined_group%||%qc_summary_combined_feature_default_group(active())else"all";paste0(qc_download_slug(active()$project_id,"summary_qc",choice,group),".csv")},content=function(file){choice<-input$summary_qc_plot_choice%||%"complexity";group<-input$summary_qc_combined_group%||%qc_summary_combined_feature_default_group(active());write.csv(qc_summary_plot_data_for_download(active(),choice,group),file,row.names=FALSE)})
  output$download_summary_qc_table<-downloadHandler(filename=function()"summary_qc_table.csv",content=function(file)write.csv(summary_data(),file,row.names=FALSE))
  output$download_summary_qc_details<-downloadHandler(filename=function()"summary_qc_detailed_metrics.csv",content=function(file){details<-attr(summary_data(),"details");if(is.null(details))details<-data.frame();write.csv(details,file,row.names=FALSE)})
  output$library_diagnostics_group_selector<-renderUI({
    choices<-qc_library_group_choices(active())
    if(!length(choices))return(NULL)
    selected<-input$library_diagnostics_group%||%"all"
    values<-unname(unlist(choices,use.names=FALSE))
    if(!selected%in%values)selected<-"all"
    selectInput("library_diagnostics_group","Sample group",choices=choices,selected=selected)
  })
  output$library_diagnostics_metaplot_selector<-renderUI({
    if(!identical(input$library_diagnostics_plot_choice%||%"midlc","metaplots"))return(NULL)
    data<-qc_library_metaplot_plot_data(active())
    available<-intersect(c("tss","tts","trna"),unique(as.character(data$feature)))
    if(!length(available))return(NULL)
    labels<-c(tss="TSS",tts="TTS",trna="tRNA")[available]
    choices<-as.list(setNames(available,labels))
    if(length(available)>1L)choices<-c(list("All panels"="all"),choices)
    selectInput("library_diagnostics_metaplot_panel","Metaplot panel",choices=choices,selected=if(length(available)>1L)"all"else available[[1]])
  })
  output$library_diagnostics_plot_controls<-renderUI({
    choice<-input$library_diagnostics_plot_choice%||%"midlc"
    bar_choices<-c("jackpot","sequence_bias")
    ytab_plot_customization_controls("library_diagnostics",
                                     include_bars = choice %in% bar_choices,
                                     include_value_labels = choice %in% bar_choices,
                                     default_height = "medium",
                                     default_show_value_labels = FALSE)
  })
  output$library_diagnostics_selected_plot<-renderUI({
    choice<-input$library_diagnostics_plot_choice%||%"midlc"
    diagnostic_height<-ytab_plot_height_px(input$library_diagnostics_plot_height%||%"medium")
    group<-input$library_diagnostics_group%||%"all"
    if(identical(choice,"jackpot"))return(tagList(ytab_plot_frame(plotOutput("library_jackpot_depth_plot",width="100%",height=diagnostic_height),input$library_diagnostics_plot_width%||%"standard","app-rendered"),qc_library_sample_key_ui(active(),qc_library_summary_plot_data(active()),group),qc_library_plot_provenance_ui("Source: LibraryDiagnostics summary table.")))
    if(identical(choice,"centromere")){
      if(nrow(qc_library_centromere_plot_data(active())))return(tagList(ytab_plot_frame(plotOutput("library_centromere_bias_plot",width="100%",height=diagnostic_height),input$library_diagnostics_plot_width%||%"standard","app-rendered"),qc_library_sample_key_ui(active(),qc_library_centromere_plot_data(active()),group),qc_library_plot_provenance_ui("Source: existing centromere_bins.tsv tables. Archived static PNGs are not shown as current outputs.")))
      return(tags$p(class="text-muted","No live centromere-bias table is available for this project. Archived static diagnostic PNGs are not shown as current outputs."))
    }
    if(identical(choice,"metaplots")){
      if(nrow(qc_library_metaplot_plot_data(active()))){
        panel<-input$library_diagnostics_metaplot_panel%||%"all"
        panels<-if(identical(panel,"all"))intersect(c("tss","tts","trna"),unique(as.character(qc_library_metaplot_plot_data(active())$feature)))else panel
        plot_ids<-c(tss="library_metaplot_tss_plot",tts="library_metaplot_tts_plot",trna="library_metaplot_trna_plot")[panels]
        return(tagList(tags$div(class="ytab-stacked-plot-panels",lapply(plot_ids,function(id)ytab_plot_frame(plotOutput(id,width="100%",height=diagnostic_height),input$library_diagnostics_plot_width%||%"standard","app-rendered"))),qc_library_sample_key_ui(active(),qc_library_metaplot_plot_data(active()),group),qc_library_plot_provenance_ui("Source: existing TSS/TTS/tRNA metaplot TSV tables. Archived static PNGs are not shown as current outputs.")))
      }
      return(tags$p(class="text-muted","No live feature-metaplot table is available for this project. Archived static diagnostic PNGs are not shown as current outputs."))
    }
    if(identical(choice,"sequence_bias"))return(tagList(ytab_plot_frame(plotOutput("library_sequence_bias_plot",width="100%",height=diagnostic_height),input$library_diagnostics_plot_width%||%"standard","app-rendered"),qc_library_sample_key_ui(active(),qc_library_sequence_bias_plot_data(active()),group),qc_library_plot_provenance_ui("Source: existing sequence-bias TSV tables.")))
    tagList(ytab_plot_frame(plotOutput("library_midlc_plot",width="100%",height=diagnostic_height),input$library_diagnostics_plot_width%||%"standard","app-rendered"),qc_library_sample_key_ui(active(),qc_library_midlc_plot_data(active()),group),qc_library_plot_provenance_ui("Source: existing MidLC CSV tables."))
  })
  output$library_midlc_plot<-renderPlot(ytab_with_plot_display_options(input,"library_diagnostics",{invisible(qc_library_plot_cache_key(active(),"midlc",input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group",input$library_diagnostics_metaplot_panel%||%"",qc_selected_samples(),diagnostic_selected_run()));plot_qc_library_midlc(active(),input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group")}))
  output$library_jackpot_depth_plot<-renderPlot(ytab_with_plot_display_options(input,"library_diagnostics",{invisible(qc_library_plot_cache_key(active(),"jackpot",input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group",input$library_diagnostics_metaplot_panel%||%"",qc_selected_samples(),diagnostic_selected_run()));plot_qc_library_jackpot_depth(active(),input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group")}))
  output$library_sequence_bias_plot<-renderPlot(ytab_with_plot_display_options(input,"library_diagnostics",{invisible(qc_library_plot_cache_key(active(),"sequence_bias",input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group",input$library_diagnostics_metaplot_panel%||%"",qc_selected_samples(),diagnostic_selected_run()));plot_qc_library_sequence_bias(active(),input$library_diagnostics_group%||%"all")}))
  output$library_centromere_bias_plot<-renderPlot(ytab_with_plot_display_options(input,"library_diagnostics",{invisible(qc_library_plot_cache_key(active(),"centromere",input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group",input$library_diagnostics_metaplot_panel%||%"",qc_selected_samples(),diagnostic_selected_run()));plot_qc_library_centromere_bias(active(),input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group")}))
  output$library_metaplot_tss_plot<-renderPlot(ytab_with_plot_display_options(input,"library_diagnostics",{invisible(qc_library_plot_cache_key(active(),"metaplot",input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group","tss",qc_selected_samples(),diagnostic_selected_run()));plot_qc_library_metaplot(active(),"tss",input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group")}))
  output$library_metaplot_tts_plot<-renderPlot(ytab_with_plot_display_options(input,"library_diagnostics",{invisible(qc_library_plot_cache_key(active(),"metaplot",input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group","tts",qc_selected_samples(),diagnostic_selected_run()));plot_qc_library_metaplot(active(),"tts",input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group")}))
  output$library_metaplot_trna_plot<-renderPlot(ytab_with_plot_display_options(input,"library_diagnostics",{invisible(qc_library_plot_cache_key(active(),"metaplot",input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group","trna",qc_selected_samples(),diagnostic_selected_run()));plot_qc_library_metaplot(active(),"trna",input$library_diagnostics_group%||%"all",input$library_diagnostics_color_by%||%"group")}))
  output$download_library_diagnostics_plot<-downloadHandler(filename=function(){choice<-input$library_diagnostics_plot_choice%||%"midlc";group<-input$library_diagnostics_group%||%"all";panel<-if(identical(choice,"metaplots"))input$library_diagnostics_metaplot_panel%||%"all"else"";paste0(qc_download_slug(active()$project_id,"library_diagnostics",choice,group,panel),".png")},content=function(file){choice<-input$library_diagnostics_plot_choice%||%"midlc";group<-input$library_diagnostics_group%||%"all";panel<-if(identical(choice,"metaplots"))input$library_diagnostics_metaplot_panel%||%"all"else"";qc_download_plot_png(file,"library_diagnostics",qc_render_current_library_plot(active(),choice,group,input$library_diagnostics_color_by%||%"group",panel),"medium")})
  output$download_library_diagnostics_plotted_data<-downloadHandler(filename=function(){choice<-input$library_diagnostics_plot_choice%||%"midlc";group<-input$library_diagnostics_group%||%"all";panel<-if(identical(choice,"metaplots"))input$library_diagnostics_metaplot_panel%||%"all"else"";paste0(qc_download_slug(active()$project_id,"library_diagnostics",choice,group,panel),".csv")},content=function(file){choice<-input$library_diagnostics_plot_choice%||%"midlc";group<-input$library_diagnostics_group%||%"all";panel<-if(identical(choice,"metaplots"))input$library_diagnostics_metaplot_panel%||%"all"else"all";write.csv(qc_library_plot_data_for_download(active(),choice,group,panel),file,row.names=FALSE)})
  diagnostic_inventory<-reactive(diagnostic_inventory_state())
  observeEvent(input$refresh_diagnostic_files,{diagnostic_inventory_state(build_diagnostic_file_inventory(active()))},ignoreInit=TRUE)
  output$diagnostic_result_selector<-renderUI({x<-diagnostic_inventory();runs<-sort(unique(x$run_id));if(length(runs)>1){labels<-ifelse(runs=="current_legacy","Legacy diagnostics result",tools::toTitleCase(gsub("_"," ",runs)));selectInput("diagnostic_run_filter","Diagnostics result",choices=c(list("All results"="All"),as.list(setNames(unname(runs),unname(labels)))),selected=diagnostic_selected_run())}})
  observeEvent(input$diagnostic_run_filter,{diagnostic_selected_run(input$diagnostic_run_filter%||%"All")},ignoreInit=TRUE)
  output$diagnostic_file_filters<-renderUI({x<-diagnostic_inventory();if(!nrow(x))return(NULL);controls<-list();samples<-sort(unique(x$sample));types<-sort(unique(x$display_type[!(x$is_plot&x$extension%in%c("png","jpg","jpeg","svg"))]));sets<-sort(unique(x$sample_set));if(length(samples)>1)controls<-c(controls,list(selectInput("diagnostic_sample_filter","Sample",c("All",samples))));if(length(types)>1)controls<-c(controls,list(selectInput("diagnostic_type_filter","Type",c("All",types))));if(length(sets)>1)controls<-c(controls,list(selectInput("diagnostic_set_filter","Sample set",c("All",sets))));tags$div(class="ytab-filter-row",controls)})
  filtered_diagnostic_inventory<-reactive({x<-filter_diagnostic_inventory(diagnostic_inventory(),sample=input$diagnostic_sample_filter%||%"All",plot_type="All",sample_set=input$diagnostic_set_filter%||%"All",run_id=input$diagnostic_run_filter%||%diagnostic_selected_run(),type=input$diagnostic_type_filter%||%"All");if(!isTRUE(input$diagnostic_show_archived_static))x<-x[!(x$is_plot&x$extension%in%c("png","jpg","jpeg","svg")),,drop=FALSE];x})
  output$diagnostic_files_empty<-renderUI({x<-filtered_diagnostic_inventory();if(!nrow(diagnostic_inventory()))return(tags$p(class="text-muted","No diagnostic files are available for the selected run."));if(!nrow(x))return(tags$p(class="text-muted","No diagnostic files match the current filters. Archived/static generated images are hidden by default."))})
  output$diagnostic_files_table<-DT::renderDT({x<-filtered_diagnostic_inventory();if(!nrow(x))return(NULL);visible<-data.frame(File=x$filename,Type=x$display_type,`Sample set`=x$sample_set,`Sample or project`=x$sample,Size=x$size_display,Modified=x$modified,Action=ifelse(x$viewable,"View · Show path","Show path"),check.names=FALSE);compact_qc_table(visible)})

  qc_eligibility<-reactive(qc_sample_eligibility(active(),sample_pipeline_status()))
  observe({choices<-qc_eligibility()$Sample;selected<-intersect(qc_selected_samples(),choices);updateSelectizeInput(session,"qc_samples",choices=choices,selected=selected,server=TRUE)})
  observeEvent(input$qc_samples,qc_selected_samples(as.character(input$qc_samples%||%character())),ignoreInit=TRUE)
  output$library_diagnostics_readiness<-renderUI(if(!any(qc_eligibility()$Eligible=="Yes"))tags$p(class="ytab-warning","Library Diagnostics requires successful CreateHitFile outputs."))
  output$qc_sample_presets<-renderUI({x<-qc_eligibility();condition<-tolower(x$Condition);tags$div(class="ytab-actions",actionButton("qc_select_all","Select all eligible",class="btn-secondary"),if(any(condition=="parent"))actionButton("qc_select_parents","Select parents",class="btn-secondary"),if(any(condition=="treated"))actionButton("qc_select_treated","Select treated",class="btn-secondary"),actionButton("qc_clear_samples","Clear selection",class="btn-secondary"))})
  observeEvent(input$qc_select_all,qc_selected_samples(qc_eligibility()$Sample[qc_eligibility()$Eligible=="Yes"]),ignoreInit=TRUE)
  observeEvent(input$qc_select_parents,{x<-qc_eligibility();qc_selected_samples(x$Sample[x$Eligible=="Yes"&tolower(x$Condition)=="parent"])},ignoreInit=TRUE)
  observeEvent(input$qc_select_treated,{x<-qc_eligibility();qc_selected_samples(x$Sample[x$Eligible=="Yes"&tolower(x$Condition)=="treated"])},ignoreInit=TRUE)
  observeEvent(input$qc_clear_samples,qc_selected_samples(character()),ignoreInit=TRUE)
  output$qc_sample_preview<-DT::renderDT({x<-qc_eligibility();selected<-qc_selected_samples();x$Selected<-ifelse(x$Sample%in%selected,"Yes","No");compact_qc_table(x[c("Sample","Condition","Background","Pool","Hit-file status","Eligible","Reason","Selected")])})
  qc_resolution<-reactive({x<-qc_eligibility();selected<-qc_selected_samples();eligible<-intersect(selected,x$Sample[x$Eligible=="Yes"]);list(selected=selected,eligible=eligible,excluded=setdiff(selected,eligible),label=if(length(eligible)&&setequal(eligible,x$Sample[x$Eligible=="Yes"]))"All eligible samples"else if(length(eligible)&&all(tolower(x$Condition[match(eligible,x$Sample)])=="parent"))"Parents"else if(length(eligible)&&all(tolower(x$Condition[match(eligible,x$Sample)])=="treated"))"Treated"else"Custom selection")})
  available_diagnostic_runs<-reactive({status_tick();discover_diagnostics_runs(active())})
  qc_result_state<-reactive({r<-qc_resolution();runs<-available_diagnostic_runs();valid<-Filter(function(x)validate_diagnostics_run_inputs(x,r$eligible,active()),runs);signature<-if(length(valid))valid[[1]]$cache_signature else"";resolve_qc_result_state(r$eligible,signature,runs,sum(qc_eligibility()$Eligible=="Yes"),if(length(r$eligible))r$label else"No selection")})
  output$qc_current_selection<-renderUI({s<-qc_result_state();tags$div(class="ytab-selection-summary",tags$b("Current selection"),tags$p(sprintf("%d of %d eligible samples selected",s$current_selection_count,s$eligible_count)),tags$div(class="ytab-project-preview",tags$span(class="ytab-preview-badge",sprintf("%d selected",s$current_selection_count)),if(s$current_selection_count)tags$span(class="ytab-preview-badge",s$current_selection_label),tags$span(class="ytab-preview-badge",sprintf("%d eligible",s$eligible_count))),if(!s$current_selection_count)tags$p(class="text-muted","No samples are currently selected."))})
  output$qc_selection_message<-renderUI({r<-qc_resolution();if(!length(r$selected))return(tags$p(class="ytab-warning","Select at least one eligible sample."));if(!length(r$eligible))return(tags$p(class="ytab-warning","None of the selected samples has a valid CreateHitFile output."));if(length(r$excluded))tags$p(class="ytab-warning",sprintf("Running diagnostics on %d eligible samples. %d selected samples were excluded because hit files were unavailable.",length(r$eligible),length(r$excluded)))else tags$p(sprintf("%s: %d eligible samples selected.",r$label,length(r$eligible)))})
  matching_cache<-reactive(qc_result_state()$matching_result)
  output$library_diagnostics_cache_ui<-renderUI({if(!identical(input$library_diagnostics_mode,"run"))return(NULL);if(!is.null(matching_cache()))tags$div(class="alert alert-info",tags$b("Cached diagnostics available"),tags$p("A matching successful result already exists for this sample set."))else{base<-file.path(active()$project_root,"manifests","library_diagnostics","runs");if(dir.exists(base)&&length(list.files(base,pattern="manifest\\.json$",recursive=TRUE)))tags$p(class="text-muted","Existing diagnostics were found, but they were generated from a different sample selection.")}})
  output$qc_execution_summary<-renderUI({r<-qc_resolution();tags$dl(class="ytab-execution-summary",tags$dt("Samples"),tags$dd(length(r$eligible)),tags$dt("Sample set"),tags$dd(if(length(r$eligible))r$label else"None"),tags$dt("Matching cache"),tags$dd(if(is.null(matching_cache()))"No"else"Yes"),tags$dt("Execution"),tags$dd(if(identical(input$library_diagnostics_mode,"run"))"Run diagnostics"else"Preview command"))})
  output$library_diagnostics_action<-renderUI({r<-qc_resolution();mode<-input$library_diagnostics_mode%||%"preview";force<-FALSE;actionButton("run_library_diagnostics",diagnostic_action_label(mode,!is.null(matching_cache()),force),class="btn-primary",disabled=if(!length(r$eligible)||jobs$job_is_running())"disabled"else NULL)})

  fitness_design<-reactive({fitness_design_tick();read_fitness_design(active())});fitness_issues<-reactive({fitness_design_tick();read_fitness_design_issues(active())})
  observeEvent(active_project_path(),{path<-active_project_path();if(!identical(path,fitness_selection_project())){design<-fitness_design();fitness_selected_comparisons(fitness_default_selection(design));fitness_selection_project(path);fitness_selected_result("");updateTextInput(session,"fitness_analysis_id",value=fitness_default_analysis_id(design))}},ignoreInit=FALSE)
  observe({design<-fitness_design();contract<-normalize_fitness_design_columns(design);if(!contract$valid)return();valid<-unname(as.character(contract$data$comparison_id[contract$data$valid%||%FALSE]));updateSelectizeInput(session,"fitness_comparisons",choices=as.list(setNames(valid,valid)),selected=unname(intersect(fitness_selected_comparisons(),valid)),server=TRUE)})
  observeEvent(input$fitness_comparisons,fitness_selected_comparisons(unname(as.character(input$fitness_comparisons%||%character()))),ignoreInit=TRUE)
  output$fitness_design_state<-renderUI({design<-fitness_design();error<-attr(design,"design_error")%||%"";if(nzchar(error))return(tags$div(class="alert alert-danger",error));if(!nrow(design))tags$p(class="text-muted","No comparison design has been generated for this project.")else tags$p(sprintf("%d valid parent-treated comparisons were detected.",sum(design$valid)))})
  output$fitness_design_cards<-renderUI({x<-fitness_design_summary(fitness_design(),fitness_issues());tags$div(class="ytab-stat-grid",tags$div(tags$b(x$valid),"Valid comparisons"),tags$div(tags$b(x$included),"Included comparisons"),tags$div(tags$b(x$incomplete),"Incomplete comparisons"),tags$div(tags$b(x$issues),"Design issues"))})
  output$fitness_design_generate_action<-renderUI(if(!file.exists(fitness_design_path(active())))actionButton("fitness_generate_design","Generate design",class="btn-primary"))
  output$fitness_design_presets<-renderUI({design<-fitness_design();if(!nrow(design))return(NULL);backgrounds<-sort(unique(design$background[design$valid]));tags$div(class="ytab-actions",actionButton("fitness_design_all","Select all valid",class="btn-secondary"),actionButton("fitness_design_none","Select none",class="btn-secondary"),lapply(backgrounds,function(bg)actionButton(paste0("fitness_background_",make.names(bg)),paste("Select",bg,"background"),class="btn-secondary")))})
  observeEvent(input$fitness_design_all,fitness_selected_comparisons(fitness_default_selection(fitness_design())),ignoreInit=TRUE);observeEvent(input$fitness_design_none,fitness_selected_comparisons(character()),ignoreInit=TRUE)
  observe({backgrounds<-sort(unique(fitness_design()$background[fitness_design()$valid]));lapply(backgrounds,function(bg)local({value<-bg;observeEvent(input[[paste0("fitness_background_",make.names(value))]],fitness_selected_comparisons(as.character(fitness_design()$comparison_id[fitness_design()$valid&fitness_design()$background==value])),ignoreInit=TRUE)}))})
  output$fitness_design_table<-DT::renderDT({x<-fitness_design_table(fitness_design(),active()$project_root);if(!nrow(x))return(NULL);compact_qc_table(x)})
  output$fitness_design_issues_table<-DT::renderDT({x<-fitness_issues();if(!nrow(x))return(DT::datatable(data.frame(Message="No design issues were detected."),rownames=FALSE,options=list(dom="t")));compact_qc_table(x)})
  output$fitness_design_input_details<-DT::renderDT({x<-fitness_design();if(!nrow(x))return(NULL);details<-data.frame(Comparison=x$comparison_id,`Parent summary`=vapply(x$parent_summary_path,relative_project_path,"",project_root=active()$project_root),`Treated summary`=vapply(x$treated_summary_path,relative_project_path,"",project_root=active()$project_root),check.names=FALSE);compact_qc_table(details)})
  output$fitness_selection_summary<-renderUI({design<-fitness_design();selected<-fitness_selected_comparisons();tags$p(sprintf("%d of %d valid comparisons selected",length(selected),sum(design$valid)))})
  observeEvent(input$fitness_select_all,fitness_selected_comparisons(fitness_default_selection(fitness_design())),ignoreInit=TRUE);observeEvent(input$fitness_clear_selection,fitness_selected_comparisons(character()),ignoreInit=TRUE)
  output$fitness_selected_preview<-DT::renderDT({resolved<-resolve_selected_fitness_comparisons(fitness_design(),fitness_selected_comparisons());if(!resolved$valid||!nrow(resolved$rows))return(NULL);x<-resolved$rows;optional<-function(name)if(name%in%names(x))x[[name]]else rep("",nrow(x));compact_qc_table(data.frame(Comparison=x$comparison_id,Parent=x$parent_sample,Treated=x$treated_sample,Background=optional("background"),Pool=optional("pool"),`Input status`=optional("status"),check.names=FALSE))})
  output$fitness_classifier_controls<-renderUI({if(!isTRUE(input$fitness_annotate_classifier))return(NULL);root<-file.path(active()$project_root,"classifier");targets<-if(dir.exists(root))basename(list.dirs(root,recursive=FALSE))else character();tags<-targets[vapply(targets,essentiality_valid_tag,logical(1))];if(!length(tags))return(tags$div(class="alert alert-warning","No classifier result is available. Disable annotation to run fitness analysis without it."));choices<-c(list("Recommended target"="recommended"),as.list(setNames(tags,tags)));tagList(selectInput("fitness_classifier_target","Classifier target",choices=choices),helpText("Classifier annotation adds parent essentiality labels to the fitness result. It does not change CPM normalization, effect sizes, z-scores, or fitness calls."))})
  fitness_results<-reactive({status_tick();discover_fitness_results(active())});fitness_cache_signature<-reactive(fitness_current_cache_signature(fitness_results(),fitness_selected_comparisons(),active(),input$fitness_analysis_id%||%fitness_default_analysis_id(fitness_design()),isTRUE(input$fitness_annotate_classifier),input$fitness_classifier_target%||%""));fitness_result_state<-reactive(resolve_fitness_result_state(fitness_selected_comparisons(),fitness_cache_signature(),fitness_results(),nrow(fitness_design())>0L));fitness_matching<-reactive(fitness_result_state()$matching_result)
  output$fitness_run_summary<-renderUI({selected<-fitness_selected_comparisons();aid<-input$fitness_analysis_id%||%fitness_default_analysis_id(fitness_design());tags$dl(class="ytab-execution-summary",tags$dt("Analysis"),tags$dd(gsub("_"," ",aid)),tags$dt("Comparison count"),tags$dd(length(selected)),tags$dt("Comparison scope"),tags$dd(if(length(selected)==sum(fitness_design()$valid))"All valid comparisons"else"Selected subset"),tags$dt("Input source"),tags$dd("Raw per-sample SummaryTable"),tags$dt("Normalization"),tags$dd("CPM inside R"),tags$dt("Classifier annotation"),tags$dd(if(isTRUE(input$fitness_annotate_classifier))"Yes"else"No"),tags$dt("Matching result"),tags$dd(if(is.null(fitness_matching()))"No"else"Yes"),tags$dt("Expected output"),tags$dd(tags$code(file.path("treated_vs_parent",aid))))})
  output$fitness_run_action<-renderUI({mode<-input$fitness_execution_mode%||%"preview";force<-FALSE;error<-validate_fitness_analysis_id(input$fitness_analysis_id%||%"");actionButton("fitness_run",fitness_action_label(mode,!is.null(fitness_matching()),force),class="btn-primary",disabled=if(!length(fitness_selected_comparisons())||nzchar(error)||jobs$job_is_running())"disabled"else NULL)})
  observeEvent(input$fitness_continue_analysis,go_to("fitness","fitness_tabs","fitness_run"),ignoreInit=TRUE);observeEvent(input$fitness_change_design,go_to("fitness","fitness_tabs","fitness_design"),ignoreInit=TRUE)
  output$essentiality_readiness<-renderUI({s<-sample_pipeline_status();if(!any(s$hit_file%in%c("success","skipped")))blocked_card("Essentiality preprocessing blocked","Raw hit files are unavailable.","Create hit files")})
  output$fitness_readiness<-renderUI({s<-sample_pipeline_status();if(!all(s$summary[s$included]%in%c("success","skipped")))blocked_card("Fitness analysis blocked","Fitness analysis requires raw SummaryTable outputs.","Complete raw SummaryTable")else if(!file.exists(file.path(active()$project_root,"config","comparison_design.csv")))blocked_card("Comparison design required","Raw summaries are ready.","Generate comparison design")})
  output$results_readiness<-renderUI({p<-active();if(!dir.exists(p$export_root)||!length(list.files(p$export_root)))tags$p(class="text-muted","Reports and stable result exports become available after preprocessing and downstream analyses.")})
  output$results_exports_roots<-renderUI({p<-active();tags$dl(class="ytab-meta",tags$dt("Project root"),tags$dd(tags$code(p$project_root)),tags$dt("Export root"),tags$dd(tags$code(p$export_root)))})
  output$results_exports_table<-DT::renderDT({p<-active();roots<-c("Generated pipeline plots"=file.path(p$project_root,"treated_vs_parent"),"Result tables"=p$project_root,"Manifests"=file.path(p$project_root,"manifests"),"Reports"=file.path(p$export_root,"report"),"Export bundles"=file.path(p$export_root,"bundles"));rows<-list();for(group in names(roots)){root<-roots[[group]];if(!dir.exists(root))next;files<-list.files(root,recursive=TRUE,full.names=TRUE);files<-files[file.exists(files)&!dir.exists(files)];if(identical(group,"Result tables"))files<-files[grepl("\\.(csv|tsv)$",files,ignore.case=TRUE)];if(identical(group,"Generated pipeline plots"))files<-files[grepl("\\.(png|pdf|svg)$",files,ignore.case=TRUE)];if(!length(files))next;info<-file.info(files);rows[[length(rows)+1L]]<-data.frame(Group=group,File=basename(files),Type=toupper(tools::file_ext(files)),Size=human_file_size(info$size),Modified=format(info$mtime,"%Y-%m-%d %H:%M"),Path=normalizePath(files,winslash="/",mustWork=FALSE),check.names=FALSE,stringsAsFactors=FALSE)};data<-if(length(rows))do.call(rbind,rows)else data.frame(Group="No files available",File="",Type="",Size="",Modified="",Path="",check.names=FALSE);data<-data[seq_len(min(nrow(data),300L)),,drop=FALSE];DT::datatable(data,rownames=FALSE,filter="top",selection="none",options=list(pageLength=12,scrollX=TRUE,columnDefs=list(list(targets=5,visible=FALSE))))})

  python_bin <- locate_python_executable
  run_cli <- function(script,args=character()) {
    checked<-validate_project_config_path(active_project_path(),repo_root);if(!checked$valid){showNotification(checked$message,type="error");return(invisible(NULL))}
    full<-c(file.path(repo_root,"scripts","local",script),"--project-config",checked$path,args)
    result<-tryCatch(run_process_sync(python_bin(),full,wd=repo_root),error=function(e)list(status=1L,stdout="",stderr=conditionMessage(e)));log_text(paste(c(format_command_for_display(python_bin(),full),result$stdout,result$stderr,paste("Exit status:",result$status)),collapse="\n"));invisible(result)
  }
  samples_arg <- function(value) if(length(value)) c("--samples",paste(value,collapse=",")) else character()
  dry_arg <- function(value) if(isTRUE(value)) "--dry-run" else character()
  start_stage_job<-function(stage,script,selected,dry,force=FALSE,keep=FALSE){if(jobs$job_is_running()){showNotification("A pipeline job is already running.",type="warning");return(FALSE)};if(!length(selected)){showNotification("No eligible samples are available for this stage.",type="warning");return(FALSE)};threads<-as.integer(input$preprocessing_threads%||%2L);if(threads<2L)threads<-2L;args<-c(file.path(repo_root,"scripts/local",script),"--project-config",active_project_path(),"--threads",as.character(threads),samples_arg(selected),dry_arg(dry),if(force)"--force",if(keep)"--keep-going");jobs$start_job(python_bin(),args,active_project_path(),stage,wd=repo_root,selected_items=selected);TRUE}
  resolve_stage_request<-function(stage,all=FALSE,force=FALSE){requested<-if(all)included_samples()else isolate(preprocessing_selected_samples());if(!all&&!length(requested)){showNotification("No samples are selected. Open the Samples tab to choose samples.",type="warning");return(NULL)};resolution<-resolve_selected_preprocessing_samples(requested,active()$samples,sample_pipeline_status(),stage,force);if(!length(resolution$eligible)){message<-if(length(resolution$already_complete)==length(resolution$included)&&length(resolution$included))switch(stage,create_hit_file="All selected samples already have successful hit files.",summary="All selected samples already have successful summary tables.",mapping="All selected samples already have successful mapping outputs.")else if(stage=="create_hit_file"&&length(resolution$blocked))"The selected samples do not yet have mapped BAM files."else if(stage=="summary"&&length(resolution$blocked))"The selected samples do not yet have hit files."else"No eligible samples are available for this stage.";showNotification(message,type="warning");return(NULL)};excluded<-length(resolution$blocked)+length(resolution$missing);if(excluded)showNotification(sprintf("Running %d eligible samples. %d requested samples were excluded because prerequisites were unavailable.",length(resolution$eligible),excluded),type="warning");resolution}
  launch_stage<-function(stage,script,all,dry,force,keep){resolution<-resolve_stage_request(stage,all,force);if(!is.null(resolution)){selected<-if(all)unique(c(resolution$eligible,resolution$already_complete))else resolution$eligible;manager_stage<-if(stage=="mapping")"mapfastq"else stage;effective_keep<-resolve_keep_going(manager_stage,length(selected),all,keep);start_stage_job(manager_stage,script,selected,dry,force,effective_keep)}}
  observeEvent(input$preview_mapfastq,launch_stage("mapping","ytab_run_mapfastq.py",FALSE,TRUE,FALSE,TRUE),ignoreInit=TRUE)
  observeEvent(input$run_mapfastq,launch_stage("mapping","ytab_run_mapfastq.py",FALSE,FALSE,FALSE,TRUE),ignoreInit=TRUE)
  observeEvent(input$run_all_mapfastq,launch_stage("mapping","ytab_run_mapfastq.py",TRUE,FALSE,FALSE,TRUE),ignoreInit=TRUE)
  observeEvent(input$preview_create_hit_file,launch_stage("create_hit_file","ytab_run_create_hit_file.py",FALSE,TRUE,FALSE,TRUE),ignoreInit=TRUE)
  observeEvent(input$run_create_hit_file,launch_stage("create_hit_file","ytab_run_create_hit_file.py",FALSE,FALSE,FALSE,TRUE),ignoreInit=TRUE)
  observeEvent(input$run_all_create_hit_file,launch_stage("create_hit_file","ytab_run_create_hit_file.py",TRUE,FALSE,FALSE,TRUE),ignoreInit=TRUE)
  observeEvent(input$preview_summary_table,launch_stage("summary","ytab_run_summary_table.py",FALSE,TRUE,FALSE,TRUE),ignoreInit=TRUE)
  observeEvent(input$run_summary_table,launch_stage("summary","ytab_run_summary_table.py",FALSE,FALSE,FALSE,TRUE),ignoreInit=TRUE)
  observeEvent(input$run_all_summary_table,launch_stage("summary","ytab_run_summary_table.py",TRUE,FALSE,FALSE,TRUE),ignoreInit=TRUE)
  launch_library_diagnostics<-function(){r<-isolate(qc_resolution());mode<-isolate(input$library_diagnostics_mode%||%"preview");force<-FALSE;start_stage_job("library_diagnostics","ytab_run_library_diagnostics.py",r$eligible,mode=="preview",force,FALSE)}
  observeEvent(input$run_library_diagnostics,{
    r<-qc_resolution();req(length(r$eligible));mode<-input$library_diagnostics_mode%||%"preview";force<-FALSE
    if(mode=="preview"||(!is.null(matching_cache())&&!force)){launch_library_diagnostics();return()}
    output_path<-file.path("library_diagnostics","runs",tolower(gsub("[^A-Za-z0-9]+","_",r$label)))
    showModal(modalDialog(title="Run Library Diagnostics?",tags$dl(class="ytab-meta",tags$dt("Sample set"),tags$dd(r$label),tags$dt("Eligible samples"),tags$dd(length(r$eligible)),tags$dt("Execution"),tags$dd("Real analysis"),tags$dt("Matching cache"),tags$dd(if(is.null(matching_cache()))"No"else"Yes"),tags$dt("Output"),tags$dd(tags$code(output_path))),tags$details(tags$summary("Selected sample details"),tags$p(paste(r$eligible,collapse=", "))),footer=tagList(modalButton("Cancel"),actionButton("confirm_library_diagnostics","Run diagnostics",class="btn-primary"))))
  },ignoreInit=TRUE)
  observeEvent(input$confirm_library_diagnostics,{removeModal();launch_library_diagnostics()},ignoreInit=TRUE)
  observeEvent(input$prepare_reference,{p<-active();if(reference_readiness(p,repo_root)$runnable){showNotification("Reference is already complete; rebuilding was skipped.");return()};args<-c(file.path(repo_root,"scripts/local/ytab_prepare_reference.py"),"--species",p$species,"--threads",as.character(p$threads),"--repo-root",repo_root);result<-run_process_sync(python_bin(),args,wd=repo_root);log_text(paste(c(format_command_for_display(python_bin(),args),result$stdout,result$stderr),collapse="\n"));state$load(active_project_path(),repo_root);status_tick(status_tick()+1L)},ignoreInit=TRUE)
  essentiality_server(input,output,session,active,active_project_path,repo_root,
                      jobs,status_tick,go_to,python_bin,log_text)
  gene_domain_explorer_server(input,output,session,active,active_project_path,
                              repo_root,python_bin,log_text)
  genome_browser_tracks <- reactive({
    status_tick()
    genome_browser_track_inventory(active())
  })
  genome_browser_selected_rows <- reactive({
    tracks <- genome_browser_tracks()
    if (!nrow(tracks)) return(data.frame())
    preset <- input$genome_browser_track_preset %||% "all"
    if (identical(preset, "custom")) {
      selected <- input$genome_browser_tracks %||% genome_browser_custom_tracks()
      return(tracks[tracks$sample %in% selected, , drop = FALSE])
    }
    genome_browser_preset_rows(preset, tracks)
  })
  observeEvent(genome_browser_tracks(), {
    tracks <- genome_browser_tracks()
    choices <- if (nrow(tracks)) as.list(stats::setNames(tracks$sample, tracks$display_label)) else list()
    presets <- genome_browser_preset_choices(tracks)
    selected_preset <- input$genome_browser_track_preset %||% "all"
    if (!selected_preset %in% unname(presets)) selected_preset <- "all"
    updateSelectInput(session, "genome_browser_track_preset", choices = presets, selected = selected_preset)
    selected <- genome_browser_preset_rows(selected_preset, tracks)$sample
    updateSelectizeInput(session, "genome_browser_tracks", choices = choices, selected = selected, server = TRUE)
    updateTextInput(session, "genome_browser_locus", value = genome_browser_default_locus(active()))
  }, ignoreInit = FALSE)
  observeEvent(input$genome_browser_track_preset, {
    tracks <- genome_browser_tracks()
    if (!nrow(tracks)) return()
    preset <- input$genome_browser_track_preset %||% "all"
    selected <- if (identical(preset, "custom")) {
      current <- genome_browser_custom_tracks()
      if (length(current)) current else tracks$sample
    } else genome_browser_preset_rows(preset, tracks)$sample
    updateSelectizeInput(session, "genome_browser_tracks", selected = selected)
  }, ignoreInit = TRUE)
  observeEvent(input$genome_browser_tracks, {
    selected <- input$genome_browser_tracks %||% character()
    preset <- input$genome_browser_track_preset %||% ""
    if (!identical(preset, "custom")) {
      expected <- genome_browser_preset_rows(preset, genome_browser_tracks())$sample
      if (!setequal(selected, expected)) updateSelectInput(session, "genome_browser_track_preset", selected = "custom")
    }
    genome_browser_custom_tracks(selected)
  }, ignoreInit = TRUE)
  output$genome_browser_status <- renderUI({
    tracks <- genome_browser_tracks()
    reference <- genome_browser_reference_config(active())
    if (is.null(reference)) return(tags$p(class = "ytab-warning", "Indexed reference FASTA is not available for IGV."))
    if (!nrow(tracks)) return(tags$p(class = "ytab-warning", "No CreateHitFile browser tracks were found for this project."))
    tags$p(class = "text-muted", sprintf("%d insertion tracks discovered.", nrow(tracks)))
  })
  output$genome_browser_preset_message <- renderUI({
    rows <- genome_browser_selected_rows()
    if (!nrow(rows)) tags$p(class = "ytab-warning", "No tracks match this preset for the current project.") else NULL
  })
  output$genome_browser_lane_key <- renderUI({
    rows <- genome_browser_selected_rows()
    if (!nrow(rows)) return(NULL)
    lane_badges <- lapply(seq_len(nrow(rows)), function(i) {
      role <- rows$role[[i]] %||% ""
      css <- if (identical(role, "treated")) {
        "ytab-igv-lane-badge ytab-igv-lane-treated"
      } else if (identical(role, "parent")) {
        "ytab-igv-lane-badge ytab-igv-lane-parent"
      } else {
        "ytab-igv-lane-badge"
      }
      tags$span(
        class = css,
        tags$span(class = "ytab-igv-lane-number", paste0(i, ".")),
        rows$display_label[[i]] %||% rows$sample[[i]]
      )
    })
    tags$div(
      class = "ytab-igv-lane-key",
      tags$span(class = "ytab-igv-lane-key-title", "Displayed lanes"),
      lane_badges,
      tags$span(class = "ytab-igv-lane-badge ytab-igv-lane-gene", "Gene annotations")
    )
  })
  output$genome_browser_track_table <- DT::renderDT({
    rows <- genome_browser_selected_rows()
    if (!nrow(rows)) return(DT::datatable(data.frame(Message = "No tracks selected."), rownames = FALSE))
    visible <- rows[, intersect(c("display_label", "role", "pool", "format", "source_file"), names(rows)), drop = FALSE]
    names(visible) <- c("Track", "Role", "Pool", "Format", "Source file")[seq_along(visible)]
    DT::datatable(visible, rownames = FALSE, selection = "none", options = list(scrollX = TRUE, pageLength = 10))
  })
  output$genome_browser_technical <- renderText({
    rows <- genome_browser_selected_rows()
    reference <- genome_browser_reference_config(active())
    paste(c(
      paste("Project:", active()$project_id),
      paste("Reference available:", if (is.null(reference)) "no" else "yes"),
      paste("Selected tracks:", nrow(rows)),
      paste("Locus:", input$genome_browser_locus %||% "")
    ), collapse = "\n")
  })
  genome_browser_load <- function() {
    if (!identical(isolate(input$workspace_tabs %||% ""), "genome_browser"))
      return(invisible(FALSE))
    project <- isolate(active())
    rows <- isolate(genome_browser_selected_rows())
    locus <- isolate(input$genome_browser_locus %||% "")
    config <- genome_browser_session_config(project, rows, locus)
    if (is.null(config)) {
      session$sendCustomMessage("ytab_igv_warning", "Genome browser cannot load because the indexed reference FASTA is unavailable.")
      return(invisible(FALSE))
    }
    session$sendCustomMessage("ytab_igv_init", config)
    invisible(TRUE)
  }
  observeEvent(active_project_path(), {
    session$onFlushed(function() genome_browser_load(), once = TRUE)
  }, ignoreInit = FALSE)
  observeEvent(input$workspace_tabs, {
    if (identical(input$workspace_tabs %||% "", "genome_browser"))
      session$onFlushed(function() genome_browser_load(), once = TRUE)
  }, ignoreInit = TRUE)
  observeEvent(input$genome_browser_reload, genome_browser_load(), ignoreInit = TRUE)
  observeEvent(input$genome_browser_go, {
    locus <- trimws(input$genome_browser_locus %||% "")
    if (nzchar(locus)) session$sendCustomMessage("ytab_igv_goto", locus)
  }, ignoreInit = TRUE)
  observeEvent(input$genome_browser_pan_left, {
    session$sendCustomMessage("ytab_igv_pan", -1L)
  }, ignoreInit = TRUE)
  observeEvent(input$genome_browser_pan_right, {
    session$sendCustomMessage("ytab_igv_pan", 1L)
  }, ignoreInit = TRUE)
  observeEvent(input$genome_browser_clear_tracks, {
    updateSelectInput(session, "genome_browser_track_preset", selected = "custom")
    genome_browser_custom_tracks(character())
    updateSelectizeInput(session, "genome_browser_tracks", selected = character())
  }, ignoreInit = TRUE)
  observeEvent(genome_browser_selected_rows(), {
    genome_browser_load()
  }, ignoreInit = TRUE)
  comparative_server(input,output,session,repo_root)
  refresh_fitness_design<-function(overwrite=TRUE){result<-run_cli("ytab_init_comparison_design.py",c(if(overwrite)"--overwrite","--print-design"));fitness_design_tick(fitness_design_tick()+1L);design<-fitness_design();fitness_selected_comparisons(fitness_default_selection(design));updateTextInput(session,"fitness_analysis_id",value=fitness_default_analysis_id(design));result}
  observeEvent(input$fitness_generate_design,refresh_fitness_design(FALSE),ignoreInit=TRUE);observeEvent(input$fitness_refresh_design,refresh_fitness_design(TRUE),ignoreInit=TRUE)
  observeEvent(input$fitness_review_issues,{x<-fitness_issues();showModal(modalDialog(title="Comparison design issues",if(nrow(x))DT::dataTableOutput("fitness_issues_modal")else"No design issues were detected.",easyClose=TRUE));if(nrow(x))output$fitness_issues_modal<-DT::renderDataTable(compact_qc_table(x))},ignoreInit=TRUE)
  start_fitness_job<-function(){resolved<-resolve_selected_fitness_comparisons(isolate(fitness_design()),isolate(fitness_selected_comparisons()));aid<-isolate(input$fitness_analysis_id%||%"");if(!resolved$valid||!nrow(resolved$rows)){showNotification(resolved$message%||%"Select at least one valid comparison.",type="error");return(FALSE)};id_error<-validate_fitness_analysis_id(aid);if(nzchar(id_error)){showNotification(id_error,type="error");return(FALSE)};missing_inputs<-resolved$rows$status%in%"Missing input";if(any(missing_inputs)){showNotification("Selected comparisons are missing required raw SummaryTable inputs.",type="error");return(FALSE)};selected<-unname(as.character(resolved$rows$comparison_id));mode<-isolate(input$fitness_execution_mode%||%"preview");annotate<-isTRUE(isolate(input$fitness_annotate_classifier));force<-FALSE;all_valid<-setequal(selected,as.character(fitness_design()$comparison_id[fitness_design()$valid]));keep<-resolve_keep_going("fitness_analysis",length(selected),all_valid,TRUE);args<-c(file.path(repo_root,"scripts/local/ytab_run_treated_vs_parent.py"),"--project-config",active_project_path(),"--analysis-id",aid,"--input-mode","raw-summary","--comparisons",paste(selected,collapse=","),if(annotate)c("--annotate-classifier","--classifier-target",isolate(input$fitness_classifier_target)),if(mode=="preview")"--dry-run",if(force)"--force",if(keep)"--keep-going");metadata<-list(selected_comparisons=selected,selected_comparison_count=length(selected),comparison_set_label=if(all_valid)"All valid comparisons"else"Selected subset",execution_mode=mode,classifier_annotation=annotate,force=force,keep_going=keep);jobs$start_job(python_bin(),args,active_project_path(),"fitness_analysis",wd=repo_root,selected_items=selected,metadata=metadata);TRUE}
  observeEvent(input$fitness_run,{mode<-input$fitness_execution_mode%||%"preview";selected<-fitness_selected_comparisons();req(length(selected));if(mode=="preview"||!is.null(fitness_matching())){start_fitness_job();return()};showModal(modalDialog(title="Run fitness analysis?",tags$dl(class="ytab-meta",tags$dt("Analysis"),tags$dd(gsub("_"," ",input$fitness_analysis_id)),tags$dt("Comparisons"),tags$dd(length(selected)),tags$dt("Input"),tags$dd("Raw SummaryTable"),tags$dt("Normalization"),tags$dd("CPM inside R"),tags$dt("Classifier annotation"),tags$dd(if(isTRUE(input$fitness_annotate_classifier))"Yes"else"No"),tags$dt("Matching result"),tags$dd(if(is.null(fitness_matching()))"No"else"Yes")),tags$details(tags$summary("Selected comparisons"),tags$p(paste(selected,collapse=", "))),footer=tagList(modalButton("Cancel"),actionButton("confirm_fitness_run","Run fitness analysis",class="btn-primary"))))},ignoreInit=TRUE)
  observeEvent(input$confirm_fitness_run,{removeModal();start_fitness_job()},ignoreInit=TRUE)
  current_fitness_result<-reactive({results<-fitness_results();if(!length(results))return(NULL);hit<-Filter(function(x)identical(x$run_id,fitness_selected_result())||identical(x$analysis_id,fitness_selected_result()),results);if(length(hit))hit[[1]]else fitness_result_state()$matching_result%||%choose_latest_fitness_result(results)})
  output$fitness_result_selector<-renderUI({results<-fitness_results();if(length(results)<=1L)return(NULL);values<-unname(vapply(results,`[[`,"","run_id"));labels<-unname(vapply(results,fitness_result_display_name,""));selectInput("fitness_result_choice","Fitness result",choices=as.list(setNames(values,labels)),selected=current_fitness_result()$run_id%||%"")})
  observeEvent(input$fitness_result_choice,fitness_selected_result(input$fitness_result_choice%||%""),ignoreInit=TRUE)
  output$fitness_results_state<-renderUI({
    state<-fitness_result_state();result<-current_fitness_result()
    if(state$status=="no_design")return(tags$p(class="ytab-warning","Generate a valid comparison design before running fitness analysis."))
    if(state$status=="no_selection")return(tags$p(class="ytab-warning","Select one or more valid comparisons."))
    if(is.null(result))return(tags$p(class="text-muted","No fitness result matches the current comparison selection."))
    tagList(
      tags$div(
        class=if(result$legacy)"ytab-result-card ytab-result-historical"else"ytab-result-card ytab-result-matching",
        title=paste("Analysis ID:",result$analysis_id,"\nRun ID:",result$run_id),
        tags$div(class="ytab-result-heading",
                 tags$h4(fitness_result_display_name(result)),
                 tags$span(class="ytab-result-badge",if(result$legacy)"Historical result"else"Matching result")),
        tags$dl(class="ytab-meta",
                tags$dt("Comparisons"),tags$dd(length(result$selected_comparisons)),
                tags$dt("Features"),tags$dd(result$feature_count),
                tags$dt("Classifier annotation"),tags$dd(if(result$classifier_annotation)"Yes"else"No"))
      ),
      if(result$legacy)tags$p(class="text-muted","This historical result lacks the complete design-aware cache signature and is not claimed as a current match.")
    )
  })
  output$fitness_smoke_warning<-renderUI(if(fitness_smoke_project(active()))tags$div(class="alert alert-warning","This project uses a shallow FASTQ subset for software validation. Its fitness calls are not final biological conclusions."))
  fitness_data<-reactive(fitness_result_data(current_fitness_result()));fitness_filtered<-reactive(filter_fitness_results(fitness_data(),input$fitness_call_filter%||%"All",input$fitness_full_gene_search%||%""))
  observe({data<-fitness_data();column<-fitness_call_column(data);calls<-if(nzchar(column))sort(unique(as.character(data[[column]])))else character();calls<-calls[!is.na(calls)&nzchar(calls)];choices<-as.list(c("All",unname(calls)));names(choices)<-c("All calls",tools::toTitleCase(gsub("_"," ",calls)));updateSelectInput(session,"fitness_call_filter",choices=choices)})
  output$fitness_summary_cards<-renderUI({data<-fitness_data();if(!nrow(data))return(NULL);counts<-fitness_call_counts(data);card<-function(key,label)tags$div(tags$b(counts[[key]]),label);tags$div(class="ytab-stat-grid",tags$div(tags$b(nrow(data)),"Total features"),card("consistently_depleted","Consistently depleted"),card("consistently_enriched","Consistently enriched"),card("single_pool_depleted","Single-pool depleted"),card("single_pool_enriched","Single-pool enriched"),card("mixed","Mixed"))})
  output$fitness_filtered_count<-renderUI(if(nrow(fitness_data()))tags$p(sprintf("Showing %d of %d features.",nrow(fitness_filtered()),nrow(fitness_data()))))
  output$fitness_results_table<-DT::renderDT({x<-fitness_result_table_data(fitness_filtered(),repo_root);if(!nrow(x))return(NULL);ytab_gene_details_datatable(x,filter="top",options=list(pageLength=10,lengthMenu=c(10,25,50),autoWidth=FALSE,ordering=TRUE,searching=TRUE,scrollX=TRUE,columnDefs=list(list(className="dt-right",targets=which(vapply(x,is.numeric,FALSE))-1L),list(className="ytab-nowrap",targets=0))))})
  fitness_ma_selected_hits <- reactive({
    result <- current_fitness_result(); if (is.null(result)) return(data.frame())
    mode <- input$fitness_ma_mode %||% "combined"; pair <- input$fitness_ma_pair %||% ""
    direction <- input$fitness_ma_hit_direction %||% "both"; n <- as.integer(input$fitness_ma_top_n %||% 10)
    min_support <- if (identical(mode, "combined")) as.integer(input$fitness_ma_min_support %||% 0) else 0L
    x <- fitness_ma_top_hits_table(result, mode, pair, direction, n, input$fitness_ma_annotation_mode %||% "top", input$fitness_ma_custom_features %||% "", min_support)
    if (nrow(x) && "feature_id" %in% names(x)) x <- ytab_join_glabrata_display(x, repo_root, "feature_id")
    x
  })
  fitness_ma_selected_hits_visible <- reactive({
    x <- fitness_ma_selected_hits()
    q <- trimws(input$fitness_gene_search %||% "")
    if (!nrow(x) || !nzchar(q)) return(x)
    fields <- intersect(c("feature_id", "original_feature_id", "cagl_display_id", "gene_display_name", "label", "hit_direction", "supporting_pool_ids"), names(x))
    if (!length(fields)) return(x)
    keep <- Reduce(`|`, lapply(x[fields], function(v) grepl(q, as.character(v), ignore.case = TRUE)))
    x[keep, , drop = FALSE]
  })
  output$fitness_selected_top_hits_table <- DT::renderDT({
    x <- fitness_ma_selected_hits_visible()
    if (!nrow(x)) return(NULL)
    details <- ytab_glabrata_gene_detail_columns(x)
    keep <- intersect(c("selected_rank", "cagl_display_id", "gene_display_name", "cg_to_sc_relationship_display", "hit_direction", "pair", "log2fc", "mean_abundance", "rank_z_strength", "supporting_pool_ids"), names(x))
    x <- x[, keep, drop = FALSE]
    for (nm in intersect(c("log2fc", "mean_abundance", "rank_lfc_strength", "rank_cpm_support", "rank_z_strength"), names(x))) x[[nm]] <- round(as.numeric(x[[nm]]), 3)
    for (nm in intersect(c("candidate_rank_order", "valid_pool_n", "rank_support_n", "selected_min_support_pools"), names(x))) x[[nm]] <- as.integer(x[[nm]])
    names(x)[names(x) == "selected_rank"] <- "Rank"; names(x)[names(x) == "cagl_display_id"] <- "CAGL ID"; names(x)[names(x) == "gene_display_name"] <- "Gene name"; names(x)[names(x) == "cg_to_sc_relationship_display"] <- "Cg-to-Sc relationship"; names(x)[names(x) == "hit_direction"] <- "Hit direction"; names(x)[names(x) == "log2fc"] <- "Directional log2FC"; names(x)[names(x) == "mean_abundance"] <- "CPM/read support"; names(x)[names(x) == "rank_z_strength"] <- "Local z-score support"; names(x)[names(x) == "supporting_pool_ids"] <- "Supporting pool IDs"
    x <- cbind(x, details)
    ytab_gene_details_datatable(x, class = "compact stripe hover", order = list(list(0, "asc")), options = list(pageLength = 10, lengthMenu = c(10, 25, 50, 100), autoWidth = FALSE, ordering = TRUE, searching = TRUE, scrollX = TRUE, columnDefs = list(list(className = "dt-right", targets = which(vapply(x, is.numeric, FALSE)) - 1L), list(className = "ytab-nowrap", targets = 0))))
  })
  output$fitness_visualization_selector<-renderUI({
    result<-current_fitness_result()
    plots<-fitness_generated_plot_inventory_for_result(active(),result)
    result_ok<-!is.null(result)&&nrow(fitness_ma_data(result,"combined"))
    control_z_choices<-fitness_control_z_scope_choices(result)
    library_size_choices<-fitness_library_size_scope_choices(result)
    mean_log2fc_ok<-nrow(fitness_mean_log2fc_data(result))>0L
    choices<-c(if(result_ok)c("MA plot"="combined_ma","Condition versus control log-log scatter"="condition_control","Top selected hits log2FC heatmap"="selected_hit_heatmap"),if(length(library_size_choices))c("Feature-level library sizes"="library_sizes"),if(mean_log2fc_ok)c("Ranked mean log2FC"="ranked_mean_log2fc","Mean log2FC distribution"="mean_log2fc_distribution"),if(length(control_z_choices))c("Control-control z histogram"="control_z"),"Effect size versus z-score"="effect_size")
    if(nrow(plots))choices<-c(choices,setNames(paste0("generated:",seq_len(nrow(plots))),tools::file_path_sans_ext(gsub("_", " ", plots$filename))))
    selectInput("fitness_visualization_choice","Plot",choices=choices,selected=if(result_ok)"combined_ma"else"effect_size")
  })
  output$fitness_plot_controls<-renderUI({
    choice<-input$fitness_visualization_choice%||%"effect_size"
    point_choices<-c("combined_ma","condition_control","ranked_mean_log2fc","effect_size")
    ytab_plot_customization_controls("fitness",
                                     include_points = choice %in% point_choices,
                                     include_bars = identical(choice,"library_sizes"),
                                     include_value_labels = identical(choice,"library_sizes"),
                                     include_labels = !startsWith(choice,"generated:"),
                                     default_height = "large",
                                     default_show_value_labels = TRUE)
  })
  output$fitness_ma_controls<-renderUI({
    if(!((input$fitness_visualization_choice%||%"") %in% c("combined_ma","condition_control","selected_hit_heatmap")))return(NULL)
    result<-current_fitness_result();pairs<-fitness_ma_pair_choices(result)
    tagList(
      selectInput("fitness_ma_mode", "Comparison view", choices=c("Combined across pools"="combined","Individual treated-control pair"="individual"), selected=input$fitness_ma_mode%||%"combined"),
      uiOutput("fitness_ma_support_control"),
      conditionalPanel("input.fitness_ma_mode == 'individual'", selectInput("fitness_ma_pair", "Pair", choices=if(length(pairs))pairs else "No pair data available")),
      selectInput("fitness_ma_hit_direction", "Hit direction", choices=c("Depleted in treatment"="depleted","Enriched in treatment"="enriched","Both depleted and enriched"="both","None"="none"), selected=input$fitness_ma_hit_direction%||%"both"),
      numericInput("fitness_ma_top_n", tagList("Number of hits",tags$span(title="Top hits are ranked by directional log2FC magnitude and CPM/read support, with local z-score used as statistical support."," ⓘ")), value=10, min=0, max=100, step=1),
      selectInput("fitness_ma_annotation_mode", "Annotation mode", choices=c("None"="none","Top highlighted hits"="top","Custom feature list"="custom","Top highlighted hits + custom feature list"="top_custom"), selected=input$fitness_ma_annotation_mode%||%"top"),
      conditionalPanel("input.fitness_ma_annotation_mode == 'custom' || input.fitness_ma_annotation_mode == 'top_custom'", textAreaInput("fitness_ma_custom_features", "Custom features", value="", placeholder="Feature or gene IDs, separated by commas or new lines", rows=3)),
      uiOutput("fitness_ma_selection_status")
    )
  })
  output$fitness_current_plot_downloads<-renderUI({
    result<-current_fitness_result();if(is.null(result))return(NULL)
    choice<-input$fitness_visualization_choice%||%"effect_size"
    tagList(
      tags$details(class="ytab-more-options", open=FALSE,
        tags$summary("Current plot downloads"),
        tags$div(class="ytab-actions",
          downloadButton("download_fitness_ma_plot", "Download current plot"),
          downloadButton("download_fitness_ma_data", "Download current plotted data")
        ),
        if(choice%in%c("combined_ma","condition_control","selected_hit_heatmap"))
          tagList(tags$hr(),tags$p("Download selected top hits"),downloadButton("download_fitness_ma_top_hits", "Download selected top hits CSV"))
      )
    )
  })
  output$fitness_library_size_controls<-renderUI({
    if(!identical(input$fitness_visualization_choice%||%"","library_sizes"))return(NULL)
    result<-current_fitness_result();choices<-fitness_library_size_scope_choices(result)
    if(!length(choices))return(tags$p(class="text-muted","Feature-level library-size data are unavailable for this result."))
    selected<-input$fitness_library_size_scope%||%"combined";if(!(selected%in%unlist(choices,use.names=FALSE)))selected<-"combined"
    tagList(
      selectInput("fitness_library_size_scope","Library-size scope",choices=choices,selected=selected),
      tags$p(class="text-muted","Default combines backgrounds; individual backgrounds remain available for background-specific library-size checks.")
    )
  })
  output$fitness_control_z_controls<-renderUI({
    if(!identical(input$fitness_visualization_choice%||%"","control_z"))return(NULL)
    result<-current_fitness_result();choices<-fitness_control_z_scope_choices(result)
    if(!length(choices))return(tags$p(class="text-muted","Control-control z-score data are unavailable for this result."))
    selected<-input$fitness_control_z_scope%||%"combined";if(!(selected%in%unlist(choices,use.names=FALSE)))selected<-"combined"
    tagList(
      selectInput("fitness_control_z_scope","Histogram scope",choices=choices,selected=selected),
      tags$p(class="text-muted","Default combines backgrounds; individual backgrounds remain available for parent-parent noise-model diagnostics.")
    )
  })
  output$fitness_ma_support_control<-renderUI({
    if(!identical(input$fitness_ma_mode%||%"combined", "combined"))return(NULL)
    result<-current_fitness_result();pools<-length(fitness_ma_available_pools(result))
    if(!pools)return(tags$p(class="text-muted","Minimum supporting pools is unavailable because by-pool z-score data are missing or incomplete."))
    choices<-c("All support classes"="0",setNames(as.character(seq_len(pools)),paste0(">= ",seq_len(pools)," supporting pool",ifelse(seq_len(pools)==1,"","s"))))
    selected<-as.character(input$fitness_ma_min_support%||%"0");if(!selected%in%unname(choices))selected<-"0"
    selectInput("fitness_ma_min_support",tagList("Concordance filter",tags$span(title="Optionally filters top-hit selection by how many valid pools show evidence in the selected direction."," ⓘ")),choices=choices,selected=selected)
  })
  output$fitness_ma_selection_status<-renderUI({
    if(!identical(input$fitness_visualization_choice%||%"", "combined_ma"))return(NULL)
    result<-current_fitness_result();if(is.null(result))return(NULL)
    mode<-input$fitness_ma_mode%||%"combined";direction<-input$fitness_ma_hit_direction%||%"both"
    threshold<-if(identical(mode,"combined"))as.integer(input$fitness_ma_min_support%||%0)else 0L
    n<-as.integer(input$fitness_ma_top_n%||%10);if(is.na(n)||n<0)n<-0L
    ranked<-fitness_ma_rank_data(result,mode,input$fitness_ma_pair%||%"",direction,threshold)
    if(!nrow(ranked))return(tags$p(class="text-muted","No MA candidates are available for the selected comparison."))
    ranked_total<-sum(!is.na(ranked$overall_candidate_rank))
    passing<-sum(!is.na(ranked$support_pass)&ranked$support_pass)
    highlighted<-fitness_ma_highlight_rows(ranked,direction,n,threshold)
    labels<-fitness_ma_label_table(ranked,highlighted,input$fitness_ma_annotation_mode%||%"top",input$fitness_ma_custom_features%||%"",mode,input$fitness_ma_pair%||%"")
    if(identical(direction,"none")||n==0L) return(tags$p(class="text-muted","Automatic hit highlighting is off; the MA plot shows the neutral feature cloud."))
    if(!passing)return(tags$p(class="text-muted",if(threshold>0L)sprintf("Ranked candidates passing log2FC/CPM filters: %d; Passing concordance filter: 0.",ranked_total)else"Ranked candidates passing log2FC/CPM filters: 0."))
    if(threshold>0L) tags$p(class="text-muted",sprintf("Ranked candidates passing log2FC/CPM filters: %d; Passing concordance filter: %d; Selected top hits: %d; Labeled: %d.",ranked_total,passing,length(highlighted),nrow(labels$data))) else tags$p(class="text-muted",sprintf("Ranked candidates passing log2FC/CPM filters: %d; Selected top hits: %d; Labeled: %d.",ranked_total,length(highlighted),nrow(labels$data)))
  })
  output$fitness_selected_top_hits_empty <- renderUI({
    if (nrow(fitness_ma_selected_hits_visible())) return(NULL)
    if (nrow(fitness_ma_selected_hits())) return(tags$p(class = "text-muted", "No selected top hits match the current search."))
    tags$p(class = "text-muted", "No selected top hits pass the current log2FC, CPM, and concordance filters.")
  })
  output$fitness_selected_top_hits_count <- renderUI({
    total <- nrow(fitness_ma_selected_hits()); visible <- nrow(fitness_ma_selected_hits_visible())
    if (!total) return(NULL)
    if (nzchar(trimws(input$fitness_gene_search %||% ""))) tags$p(class = "text-muted", sprintf("Showing %d of %d selected top hits matching search.", visible, total))
    else tags$p(class = "text-muted", sprintf("Showing %d selected top hit%s.", total, if (total == 1L) "" else "s"))
  })
  output$fitness_selected_visualization<-renderUI({
    choice<-input$fitness_visualization_choice%||%"effect_size"
    fitness_height<-ytab_plot_height_px(input$fitness_plot_height%||%"large",default="large")
    if(identical(choice,"combined_ma"))return(ytab_plot_frame(plotOutput("fitness_ma_plot",width="100%",height=fitness_height),input$fitness_plot_width%||%"standard","app-rendered"))
    if(identical(choice,"selected_hit_heatmap"))return(ytab_plot_frame(plotOutput("fitness_selected_hit_heatmap",width="100%",height=fitness_height),input$fitness_plot_width%||%"standard","app-rendered"))
    if(identical(choice,"library_sizes"))return(ytab_plot_frame(plotOutput("fitness_library_size_plot",width="100%",height=fitness_height),input$fitness_plot_width%||%"standard","app-rendered"))
    if(identical(choice,"ranked_mean_log2fc"))return(ytab_plot_frame(plotOutput("fitness_ranked_mean_log2fc_plot",width="100%",height=fitness_height),input$fitness_plot_width%||%"standard","app-rendered"))
    if(identical(choice,"mean_log2fc_distribution"))return(ytab_plot_frame(plotOutput("fitness_mean_log2fc_distribution_plot",width="100%",height=fitness_height),input$fitness_plot_width%||%"standard","app-rendered"))
    if(identical(choice,"control_z"))return(ytab_plot_frame(plotOutput("fitness_control_z_plot",width="100%",height=fitness_height),input$fitness_plot_width%||%"standard","app-rendered"))
    if(startsWith(choice,"generated:")){
      plots<-fitness_generated_plot_inventory_for_result(active(),current_fitness_result());index<-suppressWarnings(as.integer(sub("^generated:","",choice)))
      if(is.na(index)||index<1L||index>nrow(plots))return(tags$p(class="text-muted","Selected plot is unavailable."))
      return(tags$article(class="ytab-static-image-card ytab-release-card",tags$div(class="ytab-plot-card-header",tags$h4(plots$title[[index]])),tags$div(class="ytab-diagnostic-preview",tags$img(src=plots$served_url[[index]],alt=plots$filename[[index]],loading="lazy",style=paste0("width:100%;max-height:",fitness_height,";object-fit:contain"))),tags$details(tags$summary("Show filename"),tags$code(plots$filename[[index]]))))
    }
    if(identical(choice,"condition_control"))return(ytab_plot_frame(plotOutput("fitness_condition_control_plot",width="100%",height=fitness_height),input$fitness_plot_width%||%"standard","app-rendered"))
    if(identical(choice,"effect_size"))return(tagList(uiOutput("fitness_effect_plot_state"),ytab_plot_frame(plotOutput("fitness_effect_plot",width="100%",height=fitness_height),input$fitness_plot_width%||%"standard","app-rendered")))
    tags$p(class="text-muted","Selected plot is unavailable.")
  })
  output$fitness_ma_plot<-renderPlot(ytab_with_plot_display_options(input,"fitness",{
    result<-current_fitness_result();req(!is.null(result));mode<-input$fitness_ma_mode%||%"combined";pair<-input$fitness_ma_pair%||%"";direction<-input$fitness_ma_hit_direction%||%"both";n<-suppressWarnings(as.integer(input$fitness_ma_top_n%||%10));if(is.na(n))n<-10
    min_support<-if(identical(mode,"combined"))as.integer(input$fitness_ma_min_support%||%0)else 0L
    plot_fitness_ma(result,mode,pair,direction,n,input$fitness_ma_annotation_mode%||%"top",input$fitness_ma_custom_features%||%"",input$fitness_text_size%||%"medium",qc_plot_grid_enabled(),as.numeric(input$fitness_point_size%||%1.5),min_support,repo_root)
  }))
  output$fitness_condition_control_plot<-renderPlot(ytab_with_plot_display_options(input,"fitness",{
    result<-current_fitness_result();req(!is.null(result));mode<-input$fitness_ma_mode%||%"combined";pair<-input$fitness_ma_pair%||%"";direction<-input$fitness_ma_hit_direction%||%"both";n<-suppressWarnings(as.integer(input$fitness_ma_top_n%||%10));if(is.na(n))n<-10
    min_support<-if(identical(mode,"combined"))as.integer(input$fitness_ma_min_support%||%0)else 0L
    plot_fitness_condition_control_scatter(result,mode,pair,direction,n,input$fitness_ma_annotation_mode%||%"top",input$fitness_ma_custom_features%||%"",input$fitness_text_size%||%"medium",qc_plot_grid_enabled(),as.numeric(input$fitness_point_size%||%1.5),min_support,repo_root)
  }))
  output$fitness_control_z_plot<-renderPlot(ytab_with_plot_display_options(input,"fitness",{
    result<-current_fitness_result();req(!is.null(result))
    plot_fitness_control_control_z_histogram(result,input$fitness_control_z_scope%||%"combined")
  }))
  output$fitness_library_size_plot<-renderPlot(ytab_with_plot_display_options(input,"fitness",{
    result<-current_fitness_result();req(!is.null(result))
    plot_fitness_library_sizes_feature_reads(result,input$fitness_library_size_scope%||%"combined")
  }))
  output$fitness_ranked_mean_log2fc_plot<-renderPlot(ytab_with_plot_display_options(input,"fitness",{
    result<-current_fitness_result();req(!is.null(result))
    plot_fitness_ranked_mean_log2fc(result)
  }))
  output$fitness_mean_log2fc_distribution_plot<-renderPlot(ytab_with_plot_display_options(input,"fitness",{
    result<-current_fitness_result();req(!is.null(result))
    plot_fitness_mean_log2fc_distribution(result)
  }))
  output$fitness_selected_hit_heatmap<-renderPlot(ytab_with_plot_display_options(input,"fitness",{
    result<-current_fitness_result();req(!is.null(result));mode<-input$fitness_ma_mode%||%"combined";pair<-input$fitness_ma_pair%||%""
    plot_fitness_selected_hit_heatmap(result,fitness_ma_selected_hits(),mode,pair,input$fitness_text_size%||%"medium")
  }))
  fitness_effect_data<-reactive(normalize_fitness_result_columns(fitness_data())$data)
  fitness_effect_columns<-reactive({data<-fitness_effect_data();effect<-intersect(c("mean_log2fc"),names(data));z<-intersect(c("max_z","min_z"),names(data));list(effect=if(length(effect))effect[[1]]else"",z=if(length(z))z[[1]]else"")})
  output$fitness_effect_plot_state<-renderUI({columns<-fitness_effect_columns();if(!nzchar(columns$effect)||!nzchar(columns$z))tags$p(class="text-muted","This result does not contain the columns required for this plot.")})
  output$fitness_effect_plot<-renderPlot(ytab_with_plot_display_options(input,"fitness",{data<-fitness_effect_data();columns<-fitness_effect_columns();req(nrow(data),nzchar(columns$effect),nzchar(columns$z));old<-qc_plot_par(mar=c(5,5,3,1));on.exit(par(old),add=TRUE);point_opacity<-suppressWarnings(as.numeric(input$fitness_point_opacity%||%0.55));if(is.na(point_opacity))point_opacity<-0.55;point_size<-suppressWarnings(as.numeric(input$fitness_point_size%||%1.5));if(is.na(point_size))point_size<-1.5;plot(as.numeric(data[[columns$effect]]),as.numeric(data[[columns$z]]),pch=16,cex=point_size,col=rgb(0,0,0,max(0.1,min(1,point_opacity))),xlab="Mean log2 fold change",ylab=gsub("_"," ",columns$z),main="Effect size versus z-score")}))
  output$fitness_comparison_results<-DT::renderDT({result<-current_fitness_result();if(is.null(result)||!file.exists(result$comparison_summary))return(NULL);x<-read.csv(result$comparison_summary,stringsAsFactors=FALSE,check.names=FALSE);compact_qc_table(x)})
  output$fitness_result_technical<-renderUI({result<-current_fitness_result();if(is.null(result))return(NULL);tags$div(class="ytab-technical-console",tags$p("Result: ",tags$code(relative_project_path(result$table,active()$project_root))),tags$p("Manifest: ",tags$code(if(file.exists(result$manifest))relative_project_path(result$manifest,active()$project_root)else"Legacy manifest unavailable")),tags$p("Input mode: raw-summary"),tags$p("Normalization: CPM inside R"))})
  fitness_download_path<-function(path,message="No fitness result is available for download."){if(!nzchar(path%||%"")||!file.exists(path)){showNotification(message,type="warning");validate(need(FALSE,message))};path}
  output$download_fitness_full<-downloadHandler(filename=function()"treated_vs_parent_results.csv",content=function(file){result<-current_fitness_result();path<-if(is.null(result))""else result$table;file.copy(fitness_download_path(path),file,overwrite=TRUE)})
  output$download_fitness_filtered<-downloadHandler(filename=function()"treated_vs_parent_results.filtered.csv",content=function(file)write.csv(fitness_filtered(),file,row.names=FALSE))
  output$download_fitness_design<-downloadHandler(filename=function()"comparison_design.csv",content=function(file)file.copy(fitness_download_path(fitness_design_path(active()),"No comparison design is available for download."),file,overwrite=TRUE))
  output$download_fitness_comparison<-downloadHandler(filename=function()"treated_vs_parent_comparison_summary.csv",content=function(file){result<-current_fitness_result();path<-if(is.null(result))""else result$comparison_summary;file.copy(fitness_download_path(path,"No comparison-level fitness summary is available for download."),file,overwrite=TRUE)})
  output$download_fitness_manifest<-downloadHandler(filename=function()"treated_vs_parent_manifest.json",content=function(file){result<-current_fitness_result();path<-if(is.null(result))""else result$manifest;file.copy(fitness_download_path(path,"No fitness run manifest is available for download."),file,overwrite=TRUE)})
  fitness_current_plot_slug<-function(){
    choice<-input$fitness_visualization_choice%||%"effect_size"
    mode<-input$fitness_ma_mode%||%"combined"
    pair<-input$fitness_ma_pair%||%""
    scope<-if(identical(choice,"library_sizes"))input$fitness_library_size_scope%||%"combined" else if(identical(choice,"control_z"))input$fitness_control_z_scope%||%"combined" else ""
    qc_download_slug(active()$project_id,"fitness",choice,mode,pair,scope,input$fitness_ma_hit_direction%||%"")
  }
  fitness_current_ma_state<-function(){
    mode<-input$fitness_ma_mode%||%"combined"
    list(
      mode=mode,
      pair=input$fitness_ma_pair%||%"",
      direction=input$fitness_ma_hit_direction%||%"both",
      n=as.integer(input$fitness_ma_top_n%||%10),
      annotation=input$fitness_ma_annotation_mode%||%"top",
      custom=input$fitness_ma_custom_features%||%"",
      min_support=if(identical(mode,"combined"))as.integer(input$fitness_ma_min_support%||%0)else 0L
    )
  }
  fitness_current_plot_data<-function(){
    result<-current_fitness_result();req(!is.null(result))
    choice<-input$fitness_visualization_choice%||%"effect_size"
    state<-fitness_current_ma_state()
    if(identical(choice,"combined_ma")){
      data<-fitness_ma_rank_data(result,state$mode,state$pair,state$direction,state$min_support)
      hi<-fitness_ma_highlight_rows(data,state$direction,state$n,state$min_support)
      ann<-fitness_ma_annotation_rows(data,hi,state$annotation,state$custom)
      highlighted<-seq_len(nrow(data))%in%hi
      data$plotted<-TRUE;data$highlighted<-highlighted;data$labeled<-seq_len(nrow(data))%in%ann$rows
      data$annotation_source<-ifelse(data$labeled&highlighted,"top_hit",ifelse(data$labeled,"custom",ifelse(highlighted,"top_hit","none")))
      data$comparison_view<-if(identical(state$mode,"individual"))"individual_pair"else"combined";data$pair<-if(nzchar(state$pair))state$pair else NA_character_
      return(data)
    }
    if(identical(choice,"condition_control")){
      data<-fitness_condition_control_plot_data(result,state$mode,state$pair)
      ranked<-fitness_ma_rank_data(result,state$mode,state$pair,state$direction,state$min_support)
      hi<-fitness_ma_highlight_rows(ranked,state$direction,state$n,state$min_support)
      ann<-fitness_ma_annotation_rows(ranked,hi,state$annotation,state$custom,state$mode,state$pair)
      data$x_log10_control_plus_1<-log10(qc_plot_numeric(data$control_abundance)+1)
      data$y_log10_treated_plus_1<-log10(qc_plot_numeric(data$treated_abundance)+1)
      data$highlighted<-as.character(data$feature_id)%in%as.character(ranked$feature_id[hi])
      data$labeled<-as.character(data$feature_id)%in%as.character(ranked$feature_id[ann$rows])
      data$comparison_view<-if(identical(state$mode,"individual"))"individual_pair"else"combined";data$pair<-if(nzchar(state$pair))state$pair else NA_character_
      return(data)
    }
    if(identical(choice,"selected_hit_heatmap")){
      hm<-fitness_selected_hit_heatmap_data(result,fitness_ma_selected_hits(),state$mode,state$pair)
      return(hm$data)
    }
    if(identical(choice,"library_sizes")){
      data<-fitness_library_size_data(result);scope<-input$fitness_library_size_scope%||%"combined"
      if(!identical(scope,"combined"))data<-data[as.character(data$background)==scope,,drop=FALSE]
      data$plot_scope<-scope
      return(data)
    }
    if(identical(choice,"ranked_mean_log2fc")){
      data<-fitness_mean_log2fc_data(result)
      data$plot_type<-"ranked_mean_log2FC"
      return(data)
    }
    if(identical(choice,"mean_log2fc_distribution")){
      data<-fitness_mean_log2fc_data(result)
      data$plot_type<-"mean_log2FC_distribution"
      return(data)
    }
    if(identical(choice,"control_z")){
      data<-fitness_control_z_data(result);scope<-input$fitness_control_z_scope%||%"combined"
      if(!identical(scope,"combined"))data<-data[as.character(data$background)==scope,,drop=FALSE]
      data$plot_scope<-scope
      return(data)
    }
    if(startsWith(choice,"generated:")){
      plots<-fitness_generated_plot_inventory_for_result(active(),result);index<-suppressWarnings(as.integer(sub("^generated:","",choice)))
      return(data.frame(message="Current selected plot is a static generated PNG; plotted source data are not available through this button.", file=if(!is.na(index)&&index>=1L&&index<=nrow(plots))plots$file[[index]]else NA_character_, stringsAsFactors=FALSE))
    }
    data<-fitness_effect_data()
    columns<-fitness_effect_columns()
    data$effect_column<-columns$effect;data$z_column<-columns$z
    data
  }
  fitness_render_current_plot<-function(file){
    result<-current_fitness_result();req(!is.null(result))
    choice<-input$fitness_visualization_choice%||%"effect_size"
    if(startsWith(choice,"generated:")){
      plots<-fitness_generated_plot_inventory_for_result(active(),result);index<-suppressWarnings(as.integer(sub("^generated:","",choice)))
      validate(need(!is.na(index)&&index>=1L&&index<=nrow(plots)&&file.exists(plots$file[[index]]),"Selected generated plot is unavailable."))
      file.copy(plots$file[[index]],file,overwrite=TRUE);return(invisible())
    }
    dim<-qc_download_dimensions("fitness","large")
    grDevices::png(file,width=dim$width,height=dim$height,res=150)
    on.exit(grDevices::dev.off(),add=TRUE)
    ytab_with_plot_display_options(input,"fitness",{
      state<-fitness_current_ma_state()
      if(identical(choice,"combined_ma"))plot_fitness_ma(result,state$mode,state$pair,state$direction,state$n,state$annotation,state$custom,input$fitness_text_size%||%"medium",qc_plot_grid_enabled(),as.numeric(input$fitness_point_size%||%1.5),state$min_support,repo_root)
      else if(identical(choice,"condition_control"))plot_fitness_condition_control_scatter(result,state$mode,state$pair,state$direction,state$n,state$annotation,state$custom,input$fitness_text_size%||%"medium",qc_plot_grid_enabled(),as.numeric(input$fitness_point_size%||%1.5),state$min_support,repo_root)
      else if(identical(choice,"selected_hit_heatmap"))plot_fitness_selected_hit_heatmap(result,fitness_ma_selected_hits(),state$mode,state$pair,input$fitness_text_size%||%"medium")
      else if(identical(choice,"library_sizes"))plot_fitness_library_sizes_feature_reads(result,input$fitness_library_size_scope%||%"combined")
      else if(identical(choice,"ranked_mean_log2fc"))plot_fitness_ranked_mean_log2fc(result)
      else if(identical(choice,"mean_log2fc_distribution"))plot_fitness_mean_log2fc_distribution(result)
      else if(identical(choice,"control_z"))plot_fitness_control_control_z_histogram(result,input$fitness_control_z_scope%||%"combined")
      else {
        data<-fitness_effect_data();columns<-fitness_effect_columns();req(nrow(data),nzchar(columns$effect),nzchar(columns$z))
        old<-qc_plot_par(mar=c(5,5,3,1));on.exit(par(old),add=TRUE)
        point_opacity<-suppressWarnings(as.numeric(input$fitness_point_opacity%||%0.55));if(is.na(point_opacity))point_opacity<-0.55
        point_size<-suppressWarnings(as.numeric(input$fitness_point_size%||%1.5));if(is.na(point_size))point_size<-1.5
        plot(as.numeric(data[[columns$effect]]),as.numeric(data[[columns$z]]),pch=16,cex=point_size,col=rgb(0,0,0,max(0.1,min(1,point_opacity))),xlab="Mean log2 fold change",ylab=gsub("_"," ",columns$z),main="Effect size versus z-score")
      }
    })
  }
  output$download_fitness_ma_plot<-downloadHandler(filename=function(){paste0(fitness_current_plot_slug(),".png")},content=function(file)fitness_render_current_plot(file))
  output$download_fitness_ma_data <- downloadHandler(
    filename = function() paste0(fitness_current_plot_slug(), ".csv"),
    content = function(file) {
      data <- fitness_current_plot_data()
      if (!nrow(data)) data <- data.frame(no_plotted_data = "No plotted data are available for the current Fitness plot.", stringsAsFactors = FALSE)
      write.csv(data, file, row.names = FALSE)
    }
  )
  output$download_fitness_ma_top_hits <- downloadHandler(
    filename = function() paste0(active()$project_id, ".fitness.ma_top_selected_hits.csv"),
    content = function(file) {
      result <- current_fitness_result(); req(!is.null(result))
      mode <- input$fitness_ma_mode %||% "combined"; pair <- input$fitness_ma_pair %||% ""
      direction <- input$fitness_ma_hit_direction %||% "both"; n <- as.integer(input$fitness_ma_top_n %||% 10)
      min_support <- if (identical(mode, "combined")) as.integer(input$fitness_ma_min_support %||% 0) else 0L
      top <- fitness_ma_selected_hits()
      if (!nrow(top)) top <- data.frame(no_top_hits_selected = "No top hits selected under current controls", comparison_view = if (identical(mode, "individual")) "individual_pair" else "combined", hit_direction = direction, selected_min_support_pools = min_support, stringsAsFactors = FALSE)
      else {
        raw_annotation <- c("cagl_id", "gwk60_id_clean", "qng_id", "cgla_common_name",
                            "cg_to_sc_relationship", "pre_WGD_Ancestor", "scer_gene_id",
                            "scer_gene_name", "SGD_essentiality", "SGD_description",
                            "cgla_gene_name_from_deseq")
        top <- top[, setdiff(names(top), raw_annotation), drop = FALSE]
        top <- ytab_standardize_glabrata_annotation_names(top)
      }
      write.csv(top, file, row.names = FALSE)
    }
  )
  output$download_fitness_heatmap_data <- downloadHandler(
    filename = function() paste0(active()$project_id, ".fitness.selected_top_hits_log2fc_heatmap.csv"),
    content = function(file) {
      result <- current_fitness_result(); req(!is.null(result))
      mode <- input$fitness_ma_mode %||% "combined"; pair <- input$fitness_ma_pair %||% ""
      hm <- fitness_selected_hit_heatmap_data(result, fitness_ma_selected_hits(), mode, pair)
      data <- hm$data
      if (!nrow(data)) data <- data.frame(no_heatmap_data = "No selected top-hit heatmap data under current controls", stringsAsFactors = FALSE)
      write.csv(data, file, row.names = FALSE)
    }
  )

  launch_pipeline <- function(dry_override=NULL){if(jobs$job_is_running()){showNotification("A pipeline job is already running.",type="warning");return()};dry<-if(is.null(dry_override))FALSE else dry_override;args<-c(file.path(repo_root,"scripts/local/ytab_run_pipeline.py"),"--project-config",active_project_path(),"--profile",input$pipeline_profile,"--threads",as.character(input$pipeline_threads),"--target",input$pipeline_target,"--analysis-id",input$pipeline_analysis_id,"--print-plan",if(dry)"--dry-run","--keep-going");jobs$start_job(python_bin(),args,active_project_path(),input$pipeline_profile,wd=repo_root)}
  observeEvent(input$preview_pipeline,launch_pipeline(TRUE));observeEvent(input$run_pipeline,launch_pipeline(NULL));observeEvent(input$cancel_pipeline,jobs$cancel_job())
  pipeline_text <- reactive({invalidateLater(500,session);jobs$poll_job();exit<-jobs$job_exit_status();job_key<-paste(jobs$current_job()$job_id,exit,sep=":");if(!is.na(exit)&&!identical(completed_job_seen(),job_key)){completed_job_seen(job_key);if(!is.null(state$values$path)){try(run_cli("ytab_project_status.py",c("--show-next")),silent=TRUE);state$load(active_project_path(),repo_root);status_tick(status_tick()+1L);if(identical(jobs$stage(),"library_diagnostics"))diagnostic_inventory_state(build_diagnostic_file_inventory(active()));if(exit==0L)showNotification(paste(jobs$stage(),"job finished."),type="message")}};paste(c(paste("Stage:",jobs$stage()),paste("Command:",jobs$command()),sprintf("Elapsed: %.1f seconds",jobs$elapsed()),paste("Exit status:",if(is.na(exit))"running"else exit),jobs$read_job_stdout(),jobs$read_job_stderr()),collapse="\n")})
  output$pipeline_job_text<-renderText(pipeline_text());output$pipeline_job_text_logs<-renderText(pipeline_text());output$pipeline_job_text_preprocessing<-renderText(pipeline_text());output$pipeline_job_text_qc<-renderText(pipeline_text())
  output$library_diagnostics_technical<-renderText(pipeline_text())
  output$library_diagnostics_result<-renderUI({
    pipeline_text();state<-qc_result_state();r<-qc_resolution();runs<-available_diagnostic_runs();historical<-state$latest_historical_result
    result_card<-function(title,run,badge="Historical result",class="ytab-result-historical"){
      summary<-build_diagnostics_result_summary(run);tags$div(class=paste("ytab-result-card",class),tags$div(class="ytab-result-heading",tags$h4(title),tags$span(class="ytab-result-badge",badge)),tags$dl(class="ytab-meta",tags$dt("Sample set"),tags$dd(summary$sample_set),tags$dt("Samples analyzed"),tags$dd(if(is.na(summary$samples_analyzed))"Unknown"else summary$samples_analyzed),tags$dt("Output files"),tags$dd(summary$output_files),tags$dt("Completed"),tags$dd(summary$completed)),if(nzchar(summary$warning))tags$p(class="ytab-warning",summary$warning),tags$div(class="ytab-actions",actionButton("view_diagnostics_result","View results",class="btn-secondary"),actionButton("view_diagnostic_files","View diagnostic files",class="btn-secondary")))
    }
    if(!length(r$eligible))return(tagList(tags$div(class="ytab-result-card ytab-result-empty",tags$h4("No current diagnostics selection"),tags$p("Select one or more eligible samples to preview or run diagnostics."),tags$span(class="ytab-result-badge","No matching result")),if(!is.null(historical))result_card("Last completed diagnostics",historical,if(historical$is_legacy)"Legacy result"else"Historical result")))
    p<-if(identical(jobs$stage(),"library_diagnostics")&&!jobs$job_is_running())jobs$last_progress()else NULL;job_matches<-!is.null(p)&&setequal(as.character(jobs$current_job()$selected_items%||%character()),r$eligible)
    if(job_matches&&identical(p$status,"dry_run_success"))return(tags$div(class="ytab-result-card ytab-result-preview",tags$div(class="ytab-result-heading",tags$h4("Preview complete"),tags$span(class="ytab-result-badge","Preview")),tags$dl(class="ytab-meta",tags$dt("Samples validated"),tags$dd(length(r$eligible)),tags$dt("Diagnostics executed"),tags$dd("No")),tags$p("The command and inputs were validated. Library diagnostics were not executed.")))
    if(job_matches&&identical(p$status,"failed"))return(tags$div(class="ytab-result-card alert alert-danger",tags$h4("Library Diagnostics failed"),tags$p("Open Technical details for the command and error output.")))
    completed<-Filter(function(x)!x$is_legacy&&setequal(x$selected_samples,r$eligible),runs);completed<-choose_latest_diagnostics_result(completed)
    if(job_matches&&identical(p$status,"success")&&!grepl("Cached result reused",p$message%||%"",fixed=TRUE)&&!is.null(completed))return(result_card("Library Diagnostics complete",completed,"Complete","ytab-result-matching"))
    if(!is.null(state$matching_result))return(result_card("Matching cached diagnostics",state$matching_result,"Matching result · Cached","ytab-result-matching"))
    tagList(tags$div(class="ytab-result-card ytab-result-no-match",tags$div(class="ytab-result-heading",tags$h4("No matching diagnostics result"),tags$span(class="ytab-result-badge","No matching result")),tags$dl(class="ytab-meta",tags$dt("Current sample set"),tags$dd(r$label),tags$dt("Selected samples"),tags$dd(length(r$eligible))),tags$p(if(is.null(historical))"No completed diagnostics results are available for this project."else"The latest available diagnostics were generated from a different sample selection."),tags$div(class="ytab-actions",actionButton("result_preview_diagnostics","Preview diagnostics",class="btn-secondary"),actionButton("result_run_diagnostics","Run diagnostics",class="btn-primary"))),if(!is.null(historical))result_card("Latest historical result",historical,if(historical$is_legacy)"Legacy result"else"Historical result"))
  })
  observeEvent(input$view_diagnostic_files,{run<-qc_result_state()$matching_result%||%qc_result_state()$latest_historical_result;if(!is.null(run))diagnostic_selected_run(run$run_id);go_to("quality_control","qc_tabs","diagnostic_files")},ignoreInit=TRUE)
  observeEvent(input$view_diagnostics_result,{run<-qc_result_state()$matching_result%||%qc_result_state()$latest_historical_result;if(!is.null(run))diagnostic_selected_run(run$run_id);go_to("quality_control","qc_tabs","diagnostic_files")},ignoreInit=TRUE)
  observeEvent(input$result_preview_diagnostics,{updateRadioButtons(session,"library_diagnostics_mode",selected="preview")},ignoreInit=TRUE);observeEvent(input$result_run_diagnostics,{updateRadioButtons(session,"library_diagnostics_mode",selected="run")},ignoreInit=TRUE)
  observeEvent(input$refresh_project_status,{run_cli("ytab_project_status.py",c("--show-next"));state$load(active_project_path(),repo_root);showNotification("Project status refreshed.")})
  refresh_preprocessing<-function(){run_cli("ytab_project_status.py",c("--show-next"));state$load(active_project_path(),repo_root);status_tick(status_tick()+1L);showNotification("Preprocessing status refreshed.")}
  for(id in c("refresh_preprocessing_status","refresh_mapfastq_status","refresh_create_hit_file_status","refresh_summary_table_status"))local({key<-id;observeEvent(input[[key]],refresh_preprocessing(),ignoreInit=TRUE)})
  observeEvent(input$refresh_reference_status,refresh_preprocessing(),ignoreInit=TRUE)
  for(id in c("cancel_mapfastq","cancel_create_hit_file","cancel_summary_table"))local({key<-id;observeEvent(input[[key]],jobs$cancel_job(),ignoreInit=TRUE)})
  for(id in c("change_mapfastq_selection","change_create_hit_file_selection","change_summary_table_selection"))local({key<-id;observeEvent(input[[key]],go_to("preprocessing","preprocessing_tabs","samples"),ignoreInit=TRUE)})
  output$project_status_text<-renderText({p<-active();if(is.null(p$status))"status not yet calculated" else paste(vapply(p$status$stages,function(x)sprintf("%-24s %s",x$stage,x$status),""),collapse="\n")})
  observeEvent(input$build_project_report,run_cli("ytab_build_project_report.py","--force"));observeEvent(input$build_export_bundle,run_cli("ytab_export_project.py","--force"))
  observeEvent(input$refresh_exports,exports_log(paste("Project outputs:",active()$project_root,"\nExports:",active()$export_root)))
  output$exports_text<-renderText(exports_log());output$command_output<-renderText(log_text());output$import_command_output<-renderText(log_text())
  output$reference_status<-renderText({p<-active();sprintf("Selected species: %s\nReference status: %s\nBowtie2 index complete: %s",p$species,stage_status(p$status,"reference"),p$reference$bowtie2_index_complete %||% FALSE)})

}

shinyApp(ui, server)
