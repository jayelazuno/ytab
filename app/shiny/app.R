library(shiny)
library(bslib)

`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
app_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL) %||% file.path(getwd(), "app.R")
app_dir <- normalizePath(dirname(app_file), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(app_dir, "../.."), winslash = "/", mustWork = TRUE)
for (helper in c("project_discovery.R", "project_state.R", "process_helpers.R", "preprocessing_status.R", "sample_selector.R", "job_manager.R", "job_progress.R", "navigation.R", "ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R", "table_display_helpers.R", "ui_landing.R", "ui_preprocessing.R", "qc_result_state.R", "qc_plot_utils.R", "qc_mapping_stats_plot.R", "qc_summary_library_plots.R", "qc_library_diagnostics_plots.R", "ui_qc.R", "fitness_design_state.R", "fitness_result_state.R", "fitness_generated_plots.R", "fitness_condition_control_plot.R", "ui_fitness.R", "essentiality_targets.R", "essentiality_state.R", "essentiality_results.R", "essentiality_commands.R", "essentiality_generated_plots.R", "ui_essentiality.R", "essentiality_server.R", "gene_domain_explorer_state.R", "ui_gene_domain_explorer.R", "gene_domain_explorer_server.R", "comparative_resources.R", "comparative_project_data.R", "ui_comparative.R", "ui_workspace.R"))
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
    tags$link(rel="stylesheet", href="ytab-landing.css"), tags$link(rel="stylesheet", href="ytab-qc.css"), tags$link(rel="stylesheet",href="job-progress.css"), tags$link(rel="stylesheet",href="ytab_release_ui.css")),
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
    for(route in c("ytab-diagnostics","ytab-diagnostics-project","ytab-diagnostics-export","ytab-project-output"))
      suppressWarnings(try(removeResourcePath(route),silent=TRUE))
    roots<-diagnostic_resource_roots(active())
    for(route in names(roots))addResourcePath(route,roots[[route]])
    addResourcePath("ytab-project-output",active()$project_root)
  },ignoreInit=FALSE)
  sample_pipeline_status <- reactive({status_tick();build_sample_pipeline_status(active())})
  preprocessing_selected_samples <- reactiveVal(character())
  selection_project <- reactiveVal("")
  qc_selected_samples <- reactiveVal(character())
  qc_selection_project <- reactiveVal("")
  diagnostic_inventory_state <- reactiveVal(data.frame())
  diagnostic_gallery_page <- reactiveVal(1L)
  diagnostic_selected_run <- reactiveVal("All")
  fitness_selected_comparisons <- reactiveVal(character())
  fitness_selection_project <- reactiveVal("")
  fitness_design_tick <- reactiveVal(0L)
  fitness_selected_result <- reactiveVal("")
  observeEvent(active_project_path(),{
    path<-active_project_path();if(!identical(path,selection_project())){d<-active()$samples;inc<-if("include"%in%names(d))tolower(as.character(d$include))%in%c("true","1","yes")else rep(TRUE,nrow(d));preprocessing_selected_samples(as.character(d$sample[inc]));selection_project(path)}
  },ignoreInit=FALSE)
  observeEvent(active_project_path(),{
    path<-active_project_path();if(!identical(path,qc_selection_project())){eligible<-qc_sample_eligibility(active(),sample_pipeline_status());qc_selected_samples(eligible$Sample[eligible$Eligible=="Yes"]);qc_selection_project(path);diagnostic_inventory_state(build_diagnostic_file_inventory(active()));diagnostic_gallery_page(1L);diagnostic_selected_run("All")}
  },ignoreInit=FALSE)
  sample_selector_server("preprocessing_samples",sample_data=reactive(active()$samples),status_data=sample_pipeline_status,selected_state=preprocessing_selected_samples)
  job_progress_server("global_job",jobs,active,compact=TRUE)
  job_progress_server("mapping_job",jobs,active,compact=FALSE)
  job_progress_server("library_diagnostics_job",jobs,active,compact=FALSE)
  job_progress_server("fitness_job",jobs,active,compact=FALSE)
  included_samples <- reactive({d<-active()$samples;inc<-if("include"%in%names(d))tolower(as.character(d$include))%in%c("true","1","yes")else rep(TRUE,nrow(d));as.character(d$sample[inc])})
  output$active_project_badge <- renderUI({p<-active();tagList(tags$b(p$project_id),tags$span(p$species),tags$span(if(p$analysis_ready)"analysis ready" else "preprocessing"))})
  output$overview_ui <- renderUI({p<-active();s<-sample_pipeline_status();ref<-reference_readiness(p,repo_root);mc<-stage_counts(s,"mapping");hc<-stage_counts(s,"hit_file");sc<-stage_counts(s,"summary");progress<-if(sc$total)round(100*(mc$complete+hc$complete+sc$complete)/(3*sc$total))else 0;tagList(tags$div(class="ytab-stat-grid",tags$div(tags$b(p$project_id)," Project"),tags$div(tags$b(p$included_count)," Included samples"),tags$div(tags$b(ref$label)," Reference"),tags$div(tags$b(if(sc$complete==sc$total&&sc$total>0)"Ready"else"Not ready")," Analysis")),tags$div(class="progress",tags$div(class="progress-bar",style=sprintf("width:%d%%",progress),paste0(progress,"%"))),tags$p(tags$b("Next action: "),next_preprocessing_stage(p,repo_root,s)))})
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
  observeEvent(input$enter_analysis_browser,go_to("quality_control","qc_tabs","mapping_qc"),ignoreInit=TRUE);observeEvent(input$enter_analysis_browser_summary,go_to("quality_control","qc_tabs","mapping_qc"),ignoreInit=TRUE)
  observeEvent(input$overview_continue,{destination<-resolve_continue_destination(active(),repo_root,sample_pipeline_status());go_to(destination$top,"preprocessing_tabs",destination$nested)},ignoreInit=TRUE)
  output$qc_selected_samples<-renderUI(selected_names_ui())
  output$mapping_qc_readiness<-renderUI({if(!nrow(mapping_qc_data(active())))tags$p(class="text-muted","Mapping QC becomes available after MapFastq completes.")})
  output$summary_qc_readiness<-renderUI({if(!nrow(summary_qc_data(active())))tags$p(class="text-muted","Summary QC becomes available after raw SummaryTable completes.")})
  mapping_data<-reactive(mapping_qc_data(active()));summary_data<-reactive(summary_qc_data(active()))
  output$mapping_qc_selected_plot<-renderUI(ytab_plot_frame(plotOutput("mapping_qc_stats_plot",width="100%",height=ytab_plot_height_px(input$mapping_qc_plot_height%||%"medium")),input$mapping_qc_plot_width%||%"standard","app-rendered"))
  output$mapping_qc_stats_plot<-renderPlot(ytab_with_plot_display_options(input,"mapping_qc",plot_qc_mapping_stats(active())))
  output$mapping_qc_table<-DT::renderDT({data<-mapping_data();if(!nrow(data))return(NULL);compact_qc_table(data)})
  output$mapping_qc_details<-DT::renderDT({details<-attr(mapping_data(),"details");if(is.null(details)||!nrow(details))return(NULL);compact_qc_table(details)})
  output$download_mapping_qc_table<-downloadHandler(filename=function()"mapping_qc_summary.csv",content=function(file)write.csv(mapping_data(),file,row.names=FALSE))
  output$download_mapping_qc_details<-downloadHandler(filename=function()"mapping_qc_file_details.csv",content=function(file){details<-attr(mapping_data(),"details");if(is.null(details))details<-data.frame();write.csv(details,file,row.names=FALSE)})
  output$summary_qc_cards<-renderUI(summary_qc_cards(active()))
  output$summary_qc_selected_plot<-renderUI({
    choice<-input$summary_qc_plot_choice%||%"complexity"
    output_id<-switch(choice,features="summary_qc_features_plot",combined_features="summary_qc_combined_features_plot",reads_per_hit="summary_qc_reads_per_hit_plot",feature_intergenic="summary_qc_feature_intergenic_plot",genome_bins="summary_qc_genome_bins_plot",pairwise="summary_qc_pairwise_plot","summary_qc_complexity_plot")
    default_height<-if(choice%in%c("pairwise","genome_bins"))"large"else"medium"
    ytab_plot_frame(plotOutput(output_id,width="100%",height=ytab_plot_height_px(input$summary_qc_plot_height%||%default_height,default=default_height)),input$summary_qc_plot_width%||%"standard","app-rendered")
  })
  output$summary_qc_complexity_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_metric(active(),"complexity")))
  output$summary_qc_features_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_metric(active(),"features")))
  output$summary_qc_combined_features_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_combined_features_hit(active())))
  output$summary_qc_reads_per_hit_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_metric(active(),"reads_per_hit")))
  output$summary_qc_feature_intergenic_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_metric(active(),"feature_intergenic")))
  output$summary_qc_genome_bins_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_genome_bins(active())))
  output$summary_qc_pairwise_plot<-renderPlot(ytab_with_plot_display_options(input,"summary_qc",plot_qc_summary_pairwise_correlations(active())))
  output$summary_qc_table<-DT::renderDT({data<-summary_data();if(!nrow(data))return(NULL);compact_qc_table(data)})
  output$summary_qc_details<-DT::renderDT({details<-attr(summary_data(),"details");if(is.null(details)||!nrow(details))return(NULL);compact_qc_table(details)})
  output$download_summary_qc_table<-downloadHandler(filename=function()"summary_qc_table.csv",content=function(file)write.csv(summary_data(),file,row.names=FALSE))
  output$download_summary_qc_details<-downloadHandler(filename=function()"summary_qc_detailed_metrics.csv",content=function(file){details<-attr(summary_data(),"details");if(is.null(details))details<-data.frame();write.csv(details,file,row.names=FALSE)})
  output$library_diagnostics_selected_plot<-renderUI({
    choice<-input$library_diagnostics_plot_choice%||%"midlc"
    diagnostic_height<-ytab_plot_height_px(input$library_diagnostics_plot_height%||%"medium")
    if(identical(choice,"jackpot"))return(ytab_plot_frame(plotOutput("library_jackpot_depth_plot",width="100%",height=diagnostic_height),input$library_diagnostics_plot_width%||%"standard","app-rendered"))
    if(choice%in%c("centromere","metaplots")){
      plot_type<-if(identical(choice,"centromere"))"Centromere bias"else"Feature metaplots"
      plots<-qc_library_plot_inventory(active(),plot_type)
      if(!nrow(plots))return(tags$p(class="text-muted","No generated diagnostic PNGs are available for this plot type."))
      labels<-paste(plots$sample,plots$filename,sep=" — ")
      return(tagList(selectInput("library_diagnostic_png_choice","Sample plot",choices=as.list(setNames(seq_len(nrow(plots)),labels))),uiOutput("library_diagnostic_png_selected")))
    }
    if(identical(choice,"sequence_bias"))return(ytab_plot_frame(plotOutput("library_sequence_bias_plot",width="100%",height=diagnostic_height),input$library_diagnostics_plot_width%||%"standard","app-rendered"))
    ytab_plot_frame(plotOutput("library_midlc_plot",width="100%",height=diagnostic_height),input$library_diagnostics_plot_width%||%"standard","app-rendered")
  })
  output$library_diagnostic_png_selected<-renderUI({
    choice<-input$library_diagnostics_plot_choice%||%""
    plot_type<-if(identical(choice,"centromere"))"Centromere bias"else if(identical(choice,"metaplots"))"Feature metaplots"else""
    if(!nzchar(plot_type))return(NULL)
    plots<-qc_library_plot_inventory(active(),plot_type)
    index<-suppressWarnings(as.integer(input$library_diagnostic_png_choice%||%1L))
    qc_library_single_plot_card(plots,index)
  })
  output$library_midlc_plot<-renderPlot(ytab_with_plot_display_options(input,"library_diagnostics",plot_qc_library_midlc(active())))
  output$library_jackpot_depth_plot<-renderPlot(ytab_with_plot_display_options(input,"library_diagnostics",plot_qc_library_jackpot_depth(active())))
  output$library_sequence_bias_plot<-renderPlot(ytab_with_plot_display_options(input,"library_diagnostics",plot_qc_library_sequence_bias(active())))
  output$library_diagnostic_plot_gallery<-renderUI(qc_library_plot_gallery(active()))
  diagnostic_inventory<-reactive(diagnostic_inventory_state())
  observeEvent(input$refresh_diagnostic_files,{diagnostic_inventory_state(build_diagnostic_file_inventory(active()));diagnostic_gallery_page(1L)},ignoreInit=TRUE)
  output$diagnostic_result_selector<-renderUI({x<-diagnostic_inventory();runs<-sort(unique(x$run_id));if(length(runs)>1){labels<-ifelse(runs=="current_legacy","Legacy diagnostics result",tools::toTitleCase(gsub("_"," ",runs)));selectInput("diagnostic_run_filter","Diagnostics result",choices=c(list("All results"="All"),as.list(setNames(unname(runs),unname(labels)))),selected=diagnostic_selected_run())}})
  observeEvent(input$diagnostic_run_filter,{diagnostic_selected_run(input$diagnostic_run_filter%||%"All");diagnostic_gallery_page(1L)},ignoreInit=TRUE)
  output$diagnostic_file_filters<-renderUI({x<-diagnostic_inventory();if(!nrow(x))return(NULL);controls<-list();samples<-sort(unique(x$sample));types<-sort(unique(x$display_type));sets<-sort(unique(x$sample_set));if(length(samples)>1)controls<-c(controls,list(selectInput("diagnostic_sample_filter","Sample",c("All",samples))));if(identical(input$diagnostic_view_mode,"table")&&length(types)>1)controls<-c(controls,list(selectInput("diagnostic_type_filter","Type",c("All",types))));if(length(sets)>1)controls<-c(controls,list(selectInput("diagnostic_set_filter","Sample set",c("All",sets))));plot_types<-sort(unique(x$plot_type[x$is_plot]));if(identical(input$diagnostic_view_mode,"gallery")&&length(plot_types)>1)controls<-c(controls,list(selectInput("diagnostic_plot_type_filter","Plot type",c("All",plot_types))));tags$div(class="ytab-filter-row",controls)})
  filtered_diagnostic_inventory<-reactive({gallery<-identical(input$diagnostic_view_mode,"gallery");filter_diagnostic_inventory(diagnostic_inventory(),sample=input$diagnostic_sample_filter%||%"All",plot_type=if(gallery)input$diagnostic_plot_type_filter%||%"All"else"All",sample_set=input$diagnostic_set_filter%||%"All",run_id=input$diagnostic_run_filter%||%diagnostic_selected_run(),type=if(gallery)"All"else input$diagnostic_type_filter%||%"All")})
  gallery_page_data<-reactive({x<-filtered_diagnostic_inventory();x<-x[x$is_plot,,drop=FALSE];paginate_diagnostic_inventory(x,diagnostic_gallery_page(),8L)})
  output$diagnostic_files_empty<-renderUI({x<-filtered_diagnostic_inventory();if(!nrow(diagnostic_inventory()))return(tags$p(class="text-muted","No diagnostic files are available for the selected run."));if(!nrow(x))return(tags$p(class="text-muted","No diagnostic files match the current filters."));if(identical(input$diagnostic_view_mode,"gallery")&&!any(x$is_plot))tags$p(class="text-muted","No plot files are available. Switch to File table to view diagnostic tables.")})
  output$diagnostic_files_table<-DT::renderDT({x<-filtered_diagnostic_inventory();if(!nrow(x))return(NULL);visible<-data.frame(File=x$filename,Type=x$display_type,`Sample set`=x$sample_set,`Sample or project`=x$sample,Size=x$size_display,Modified=x$modified,Action=ifelse(x$viewable,"View · Show path","Show path"),check.names=FALSE);compact_qc_table(visible)})
  output$diagnostic_plot_cards<-renderUI({
    plots<-gallery_page_data()$data
    if(!nrow(plots))return(NULL)
    tags$div(class="ytab-plot-grid ytab-diagnostic-gallery-grid",lapply(seq_len(nrow(plots)),function(i){
      preview<-if(isTRUE(plots$preview_available[[i]]))tagList(
        tags$img(src=plots$served_url[[i]],alt=plots$filename[[i]],loading="lazy",style="width:100%;max-height:300px;object-fit:contain",onerror="this.style.display='none';this.nextElementSibling.style.display='block';"),
        tags$p(class="ytab-preview-unavailable",style="display:none",plots$preview_message[[i]],tags$br(),tags$code(plots$relative_path[[i]]))
      )else tags$p(class="ytab-preview-unavailable",plots$preview_message[[i]],tags$br(),tags$code(plots$relative_path[[i]]))
      tags$article(class="ytab-static-image-card ytab-release-card",title=plots$filename[[i]],tags$div(class="ytab-plot-card-header",tags$div(tags$h4(plots$plot_type[[i]]),tags$span(class="ytab-status-badge","Static generated image")),if(isTRUE(plots$preview_available[[i]]))tags$a(class="btn btn-secondary btn-sm",href=plots$served_url[[i]],target="_blank","Open original")),tags$p(class="ytab-plot-sample",tags$b("Sample: "),plots$sample[[i]]),tags$p(class="text-muted",tags$b("Sample set: "),plots$sample_set[[i]]),
        tags$div(class="ytab-diagnostic-preview",preview),
        actionButton(paste0("view_diagnostic_plot_",i),"View larger",class="btn-secondary",disabled=if(isTRUE(plots$preview_available[[i]]))NULL else"disabled"),
        tags$details(tags$summary("Show path / filename"),tags$code(plots$relative_path[[i]]),tags$br(),tags$small(plots$filename[[i]])))
    }))
  })
  output$diagnostic_gallery_pagination<-renderUI({page<-gallery_page_data();if(page$total<=8L)return(NULL);tags$div(class="ytab-gallery-pagination",actionButton("diagnostic_gallery_previous","Previous",disabled=if(page$page<=1L)"disabled"else NULL),tags$span(sprintf("Page %d of %d",page$page,page$pages)),actionButton("diagnostic_gallery_next","Next",disabled=if(page$page>=page$pages)"disabled"else NULL))})
  observeEvent(input$diagnostic_gallery_previous,diagnostic_gallery_page(max(1L,diagnostic_gallery_page()-1L)),ignoreInit=TRUE);observeEvent(input$diagnostic_gallery_next,diagnostic_gallery_page(min(gallery_page_data()$pages,diagnostic_gallery_page()+1L)),ignoreInit=TRUE)
  for(i in seq_len(8L))local({index<-i;observeEvent(input[[paste0("view_diagnostic_plot_",index)]],{plots<-gallery_page_data()$data;if(index>nrow(plots))return();row<-plots[index,,drop=FALSE];if(!isTRUE(row$preview_available[[1]]))return();showModal(modalDialog(title=row$display_title[[1]],tags$img(src=row$served_url[[1]],alt=row$filename[[1]],style="width:100%;max-height:70vh;object-fit:contain",onerror="this.style.display='none';this.nextElementSibling.style.display='block';"),tags$p(class="ytab-preview-unavailable",style="display:none",row$preview_message[[1]],tags$br(),tags$code(row$relative_path[[1]])),tags$dl(class="ytab-meta",tags$dt("Sample"),tags$dd(row$sample[[1]]),tags$dt("Plot type"),tags$dd(row$plot_type[[1]]),tags$dt("Result"),tags$dd(row$sample_set[[1]]),tags$dt("Filename"),tags$dd(row$filename[[1]]),tags$dt("Path"),tags$dd(tags$code(row$relative_path[[1]]))),easyClose=TRUE,footer=modalButton("Close")))},ignoreInit=TRUE)})

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
  fitness_data<-reactive(fitness_result_data(current_fitness_result()));fitness_filtered<-reactive(filter_fitness_results(fitness_data(),input$fitness_call_filter%||%"All",input$fitness_gene_search%||%""))
  observe({data<-fitness_data();column<-fitness_call_column(data);calls<-if(nzchar(column))sort(unique(as.character(data[[column]])))else character();calls<-calls[!is.na(calls)&nzchar(calls)];choices<-as.list(c("All",unname(calls)));names(choices)<-c("All calls",tools::toTitleCase(gsub("_"," ",calls)));updateSelectInput(session,"fitness_call_filter",choices=choices)})
  output$fitness_summary_cards<-renderUI({data<-fitness_data();if(!nrow(data))return(NULL);counts<-fitness_call_counts(data);card<-function(key,label)tags$div(tags$b(counts[[key]]),label);tags$div(class="ytab-stat-grid",tags$div(tags$b(nrow(data)),"Total features"),card("consistently_depleted","Consistently depleted"),card("consistently_enriched","Consistently enriched"),card("single_pool_depleted","Single-pool depleted"),card("single_pool_enriched","Single-pool enriched"),card("mixed","Mixed"))})
  output$fitness_filtered_count<-renderUI(if(nrow(fitness_data()))tags$p(sprintf("Showing %d of %d features.",nrow(fitness_filtered()),nrow(fitness_data()))))
  output$fitness_results_table<-DT::renderDT({x<-fitness_result_table_data(fitness_filtered());if(!nrow(x))return(NULL);compact_qc_table(x)})
  output$fitness_visualization_selector<-renderUI({
    result<-current_fitness_result()
    ma_path<-fitness_combined_ma_plot_file(result)
    plots<-fitness_generated_plot_inventory_for_result(active(),result)
    choices<-c(if(nzchar(ma_path))c("Combined treated-versus-parent MA plot"="combined_ma"),"Condition versus control log-log scatter"="condition_control","Fitness call distribution"="call_distribution","Effect size versus z-score"="effect_size")
    if(nrow(plots))choices<-c(choices,setNames(paste0("generated:",seq_len(nrow(plots))),paste("Generated:",plots$title)))
    selectInput("fitness_visualization_choice","Plot",choices=choices,selected=if(nzchar(ma_path))"combined_ma"else"call_distribution")
  })
  output$fitness_selected_visualization<-renderUI({
    choice<-input$fitness_visualization_choice%||%"call_distribution"
    fitness_height<-ytab_plot_height_px(input$fitness_plot_height%||%"large",default="large")
    if(identical(choice,"combined_ma"))return(fitness_combined_ma_plot_card(current_fitness_result(),active()))
    if(startsWith(choice,"generated:")){
      plots<-fitness_generated_plot_inventory_for_result(active(),current_fitness_result());index<-suppressWarnings(as.integer(sub("^generated:","",choice)))
      if(is.na(index)||index<1L||index>nrow(plots))return(tags$p(class="text-muted","Selected plot is unavailable."))
      return(tags$article(class="ytab-static-image-card ytab-release-card",tags$div(class="ytab-plot-card-header",tags$h4(paste("Generated:",plots$title[[index]])),tags$span(class="ytab-status-badge","Static generated image")),tags$div(class="ytab-diagnostic-preview",tags$img(src=plots$served_url[[index]],alt=plots$filename[[index]],loading="lazy",style=paste0("width:100%;max-height:",fitness_height,";object-fit:contain"))),tags$details(tags$summary("Show filename"),tags$code(plots$filename[[index]]))))
    }
    if(identical(choice,"condition_control"))return(fitness_condition_control_plot_card(current_fitness_result(),active()))
    if(identical(choice,"effect_size"))return(tagList(uiOutput("fitness_effect_plot_state"),ytab_plot_frame(plotOutput("fitness_effect_plot",width="100%",height=fitness_height),input$fitness_plot_width%||%"standard","app-rendered")))
    ytab_plot_frame(plotOutput("fitness_call_plot",width="100%",height=fitness_height),input$fitness_plot_width%||%"standard","app-rendered")
  })
  output$fitness_call_plot<-renderPlot(ytab_with_plot_display_options(input,"fitness",{data<-fitness_data();req(nrow(data));counts<-fitness_call_counts(data);counts<-counts[counts>0];if(!length(counts)){qc_plot_empty("No fitness-call values are available.");return()};horizontal<-qc_plot_bar_horizontal(names(counts));old<-qc_plot_par(mar=if(horizontal)c(5,12,3,1)else c(9,5,3,1));on.exit(par(old),add=TRUE);barplot(counts,las=if(horizontal)1 else qc_plot_label_las(2),horiz=horizontal,col=qc_plot_fill,border=qc_plot_border,lwd=qc_plot_lwd,main=sprintf("Fitness calls (%d features)",nrow(data)),xlab=if(horizontal)"Features"else"",ylab=if(horizontal)""else"Features")}))
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
