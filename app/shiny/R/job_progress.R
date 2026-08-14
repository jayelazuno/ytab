read_progress_file <- function(path) read_json_safely(path)

essentiality_job_matches <- function(job, expected_stage) {
  expected_stage <- as.character(expected_stage %||% "")
  if (length(expected_stage) != 1L || !nzchar(expected_stage)) return(TRUE)
  if (is.null(job) || !is.list(job)) return(FALSE)
  stage <- as.character(job$stage %||% "")
  length(stage) == 1L && identical(stage, expected_stage)
}

essentiality_stage_job_history <- function(project_root, expected_stage) {
  empty <- data.frame(
    Time = character(), Stage = character(), Target = character(),
    Mode = character(), Status = character(), Elapsed = character(),
    Action = character(), stringsAsFactors = FALSE, check.names = FALSE
  )
  expected_stage <- as.character(expected_stage %||% "")
  allowed <- c("sample_normalization", "summary_normalized", "combined_hits",
               "summary_combined", "classifier")
  if (length(expected_stage) != 1L || !expected_stage %in% allowed) return(empty)
  root <- file.path(project_root, "manifests", "jobs")
  if (!dir.exists(root)) return(empty)
  paths <- list.files(root, pattern = "\\.progress\\.json$", full.names = TRUE)
  rows <- lapply(paths, function(path) {
    progress <- read_progress_file(path)
    if (!essentiality_job_matches(progress, expected_stage)) return(NULL)
    scalar <- function(value, default = "") {
      value <- unlist(value %||% default, use.names = FALSE)
      if (length(value)) as.character(value[[1L]]) else default
    }
    elapsed <- suppressWarnings(as.numeric(progress$job_elapsed_seconds %||% NA_real_))
    data.frame(
      Time = scalar(progress$updated_at %||% progress$job_started_at, "Unknown"),
      Stage = progress_stage_label(expected_stage),
      Target = scalar(progress$target_tag %||% progress$target, "Not recorded"),
      Mode = scalar(progress$execution_mode,
                    if (isTRUE(progress$dry_run)) "preview" else "run"),
      Status = scalar(progress$status, "unknown"),
      Elapsed = if (is.finite(elapsed)) format_duration(elapsed) else "Unavailable",
      Action = "Review",
      stringsAsFactors = FALSE, check.names = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(empty)
  result <- do.call(rbind, rows)
  parsed <- suppressWarnings(as.POSIXct(result$Time, tz = "UTC"))
  order_key <- ifelse(is.na(parsed), 0, as.numeric(parsed))
  result[order(order_key, decreasing = TRUE), , drop = FALSE]
}

format_duration <- function(seconds) {
  if(is.null(seconds)||!length(seconds)||is.na(seconds))return("unavailable")
  seconds<-max(0,as.numeric(seconds));if(seconds<300)sprintf("%dm %02ds",floor(seconds/60),round(seconds%%60/10)*10)else if(seconds<3600)sprintf("%d min",round(seconds/60))else sprintf("%.1f hr",round(seconds/1800)/2)
}
format_eta <- function(progress) {
  if(is.null(progress)||is.null(progress$eta_seconds))return("Estimating after the first sample completes.")
  paste("Approx.",format_duration(progress$eta_seconds),"remaining")
}
progress_is_stale <- function(progress,threshold=60) {
  if(is.null(progress)||is.null(progress$updated_at))return(FALSE)
  updated<-tryCatch(as.POSIXct(progress$updated_at,tz="UTC"),error=function(e)NA);!is.na(updated)&&as.numeric(difftime(Sys.time(),updated,units="secs"))>threshold
}
progress_status_label <- function(status) if(identical(status,"dry_run_success"))"Dry run complete"else tools::toTitleCase(gsub("_"," ",status%||%"unknown"))
progress_status_badge <- function(status) tags$span(class=paste("ytab-status-badge",paste0("status-",status%||%"unknown")),progress_status_label(status))
progress_stage_label <- function(stage) {
  labels <- c(sample_normalization="Parent MidLC normalization",
              summary_normalized="Normalized SummaryTable",
              combined_hits="Combine Parent Libraries",
              summary_combined="Combined SummaryTable",
              classifier="Essentiality Classifier")
  if (length(stage) == 1L && stage %in% names(labels)) unname(labels[[stage]]) else stage
}
progress_display_metrics <- function(progress) {dry<-isTRUE(progress$dry_run)||identical(progress$status,"dry_run_success");list(validated=if(dry)progress$processed_items%||%0 else 0,scientific_completed=if(dry)0 else (progress$successful_items%||%0)+(progress$skipped_items%||%0),executed=if(dry)0 else progress$successful_items%||%0,progress_label=if(dry)"Dry-run validation progress"else"Confirmed processing progress",cancel_active=progress$status%in%c("queued","starting","running"),dismiss_available=!progress$status%in%c("queued","starting","running"))}
tail_job_log <- function(path,max_lines=100) {if(!nzchar(path%||%"")||!file.exists(path))return(character());lines<-readLines(path,warn=FALSE);tail(lines,max_lines)}
as_time_or_na <- function(value) {if(is.null(value)||!length(value)||!nzchar(value))return(as.POSIXct(NA));tryCatch(as.POSIXct(value,tz="UTC"),error=function(e)as.POSIXct(NA))}

job_progress_ui <- function(id,compact=FALSE) {ns<-NS(id);tags$div(class=if(compact)"ytab-job-progress compact"else"ytab-job-progress",uiOutput(ns("body")))}
job_progress_server <- function(id,job_manager,project,compact=FALSE,expected_stage=NULL) moduleServer(id,function(input,output,session){
  progress<-reactive({running<-job_manager$job_is_running();if(running)invalidateLater(if(compact)1500 else 1000,session);job_manager$poll_job();p<-job_manager$current_progress();if(!essentiality_job_matches(p,expected_stage))return(NULL);p})
  output$body<-renderUI({p<-progress();job<-job_manager$current_job();if(is.null(p)){if(compact)return(NULL);return(tags$p(class="text-muted","No job progress is available."))}
    if(compact&&!p$status%in%c("queued","starting","running"))return(NULL)
    now<-Sys.time();started<-as_time_or_na(p$job_started_at);item_started<-as_time_or_na(p$current_item_started_at);live_job_elapsed<-if(!is.na(started)&&p$status%in%c("queued","starting","running"))as.numeric(difftime(now,started,units="secs"))else p$job_elapsed_seconds;live_item_elapsed<-if(!is.na(item_started)&&p$status=="running")as.numeric(difftime(now,item_started,units="secs"))else p$current_item_elapsed_seconds
    if(!is.null(p$eta_seconds)&&p$status=="running"&&!is.null(p$updated_at)){updated<-as_time_or_na(p$updated_at);if(!is.na(updated))p$eta_seconds<-max(0,p$eta_seconds-as.numeric(difftime(now,updated,units="secs")))}
    current<-p$current_item%||%"No active item";index<-p$current_item_index%||%0;total<-p$total_items%||%0;processed<-p$processed_items%||%0;pct<-p$progress_percent%||%0
    metrics<-progress_display_metrics(p);stale<-progress_is_stale(p)&&job_manager$job_is_running();eta<-format_eta(p);finish_time<-as_time_or_na(p$estimated_completion_time);finish<-if(is.na(finish_time))"unavailable"else format(finish_time,"%I:%M %p")
    if(compact)return(tags$div(class="ytab-job-banner-content",progress_status_badge(p$status),tags$b(progress_stage_label(p$stage)),tags$span(current),tags$span(sprintf("Item %s of %s",index,total)),tags$span(sprintf("%d of %d processed",processed,total)),tags$span(paste("Elapsed",format_duration(live_job_elapsed))),tags$span(eta),actionButton(session$ns("cancel"),"Cancel")))
    heading<-tags$div(class="ytab-job-heading",tags$h4(paste(if(metrics$cancel_active)"Active"else"Last",progress_stage_label(p$stage),"job")),progress_status_badge(p$status))
    if((p$stage%||%"")%in%c("sample_normalization","summary_normalized","combined_hits","summary_combined","classifier")){
      item_label<-switch(p$stage,sample_normalization="Normalization sweep",summary_normalized="Target evaluation",combined_hits="Target",summary_combined="Target",classifier="Target")
      return(tagList(heading,tags$div(class="progress",tags$div(class=paste("progress-bar",if(p$status=="running")"progress-bar-striped progress-bar-animated"),style=sprintf("width:%s%%",pct),paste0(pct,"%"))),tags$div(class="ytab-stat-grid",tags$div(tags$b(p$execution_mode%||%if(isTRUE(p$dry_run))"preview"else"run"),"Execution mode"),tags$div(tags$b(current),item_label),tags$div(tags$b(p$parent_count%||%length(p$parent_samples%||%character())),"Parent libraries"),tags$div(tags$b(p$processed_items%||%0),"Completed"),tags$div(tags$b(p$skipped_items%||%0),"Skipped"),tags$div(tags$b(p$failed_items%||%0),"Failed"),tags$div(tags$b(p$remaining_items%||%0),"Remaining"),tags$div(tags$b(format_duration(live_job_elapsed)),"Elapsed")),tags$p(tags$b("Current phase: "),p$current_phase%||%"unknown"),tags$p(p$message%||%""),if(stale)tags$div(class="alert alert-warning","No progress update has been received for more than one minute. The process is still running."),if(metrics$cancel_active)actionButton(session$ns("cancel"),"Cancel job")else actionButton(session$ns("dismiss"),"Dismiss job summary"),tags$details(open=if(identical(p$status,"failed"))NA else NULL,tags$summary("Recent stdout/stderr"),tags$pre(paste(c(job_manager$current_stdout(),job_manager$current_stderr()),collapse="\n"))),tags$p("Progress file: ",tags$code(job$progress_file%||%"unavailable"))))
    }
    if(identical(tolower(p$stage%||%""),"library diagnostics"))return(tagList(heading,tags$div(class="ytab-stat-grid",tags$div(tags$b(p$execution_mode%||%if(isTRUE(p$dry_run))"preview"else"run"),"Execution mode"),tags$div(tags$b(p$sample_set_label%||%current),"Sample set"),tags$div(tags$b(p$selected_sample_count%||%length(p$selected_items%||%character())),"Selected samples"),tags$div(tags$b(format_duration(live_job_elapsed)),"Elapsed")),tags$p(tags$b("Current phase: "),p$current_phase%||%"unknown"),tags$p(p$message%||%""),if(stale)tags$div(class="alert alert-warning","No progress update has been received for more than one minute. The process is still running."),if(metrics$cancel_active)actionButton(session$ns("cancel"),"Cancel job")else actionButton(session$ns("dismiss"),"Dismiss job summary"),tags$details(tags$summary("Recent stdout/stderr"),tags$pre(paste(c(job_manager$current_stdout(),job_manager$current_stderr()),collapse="\n")))))
    if(identical(tolower(p$stage%||%""),"fitness analysis")){message<-fitness_job_status_message(p$status%||%"stale");failed<-identical(p$status,"failed");error<-p$items[[1]]$error_message%||%p$error_message%||%"Review recent stdout/stderr for the failure.";return(tagList(heading,tags$div(class="ytab-stat-grid",tags$div(tags$b(p$analysis_id%||%current),"Analysis ID"),tags$div(tags$b(p$selected_comparison_count%||%p$comparison_count%||%length(p$selected_comparisons%||%character())),"Selected comparisons"),tags$div(tags$b(p$execution_mode%||%if(isTRUE(p$dry_run))"preview"else"run"),"Execution mode"),tags$div(tags$b(if(isTRUE(p$classifier_annotation))"Yes"else"No"),"Classifier annotation"),tags$div(tags$b(format_duration(live_job_elapsed)),"Elapsed"),if(failed)tags$div(tags$b(job_manager$job_exit_status()%||%"nonzero"),"Exit status")),tags$p(tags$b("Current phase: "),p$current_phase%||%if(failed)"runner_failed"else"unknown"),tags$p(class=if(failed)"text-danger"else NULL,tags$b(message)),if(failed)tags$p(error),if(stale)tags$div(class="alert alert-warning","No progress update has been received for more than one minute. The process is still running."),if(metrics$cancel_active)actionButton(session$ns("cancel"),"Cancel job")else actionButton(session$ns("dismiss"),"Dismiss job summary"),tags$details(open=if(failed)NA else NULL,tags$summary("Recent stdout/stderr"),tags$pre(paste(c(job_manager$current_stdout(),job_manager$current_stderr()),collapse="\n"))))) }
    tagList(heading,tags$p(tags$b(metrics$progress_label)),tags$div(class="progress",tags$div(class=paste("progress-bar",if(p$status=="running")"progress-bar-striped progress-bar-animated"),style=sprintf("width:%s%%",pct),paste0(pct,"%"))),tags$div(class="ytab-stat-grid",tags$div(tags$b(current),"Current sample"),tags$div(tags$b(sprintf("%s / %s",index,total)),"Item"),tags$div(tags$b(metrics$validated),"Validated"),tags$div(tags$b(metrics$scientific_completed),"Scientifically completed"),tags$div(tags$b(metrics$executed),"Executed"),tags$div(tags$b(p$failed_items%||%0),"Failed"),tags$div(tags$b(p$remaining_items%||%0),"Remaining"),tags$div(tags$b(format_duration(live_job_elapsed)),"Elapsed")),tags$p(tags$b("Current phase: "),p$current_phase%||%"unknown"),tags$p(tags$b("Current-item elapsed: "),format_duration(live_item_elapsed)),tags$p(tags$b("Approx. remaining: "),eta),tags$p(tags$b("Expected finish: "),finish),tags$p(tags$b("ETA confidence: "),p$eta_confidence%||%"unavailable"),tags$p(tags$b("Last update: "),p$updated_at%||%"unavailable"),if(stale)tags$div(class="alert alert-warning","No progress update has been received for more than one minute. The process is still running."),tags$p(p$message%||%""),tags$p(class="text-muted","Completion estimates are based on durations and input sizes completed during the current run. Actual times can vary with storage and system load."),if(metrics$cancel_active)actionButton(session$ns("cancel"),"Cancel job")else actionButton(session$ns("dismiss"),"Dismiss job summary"),tags$details(tags$summary("Recent stdout"),tags$pre(paste(job_manager$current_stdout(),collapse="\n"))),tags$details(tags$summary("Recent stderr"),tags$pre(paste(job_manager$current_stderr(),collapse="\n"))),tags$p("Progress file: ",tags$code(job$progress_file%||%"unavailable")))})
  observeEvent(input$cancel,job_manager$cancel_job(),ignoreInit=TRUE);observeEvent(input$dismiss,job_manager$dismiss_job(),ignoreInit=TRUE);progress
})
