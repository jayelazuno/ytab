json_safe_value <- function(value) {
  if(is.data.frame(value))return(lapply(value,function(x)unname(x)))
  if(is.list(value))return(lapply(value,json_safe_value))
  if(is.atomic(value)&&length(names(value))&&any(nzchar(names(value))))return(lapply(as.list(value),json_safe_value))
  if(is.atomic(value))return(unname(value))
  value
}
atomic_json_write <- function(path, value) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile(pattern=paste0(".",basename(path),"."),tmpdir=dirname(path))
  writeLines(jsonlite::toJSON(json_safe_value(value),auto_unbox=TRUE,pretty=TRUE,null="null"),tmp,useBytes=TRUE)
  if(!file.rename(tmp,path)){unlink(tmp);stop("Could not atomically write job state: ",path)}
}

read_json_safely <- function(path) {
  if(!nzchar(path %||% "")||!file.exists(path))return(NULL)
  tryCatch(jsonlite::fromJSON(path,simplifyVector=FALSE),error=function(e)NULL)
}

job_terminal_statuses <- function() c(
  "success", "cached", "cached_success", "dry_run_success", "preview_success",
  "partial", "failed", "cancelled", "stale"
)
job_active_statuses <- function() c("queued", "starting", "running")
job_status_is_terminal <- function(status) as.character(status %||% "") %in% job_terminal_statuses()
job_status_is_active <- function(status) as.character(status %||% "") %in% job_active_statuses()
job_process_exists <- function(pid) {
  pid <- suppressWarnings(as.integer(pid %||% NA_integer_))
  if (length(pid) != 1L || is.na(pid) || pid <= 0L ||
      !requireNamespace("processx", quietly = TRUE)) return(FALSE)
  exists_function <- if ("process_exists" %in% getNamespaceExports("processx"))
    getExportedValue("processx", "process_exists") else
    getFromNamespace("process__exists", "processx")
  isTRUE(exists_function(pid))
}

job_heartbeat_is_stale <- function(progress, now = Sys.time(), threshold = 120) {
  updated <- as.character(progress$updated_at %||% "")
  if (!nzchar(updated)) return(TRUE)
  stamp <- suppressWarnings(as.POSIXct(updated, tz = "UTC"))
  is.na(stamp) || as.numeric(difftime(now, stamp, units = "secs")) > threshold
}

reconcile_job_status <- function(current_status, progress_status = "",
                                 process_alive = FALSE, valid_output = FALSE,
                                 cache_reused = FALSE, heartbeat_stale = FALSE) {
  if (job_status_is_terminal(progress_status)) return(as.character(progress_status))
  if (isTRUE(process_alive)) return("running")
  if (isTRUE(valid_output)) return(if (isTRUE(cache_reused)) "cached" else "success")
  if (job_status_is_active(current_status) && (isTRUE(heartbeat_stale) || !isTRUE(process_alive)))
    return("stale")
  if (job_status_is_terminal(current_status)) return(as.character(current_status))
  as.character(current_status %||% "unknown")
}

job_stage_has_valid_output <- function(project_root, stage) {
  patterns <- list(
    mapfastq = c("mapfastq", "\\.sorted\\.bam$"),
    create_hit_file = c("create_hit_file", "_hits\\.txt$"),
    summary = c("summary", "feature_table.*\\.(csv|tsv|txt)$"),
    library_diagnostics = c("library_diagnostics", "\\.(png|csv|tsv)$"),
    sample_normalization = c("sample_normalization", "_normalized_hits\\.txt$"),
    summary_normalized = c("sample_normalization", "normalization_target_evaluation\\.csv$"),
    combined_hits = c("combined_hits", "combined_parent_hits.*\\.txt$"),
    summary_combined = c("summary_combined", "combined_feature_table.*\\.txt$"),
    classifier = c("classifier", "essentiality_predictions.*\\.csv$"),
    fitness_analysis = c("treated_vs_parent", "treated_vs_parent_results\\.csv$")
  )
  spec <- patterns[[as.character(stage %||% "")]]
  if (is.null(spec)) return(FALSE)
  root <- file.path(project_root, spec[[1L]])
  outputs <- if (dir.exists(root)) list.files(root, pattern = spec[[2L]], recursive = TRUE,
                                              full.names = TRUE, ignore.case = TRUE) else character()
  outputs <- outputs[file.exists(outputs) & !dir.exists(outputs) & file.info(outputs)$size > 0]
  if (!length(outputs)) return(FALSE)
  manifest_stage <- switch(as.character(stage),
    mapfastq = "mapfastq", create_hit_file = "create_hit_file", summary = "summary_table",
    library_diagnostics = "library_diagnostics", sample_normalization = "sample_normalization",
    summary_normalized = "summary_normalized", combined_hits = "combined_hits",
    summary_combined = "summary_combined", classifier = "classifier",
    fitness_analysis = "treated_vs_parent", "")
  manifest_root <- file.path(project_root, "manifests", manifest_stage)
  manifests <- if (nzchar(manifest_stage) && dir.exists(manifest_root))
    list.files(manifest_root, pattern = "\\.json$", recursive = TRUE, full.names = TRUE) else character()
  any(vapply(manifests, function(path) {
    value <- read_json_safely(path)
    as.character(value$status %||% "") %in% c("success", "cached", "skipped")
  }, logical(1)))
}

new_job_manager <- function() {
  state <- reactiveValues(process=NULL,command="",args=character(),display_command="",started=NULL,stdout=character(),stderr=character(),project_config="",stage="",selected_items=character(),job_id="",progress_file="",status="idle",pid=NA_integer_,recovered=FALSE,dismissed=FALSE)
  current_job_path <- function() if(nzchar(state$project_config))file.path(dirname(dirname(state$project_config)),"manifests","orchestrator","current_job.json")else""
  raw_progress <- function() read_json_safely(state$progress_file)
  progress <- function() if(isTRUE(state$dismissed))NULL else raw_progress()
  persist_current <- function(message="") {
    if(!nzchar(state$project_config))return(invisible(NULL))
    atomic_json_write(current_job_path(),list(schema_version=2,job_id=state$job_id,stage=state$stage,status=state$status,active=job_status_is_active(state$status),pid=state$pid,project_config=state$project_config,progress_file=state$progress_file,selected_items=state$selected_items,display_command=state$display_command,started_at=if(is.null(state$started))NULL else format(state$started,tz="UTC",usetz=TRUE),updated_at=format(Sys.time(),tz="UTC",usetz=TRUE),cancel_requested=identical(state$status,"cancelled"),message=message))
  }
  list(
    start_job=function(command,args,project_config,stage_or_profile,wd=NULL,selected_items=character(),job_id=NULL,progress_file=NULL,metadata=list()) {
      if (!is.null(state$process) && state$process$is_alive()) stop("A pipeline job is already running.")
      recovered_progress<-raw_progress();if(is.null(state$process)&&!is.null(recovered_progress)&&job_status_is_active(recovered_progress$status)&&job_process_exists(state$pid))stop("A pipeline job is already running.")
      stamp<-format(Sys.time(),"%Y%m%dT%H%M%SZ",tz="UTC");random<-paste(sample(c(letters,0:9),6,replace=TRUE),collapse="");state$job_id<-job_id%||%paste(stage_or_profile,stamp,random,sep="_")
      project_root<-dirname(dirname(project_config));state$progress_file<-progress_file%||%file.path(project_root,"manifests","jobs",paste0(state$job_id,".progress.json"))
      progress_args<-if(stage_or_profile%in%c("mapfastq","create_hit_file","summary","library_diagnostics","fitness_analysis","sample_normalization","summary_normalized","combined_hits","summary_combined","classifier"))c("--job-id",state$job_id,"--progress-file",state$progress_file)else character()
      state$command<-command;state$args<-c(args,progress_args);state$display_command<-format_command_for_display(command,state$args);state$started<-Sys.time();state$stdout<-state$stderr<-character();state$project_config<-project_config;state$stage<-stage_or_profile;state$selected_items<-selected_items;state$status<-"starting";state$recovered<-FALSE;state$dismissed<-FALSE
      dir.create(dirname(state$progress_file),recursive=TRUE,showWarnings=FALSE)
      initial<-c(list(schema_version=1,job_id=state$job_id,project_id=basename(project_root),stage=state$stage,status="queued",selected_items=unname(as.character(selected_items)),total_items=length(selected_items),processed_items=0,successful_items=0,skipped_items=0,failed_items=0,remaining_items=length(selected_items),current_item=NULL,current_item_index=NULL,current_phase="queued",job_started_at=format(Sys.time(),tz="UTC",usetz=TRUE),updated_at=format(Sys.time(),tz="UTC",usetz=TRUE),job_elapsed_seconds=0,current_item_elapsed_seconds=0,progress_fraction=0,progress_percent=0,eta_seconds=NULL,estimated_completion_time=NULL,eta_confidence="unavailable",message="Command submitted",cancel_requested=FALSE,items=lapply(unname(as.character(selected_items)),function(x)list(item=x,status="queued"))),json_safe_value(metadata))
      atomic_json_write(state$progress_file,initial);state$process<-start_process_async(command,state$args,wd=wd);state$pid<-state$process$get_pid();state$status<-"running";persist_current("Process launched");invisible(state$process)
    },
    poll_job=function(){
      p<-raw_progress()
      if(!is.null(p)&&job_status_is_terminal(p$status)&&
         !identical(as.character(state$status),as.character(p$status))){
        state$status<-as.character(p$status);persist_current("Progress reached a terminal state")
      }
      if(!is.null(state$process)){
        state$process$poll_io(0);state$stdout<-tail(c(state$stdout,state$process$read_output_lines()),500);state$stderr<-tail(c(state$stderr,state$process$read_error_lines()),500)
        alive<-state$process$is_alive()
        if(!alive&&job_status_is_active(state$status)){
          exit<-state$process$get_exit_status();state$status<-if(identical(exit,0L))"success"else"failed"
          p<-raw_progress();if(!is.null(p)&&job_status_is_active(p$status)){p$status<-state$status;p$message<-"Process ended before a final progress event.";p$updated_at<-format(Sys.time(),tz="UTC",usetz=TRUE);atomic_json_write(state$progress_file,p)}
          persist_current(paste("Process exited",exit))
        }else if(alive&&job_status_is_active(state$status))persist_current("Process running")
      }
      invisible(NULL)
    },
    recover_project=function(project_config){
      state$project_config<-project_config;job<-read_json_safely(current_job_path());if(is.null(job))return(FALSE)
      state$job_id<-job$job_id%||%"";state$stage<-job$stage%||%"";state$status<-job$status%||%"unknown";state$pid<-as.integer(job$pid%||%NA);state$progress_file<-job$progress_file%||%"";state$selected_items<-unlist(job$selected_items%||%character());state$display_command<-job$display_command%||%"";state$recovered<-TRUE;state$dismissed<-FALSE
      p<-raw_progress();alive<-job_process_exists(state$pid);project_root<-dirname(dirname(project_config));valid<-job_stage_has_valid_output(project_root,state$stage);cache_reused<-grepl("cache",as.character(p$message%||%job$message%||%""),ignore.case=TRUE)
      reconciled<-reconcile_job_status(state$status,p$status%||%"",alive,valid,cache_reused,if(is.null(p))TRUE else job_heartbeat_is_stale(p))
      if(!identical(state$status,reconciled)){state$status<-reconciled;try(persist_current(if(identical(reconciled,"stale"))"Recorded process is no longer alive"else"Recovered reconciled job state"),silent=TRUE)}else if(job_status_is_terminal(reconciled)&&(!identical(job$active,FALSE)||as.integer(job$schema_version%||%1L)<2L))try(persist_current("Marked terminal job inactive"),silent=TRUE)
      TRUE
    },
    cancel_job=function(){if(!is.null(state$process)&&state$process$is_alive()){p<-progress();if(!is.null(p)){p$cancel_requested<-TRUE;p$status<-"cancelled";p$message<-sprintf("Cancelled after %d of %d items. Completed items can be reused when the job resumes.",p$processed_items%||%0,p$total_items%||%0);atomic_json_write(state$progress_file,p)};state$status<-"cancelled";persist_current("Cancellation requested");state$process$interrupt();Sys.sleep(.2);if(state$process$is_alive())state$process$kill_tree();persist_current("Process cancelled")};invisible(NULL)},
    job_is_running=function(){p<-raw_progress();if(!is.null(p)&&job_status_is_terminal(p$status))return(FALSE);if(!is.null(state$process))return(state$process$is_alive()&&job_status_is_active(state$status));job_status_is_active(state$status)&&job_process_exists(state$pid)},job_exit_status=function()if(is.null(state$process)||state$process$is_alive())NA_integer_ else state$process$get_exit_status(),
    dismiss_job=function(){state$dismissed<-TRUE;invisible(NULL)},clear_display_state=function(){state$dismissed<-TRUE;state$stdout<-state$stderr<-character();invisible(NULL)},current_job=function()reactiveValuesToList(state),current_progress=progress,last_progress=raw_progress,read_job_stdout=function()state$stdout,read_job_stderr=function()state$stderr,current_stdout=function()tail(state$stdout,100),current_stderr=function()tail(state$stderr,100),command=function()state$display_command,args=function()state$args,stage=function()state$stage,elapsed=function()if(is.null(state$started))0 else as.numeric(difftime(Sys.time(),state$started,units="secs")),job_elapsed=function()if(is.null(state$started))0 else as.numeric(difftime(Sys.time(),state$started,units="secs"))
  )
}
