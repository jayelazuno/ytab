read_diagnostics_run_metadata <- function(path, project_root, legacy=FALSE) {
  manifest<-tryCatch(jsonlite::fromJSON(path,simplifyVector=FALSE),error=function(e)NULL);if(is.null(manifest))return(NULL)
  selected<-as.character(unlist(manifest$selected_samples%||%character()));inputs<-as.character(unlist(manifest$input_files%||%manifest$input_hit_files%||%character()))
  if(!length(selected)&&length(inputs))selected<-unique(sub("_hits\\.txt$","",basename(inputs)))
  output<-as.character(manifest$output_dir%||%file.path(project_root,"library_diagnostics"));detected<-as.character(unlist(manifest$detected_outputs%||%character()));count<-if(length(detected))sum(file.exists(detected))else if(dir.exists(output))length(list.files(output,recursive=TRUE,full.names=TRUE))else 0L
  run_id<-as.character(manifest$run_id%||%if(legacy)"current_legacy"else basename(dirname(dirname(path))));label<-as.character(manifest$sample_set_label%||%if(legacy)"Legacy diagnostics result"else tools::toTitleCase(gsub("_"," ",run_id)))
  list(run_id=run_id,label=label,selected_samples=selected,sample_count=length(selected),cache_signature=as.character(manifest$cache_signature%||%""),status=as.character(manifest$status%||%"unknown"),completed=as.character(manifest$end_time%||%manifest$completed_at%||%"Unknown"),output_count=as.integer(count),output_dir=output,manifest_path=path,is_legacy=isTRUE(legacy)||!nzchar(as.character(manifest$cache_signature%||%"")),metadata_complete=length(selected)>0L,elapsed_seconds=as.numeric(manifest$elapsed_seconds%||%NA_real_),species=as.character(manifest$species%||%""),input_hashes=unlist(manifest$input_hashes%||%manifest$input_hit_file_hashes%||%list()),reference_hashes=unlist(manifest$reference_hashes%||%list()),script_hash=as.character(manifest$library_diagnostics_script_hash%||%""),scientific_parameters=manifest$scientific_parameters%||%list(),project_root=project_root)
}

validate_diagnostics_run_inputs <- function(run,current_selection,project) {
  if(is.null(run)||isTRUE(run$is_legacy)||!nzchar(run$cache_signature)||!setequal(run$selected_samples,current_selection)||!identical(run$species,as.character(project$species)))return(FALSE)
  check_hashes<-function(values){if(!length(values)||is.null(names(values))||any(!nzchar(names(values))))return(FALSE);all(vapply(seq_along(values),function(i){path<-names(values)[[i]];file.exists(path)&&!dir.exists(path)&&identical(unname(digest::digest(file=path,algo="sha256",serialize=FALSE)),unname(as.character(values[[i]])))},FALSE))}
  if(!check_hashes(run$input_hashes)||!check_hashes(run$reference_hashes))return(FALSE)
  repo_root<-dirname(dirname(dirname(project$project_root)));script<-file.path(repo_root,"src","ytab","qc","LibraryDiagnostics.py");if(!file.exists(script)||!identical(digest::digest(file=script,algo="sha256",serialize=FALSE),run$script_hash))return(FALSE)
  stored_threads<-run$scientific_parameters$threads%||%NULL;is.null(stored_threads)||identical(as.integer(stored_threads),as.integer(project$threads))
}

discover_diagnostics_runs <- function(project) {
  base<-file.path(project$project_root,"manifests","library_diagnostics");paths<-if(dir.exists(file.path(base,"runs")))list.files(file.path(base,"runs"),pattern="manifest\\.json$",recursive=TRUE,full.names=TRUE)else character();runs<-lapply(paths,read_diagnostics_run_metadata,project_root=project$project_root,legacy=FALSE)
  legacy<-file.path(base,"library_diagnostics_manifest.json");if(file.exists(legacy))runs<-c(runs,list(read_diagnostics_run_metadata(legacy,project$project_root,legacy=TRUE)))
  Filter(function(x)!is.null(x)&&x$status%in%c("success","cached","skipped"),runs)
}

choose_latest_diagnostics_result <- function(runs) {
  if(!length(runs))return(NULL);times<-vapply(runs,function(x){value<-tryCatch(as.POSIXct(x$completed,tz="UTC"),error=function(e)as.POSIXct(NA));if(is.na(value))0 else as.numeric(value)},0);runs[[which.max(times)]]
}
choose_matching_diagnostics_result <- function(runs,current_selection,current_cache_signature="") {
  if(!length(current_selection)||!nzchar(current_cache_signature))return(NULL)
  candidates<-Filter(function(x)!x$is_legacy&&nzchar(x$cache_signature)&&identical(x$cache_signature,current_cache_signature)&&setequal(x$selected_samples,current_selection),runs)
  choose_latest_diagnostics_result(candidates)
}
format_diagnostics_sample_set <- function(selection,eligible=data.frame()) {
  selection<-as.character(selection);if(!length(selection))return("No selection")
  if(nrow(eligible)&&setequal(selection,eligible$Sample[eligible$Eligible=="Yes"]))return("All eligible samples")
  condition<-if(nrow(eligible))tolower(eligible$Condition[match(selection,eligible$Sample)])else character()
  if(length(condition)&&all(condition=="parent"))"Parents"else if(length(condition)&&all(condition=="treated"))"Treated"else"Custom selection"
}
resolve_qc_result_state <- function(current_selection,current_cache_signature="",available_runs=list(),eligible_count=0L,current_label=NULL) {
  current_selection<-as.character(current_selection);matching<-choose_matching_diagnostics_result(available_runs,current_selection,current_cache_signature);latest<-choose_latest_diagnostics_result(available_runs);count<-length(current_selection);label<-current_label%||%if(count)"Custom selection"else"No selection"
  status<-if(!count)"no_selection"else if(!length(available_runs))"no_result"else if(!is.null(matching)&&matching$status=="cached")"matching_cached_result"else if(!is.null(matching))"matching_completed_result"else"historical_result_only"
  message<-switch(status,no_selection="No samples are currently selected.",no_result="No completed diagnostics results are available for this project.",matching_cached_result="A matching successful diagnostics result exists for the current sample selection.",matching_completed_result="A matching completed diagnostics result exists for the current sample selection.",historical_result_only="The latest available diagnostics were generated from a different sample selection.")
  list(current_selection_count=count,current_selection_label=label,eligible_count=as.integer(eligible_count),matching_result=matching,latest_historical_result=latest,match_status=status,message=message)
}

build_diagnostics_result_summary <- function(run) {
  if(is.null(run))return(NULL);list(sample_set=if(run$is_legacy)"Legacy diagnostics result"else run$label,samples_analyzed=if(run$metadata_complete)run$sample_count else NA_integer_,output_files=run$output_count,completed=run$completed,run_id=run$run_id,is_legacy=run$is_legacy,warning=if(run$is_legacy&&!run$metadata_complete)"Detailed sample-set metadata were not recorded for this legacy result."else"")
}
