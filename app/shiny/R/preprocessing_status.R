resolve_reference_path <- function(path, repo_root) {
  if (is.null(path) || !length(path) || !nzchar(path)) return("")
  ytab_resolve_path(path, repo_root)
}

reference_readiness <- function(project, repo_root) {
  ref <- project$reference %||% list()
  fasta <- resolve_reference_path(ref$fasta %||% "", repo_root)
  annotation_values <- c(ref$feature_table %||% "", ref$gff %||% "", ref$gtf %||% "")
  annotations <- vapply(annotation_values[nzchar(annotation_values)], resolve_reference_path, "", repo_root=repo_root)
  prefix <- resolve_reference_path(ref$bowtie2_index_prefix %||% "", repo_root)
  index_complete <- nzchar(prefix) && any(vapply(c("bt2","bt2l"), function(ext)
    all(file.exists(paste0(prefix,c(".1.",".2.",".3.",".4.",".rev.1.",".rev.2."),ext))), FALSE))
  fasta_found <- nzchar(fasta) && file.exists(fasta) && !dir.exists(fasta)
  annotation_found <- length(annotations) > 0L && any(file.exists(annotations))
  runnable <- fasta_found && annotation_found && index_complete
  list(runnable=runnable, can_prepare=fasta_found && annotation_found, fasta=fasta,
    fasta_found=fasta_found, annotations=annotations, annotation_found=annotation_found,
    index_prefix=prefix, index_complete=index_complete,
    label=if(runnable)"ready" else if(fasta_found)"preparation required" else "files required")
}

read_stage_rows <- function(project_root, stage) {
  rel <- c(mapping="mapfastq/mapfastq_status.csv",hit_file="create_hit_file/create_hit_file_status.csv",summary="summary/summary_status.csv")[[stage]]
  path <- file.path(project_root,"manifests",rel)
  if (!file.exists(path)) return(data.frame())
  tryCatch(read.csv(path,stringsAsFactors=FALSE,check.names=FALSE),error=function(e)data.frame())
}

complete_stage_statuses <- c("success","skipped","cached","imported_success","imported_or_not_required")

project_starts_from_hit_files <- function(project) {
  identical(tolower(as.character(project$start_stage %||% "fastq")), "create_hit_file")
}

manifest_sample_status <- function(rows, sample, prerequisite=TRUE) {
  if(!prerequisite)return("blocked")
  if(!nrow(rows)||!"sample"%in%names(rows)||!sample%in%rows$sample)return("ready")
  hit<-tail(rows[rows$sample==sample,,drop=FALSE],1);value<-tolower(as.character(hit$status[[1]]))
  message<-if("message"%in%names(hit))tolower(as.character(hit$message[[1]]))else""
  if(value=="skipped"&&grepl("dry run",message,fixed=TRUE))return("ready")
  if(value%in%c(complete_stage_statuses,"failed","running","stale","not_required"))value else "not started"
}

nonempty_file <- function(path) nzchar(path %||% "") && file.exists(path) && !dir.exists(path) && isTRUE(file.info(path)$size > 0)

scientific_status_row <- function(rows,sample) {
  if(!nrow(rows)||!"sample"%in%names(rows)||!sample%in%rows$sample)return(NULL)
  candidates<-rows[rows$sample==sample,,drop=FALSE]
  messages<-if("message"%in%names(candidates))tolower(as.character(candidates$message))else rep("",nrow(candidates))
  dry<-grepl("dry run; command was not executed",messages,fixed=TRUE)
  if("dry_run"%in%names(candidates))dry<-dry|tolower(as.character(candidates$dry_run))%in%c("true","1","yes")
  if("status"%in%names(candidates))dry<-dry|tolower(as.character(candidates$status))=="dry_run_success"
  candidates<-candidates[!dry,,drop=FALSE];if(!nrow(candidates))NULL else tail(candidates,1)
}

row_value <- function(row,column,default="")if(!is.null(row)&&column%in%names(row))as.character(row[[column]][[1]])else default

sample_value <- function(samples, i, column, default="") {
  if(column %in% names(samples)) as.character(samples[[column]][[i]]) else default
}

resolve_project_file <- function(path, repo_root) {
  if(!nzchar(path %||% ""))return("")
  ytab_resolve_path(path, repo_root)
}

build_sample_pipeline_status <- function(project, samples=project$samples, project_status=project$status) {
  if(!nrow(samples))return(data.frame())
  map_rows<-read_stage_rows(project$project_root,"mapping");hit_rows<-read_stage_rows(project$project_root,"hit_file");sum_rows<-read_stage_rows(project$project_root,"summary")
  included<-if("include"%in%names(samples))tolower(as.character(samples$include))%in%c("true","1","yes")else rep(TRUE,nrow(samples))
  result<-lapply(seq_len(nrow(samples)),function(i){
    s<-as.character(samples$sample[[i]]);map_row<-scientific_status_row(map_rows,s);hit_row<-scientific_status_row(hit_rows,s);sum_row<-scientific_status_row(sum_rows,s)
    imported<-project_starts_from_hit_files(project);repo_root<-project$repo_root %||% normalizePath(file.path(project$project_root,"../../.."),winslash="/",mustWork=FALSE)
    fastq_1<-sample_value(samples,i,"fastq_1","")
    sample_hit<-sample_value(samples,i,"hit_file","")
    if(nzchar(sample_hit))sample_hit<-resolve_project_file(sample_hit,repo_root)
    bam<-file.path(project$project_root,"mapfastq",s,paste0(s,".sorted.bam"));bai<-c(paste0(bam,".bai"),sub("\\.bam$",".bai",bam));stats<-file.path(project$project_root,"mapfastq",s,paste0(s,".mapping_stats.csv"));hits<-if(nzchar(sample_hit))sample_hit else file.path(project$project_root,"create_hit_file",s,paste0(s,"_hits.txt"));feature_files<-if(dir.exists(file.path(project$project_root,"summary",s)))list.files(file.path(project$project_root,"summary",s),pattern="\\.feature_table.*\\.(csv|tsv|txt)$",full.names=TRUE)else character()
    mapped_ok<-nonempty_file(bam)&&any(vapply(bai,nonempty_file,FALSE));hits_ok<-nonempty_file(hits);summary_ok<-length(feature_files)>0&&any(vapply(feature_files,nonempty_file,FALSE))
    map_value<-tolower(row_value(map_row,"status",""));hit_value<-tolower(row_value(hit_row,"status",""));sum_value<-tolower(row_value(sum_row,"status",""))
    mapping<-if(map_value%in%c("failed","running"))map_value else if(imported)"imported_or_not_required"else if(mapped_ok&&nonempty_file(stats))"success"else if(nonempty_file(fastq_1))"ready"else"blocked"
    hit<-if(hit_value%in%c("failed","running"))hit_value else if(imported&&hits_ok)"imported_success"else if(hits_ok)"success"else if(mapped_ok)"ready"else"blocked"
    summary<-if(sum_value%in%c("failed","running"))sum_value else if(summary_ok)"success"else if(hits_ok)"ready"else"blocked"
    data.frame(sample=s,included=included[[i]],fastq=if(imported)"not_required"else if(nonempty_file(fastq_1))"success"else"missing",mapping=mapping,bam_index=if(imported)"not_required"else if(any(vapply(bai,nonempty_file,FALSE)))"success"else if(nonempty_file(bam))"ready"else"blocked",mapped_bam=bam,hit_file=hit,hits_path=hits,summary=summary,feature_table=if(length(feature_files))paste(feature_files,collapse=";")else"",mapping_elapsed=row_value(map_row,"elapsed_seconds"),hit_elapsed=row_value(hit_row,"elapsed_seconds"),summary_elapsed=row_value(sum_row,"elapsed_seconds"),mapping_message=if(imported)"MapFastq is not required for imported hit-file projects."else if(mapping=="success")"Scientific mapping outputs available."else row_value(map_row,"message"),hit_message=if(imported&&hits_ok)"Existing CreateHitFile output imported."else if(hit=="ready")"Mapped BAM and index available; hit-file creation is ready."else row_value(hit_row,"message"),summary_message=if(summary=="ready")"Hit file available; SummaryTable is ready."else row_value(sum_row,"message"),stringsAsFactors=FALSE)
  })
  do.call(rbind,result)
}

resolve_selected_preprocessing_samples <- function(selected_samples,sample_sheet,stage_status,stage,force=FALSE) {
  requested<-unique(as.character(selected_samples));known<-as.character(sample_sheet$sample);missing<-setdiff(requested,known)
  included_names<-as.character(stage_status$sample[stage_status$included]);included<-intersect(requested,included_names);rows<-stage_status[match(included,stage_status$sample),,drop=FALSE]
  column<-c(mapping="mapping",create_hit_file="hit_file",summary="summary")[[stage]];prereq<-switch(stage,mapping=rows$fastq=="success",create_hit_file=rows$mapping%in%c("success","imported_or_not_required")&rows$bam_index%in%c("success","not_required"),summary=rows$hit_file%in%c("success","imported_success"))
  complete<-rows[[column]]%in%complete_stage_statuses;eligible<-included[prereq&(force|!complete)];already_complete<-included[complete&!force];blocked<-included[!prereq]
  reasons<-setNames(if(stage=="create_hit_file")rep("mapped BAM or BAM index unavailable",length(blocked))else if(stage=="summary")rep("hit file unavailable",length(blocked))else rep("FASTQ or reference unavailable",length(blocked)),blocked)
  warnings<-c(if(length(missing))paste("Unknown samples:",paste(missing,collapse=", ")),if(length(blocked))paste("Blocked:",paste(sprintf("%s (%s)",blocked,reasons),collapse=", ")))
  list(requested=requested,included=included,eligible=eligible,already_complete=already_complete,blocked=blocked,missing=missing,reasons=reasons,warnings=warnings)
}

stage_counts <- function(status, column, included_only=TRUE) {
  rows<-if(included_only)status[status$included,,drop=FALSE]else status; values<-rows[[column]]
  list(total=nrow(rows),complete=sum(values%in%complete_stage_statuses),failed=sum(values=="failed"),remaining=sum(!values%in%complete_stage_statuses))
}

next_preprocessing_stage <- function(project, repo_root, status=build_sample_pipeline_status(project)) {
  ref<-reference_readiness(project,repo_root);included<-status[status$included,,drop=FALSE]
  if(!ref$runnable && !(project_starts_from_hit_files(project) && ref$can_prepare))return("Reference")
  if(any(!included$mapping%in%complete_stage_statuses))return("Mapping")
  if(any(!included$hit_file%in%complete_stage_statuses))return("Hit Files")
  if(any(!included$summary%in%complete_stage_statuses))return("Summary Tables")
  "Complete"
}
