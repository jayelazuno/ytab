#!/usr/bin/env Rscript
# Run > Mapping runs MapFastq, CreateHitFile and SummaryTable as one job (ytab_run_preprocessing.py). The merged
# subtabs are gone, old links open Mapping, and the sample resolution sends only samples that still need a stage.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
repo_root<-normalizePath(args[[1L]],mustWork=TRUE)
library(shiny);library(bslib)
`%||%`<-function(x,y)if(is.null(x)||!length(x))y else x
for(name in c("ui_helpers.R","ui_components.R","job_progress.R","sample_selector.R","preprocessing_status.R","navigation.R","ui_preprocessing.R"))
  if(file.exists(file.path(repo_root,"app/shiny/R",name)))source(file.path(repo_root,"app/shiny/R",name))
html<-as.character(run_module_ui(4L))
values<-regmatches(html,gregexpr('data-value="[a-z_]+"',html))[[1]]
stopifnot(setequal(unique(sub('data-value="([a-z_]+)"',"\\1",values)),c("samples","reference","mapping","profiles")),
          grepl("run_all_preprocess",html,fixed=TRUE),grepl("preview_preprocess",html,fixed=TRUE),grepl("mapping_status_table",html,fixed=TRUE),
          !grepl("run_create_hit_file|run_summary_table|hits_status_table|summary_status_table",html))
for(old in c("hit_files","summary_tables","status"))stopifnot(navigation_translate_legacy("run",old)$nested_value=="mapping")
stopifnot(resolve_continue_destination(list(),repo_root,NULL,reviewed=FALSE)$nested=="samples",
          keep_going_supported_stage("preprocessing"))
js<-readLines(file.path(repo_root,"app/shiny/www/ytab_nav_menus.js"))
stopifnot(!any(grepl('"hit_files"|"summary_tables"|\\["Status", "status"\\]',js)))
status<-data.frame(sample=c("done","hpc_hits","fresh","no_input"),included=TRUE,fastq=c("missing","not_required","success","missing"),
                   mapping=c("success","success","ready","blocked"),bam_index=c("success","not_local","blocked","blocked"),
                   hit_file=c("success","success","blocked","blocked"),summary=c("success","ready","blocked","blocked"),stringsAsFactors=FALSE)
r<-resolve_selected_preprocessing_samples(status$sample,data.frame(sample=status$sample),status,"preprocess")
stopifnot(setequal(r$already_complete,"done"),setequal(r$eligible,c("hpc_hits","fresh")),identical(r$blocked,"no_input"),grepl("FASTQ unavailable",r$reasons[["no_input"]]))
# The combined runner on a project whose samples are all done runs nothing and writes nothing to the project.
python<-unname(Sys.which("python"));config<-file.path(repo_root,"output/projects/Gale_2020/config/project.yaml")
if(nzchar(python)&&file.exists(config)){
  stamp<-tempfile();file.create(stamp);Sys.sleep(1.1);progress<-file.path(tempdir(),"pp_smoke.progress.json")
  out<-system2(python,c(file.path(repo_root,"scripts/local/ytab_run_preprocessing.py"),"--project-config",config,"--dry-run","--progress-file",progress),stdout=TRUE,stderr=TRUE)
  touched<-list.files(file.path(repo_root,"output/projects/Gale_2020"),recursive=TRUE,full.names=TRUE);touched<-touched[file.mtime(touched)>file.mtime(stamp)]
  p<-jsonlite::fromJSON(progress);stopifnot(!length(touched),all(p$items$status=="skipped"),any(grepl("Mapping: already done",out,fixed=TRUE)))
}
cat("PASS: Run > Mapping runs mapping, hit files and SummaryTable together; done samples are skipped untouched\n")
