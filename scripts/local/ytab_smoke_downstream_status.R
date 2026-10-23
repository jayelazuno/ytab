#!/usr/bin/env Rscript
# A stage whose own outputs are absent counts as done when a downstream output exists (e.g. mapping ran on the HPC and
# only hit files or summary tables were copied back). Synthetic project: one sample per evidence pattern.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
repo_root<-normalizePath(args[[1L]],mustWork=TRUE)
library(shiny);library(bslib)
`%||%`<-function(x,y)if(is.null(x)||!length(x))y else x
for(name in c("ui_helpers.R","ui_components.R","reference_readiness.R","preprocessing_status.R","sample_selector.R","navigation.R","ui_workspace.R"))
  if(file.exists(file.path(repo_root,"app/shiny/R",name)))source(file.path(repo_root,"app/shiny/R",name))
root<-tempfile("downstream_status_");dir.create(root)
put<-function(...){path<-file.path(root,...);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);writeLines("x",path)}
fastq<-file.path(root,"fastq/raw.fastq.gz");put("fastq/raw.fastq.gz")
put("create_hit_file/hits_only/hits_only_hits.txt")                      # hit file, no BAM
put("summary/summary_only/summary_only.feature_table.csv")             # summary, no BAM, no hit file
samples<-data.frame(sample=c("hits_only","summary_only","fastq_only","nothing"),include="true",
                    fastq_1=c("","","",""),stringsAsFactors=FALSE);samples$fastq_1[3]<-fastq
project<-list(project_root=root,repo_root=repo_root,samples=samples,start_stage="fastq")
s<-build_sample_pipeline_status(project);rownames(s)<-s$sample
stopifnot(s["hits_only","mapping"]=="success",s["hits_only","hit_file"]=="success",s["hits_only","fastq"]=="not_required",s["hits_only","bam_index"]=="not_local",
          s["summary_only","mapping"]=="success",s["summary_only","hit_file"]=="success",s["summary_only","summary"]=="success",
          s["fastq_only","mapping"]=="ready",s["fastq_only","hit_file"]=="blocked",
          s["nothing","mapping"]=="blocked",grepl("FASTQ file not found",s["nothing","mapping_message"]),s["nothing","fastq"]=="missing",
          all(c(s["hits_only","mapping"],s["summary_only","hit_file"])%in%complete_stage_statuses),
          grepl("elsewhere",s["hits_only","mapping_message"]),grepl("SummaryTable exists",s["summary_only","hit_message"]))
# Done samples are reported as done, not blocked; forcing a rerun without the BAM is blocked.
sel<-resolve_selected_preprocessing_samples(c("hits_only","summary_only"),samples,s,"create_hit_file")
stopifnot(setequal(sel$already_complete,c("hits_only","summary_only")),!length(sel$blocked),!length(sel$eligible))
forced<-resolve_selected_preprocessing_samples("hits_only",samples,s,"create_hit_file",force=TRUE);stopifnot(identical(forced$blocked,"hits_only"))
# Project > Overview is retired; old links land on Samples.
html<-as.character(project_module_ui())
stopifnot(!grepl('data-value="overview"',html,fixed=TRUE),grepl('data-value="samples"',html,fixed=TRUE),
          navigation_translate_legacy("project","overview")$nested_value=="samples",navigation_state_defaults("project")$active_project_subtab=="samples")
# Real data: Gale 2020 was mapped on the HPC; only hit files and summaries are local.
gale<-file.path(repo_root,"output/projects/Gale_2020/config/project.yaml")
if(file.exists(gale)&&exists("read_project_summary")){g<-build_sample_pipeline_status(read_project_summary(gale,repo_root));stopifnot(all(g$mapping[g$included]%in%complete_stage_statuses))}
cat("PASS: downstream outputs mark upstream stages success; Overview subtab retired\n")
