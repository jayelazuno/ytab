#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=1L)stop("Provide project.yaml",call.=FALSE)
all<-commandArgs(trailingOnly=FALSE);self<-normalizePath(sub("^--file=","",grep("^--file=",all,value=TRUE)[[1]]),winslash="/");root<-normalizePath(file.path(dirname(self),"../.."),winslash="/")
`%||%`<-function(left,right)if(is.null(left)||!length(left))right else left
source(file.path(root,"app/shiny/R/project_discovery.R"),local=TRUE);source(file.path(root,"app/shiny/R/preprocessing_status.R"),local=TRUE);source(file.path(root,"app/shiny/R/process_helpers.R"),local=TRUE)
checked<-validate_project_config_path(args[[1]],root);stopifnot(checked$valid);project<-read_project_summary(checked$path,root);stopifnot(project$species=="glabrata")
manifest<-file.path(project$project_root,"config","reference_resolved.json");stopifnot(file.exists(manifest));r<-reference_readiness(project,root)
stopifnot(r$fasta_found,r$annotation_found,r$index_complete,r$runnable,identical(r$label,"ready"),next_preprocessing_stage(project,root)=="Mapping",startsWith(r$fasta,paste0(root,"/")))
status_args<-c(file.path(root,"scripts/local/ytab_project_status.py"),"--project-config",checked$path,"--show-next");result<-run_process_sync(locate_python_executable(),status_args,wd=root);stopifnot(result$status==0L);project<-read_project_summary(checked$path,root);stopifnot(reference_readiness(project,root)$runnable)
cat("PASS\n")
