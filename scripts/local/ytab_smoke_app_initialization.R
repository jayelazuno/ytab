#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE); execute <- "--execute" %in% args; args <- setdiff(args,"--execute")
if(length(args)!=1L)stop("Usage: Rscript scripts/local/ytab_smoke_app_initialization.R <project.yaml> [--execute]",call.=FALSE)
all_args<-commandArgs(trailingOnly=FALSE);file_arg<-grep("^--file=",all_args,value=TRUE);self<-normalizePath(sub("^--file=","",file_arg[[1]]),winslash="/",mustWork=TRUE)
repo_root<-normalizePath(file.path(dirname(self),"../.."),winslash="/",mustWork=TRUE);`%||%`<-function(left,right)if(is.null(left)||!length(left))right else left
source(file.path(repo_root,"app/shiny/R/project_discovery.R"),local=TRUE);source(file.path(repo_root,"app/shiny/R/process_helpers.R"),local=TRUE)
source_info<-read_project_summary(args[[1]],repo_root);fastq_dir<-normalize_local_directory(source_info$fastq_directory)
project_id<-if(execute)sprintf("app_init_command_smoke_%s",format(Sys.time(),"%Y%m%d%H%M%S"))else"app_init_command_smoke"
init_script<-file.path(repo_root,"scripts/local/ytab_init_local_project.py")
init_args<-c(init_script,"--project-id",project_id,"--fastq-dir",fastq_dir,"--species","glabrata","--threads","2","--repo-root",repo_root)
python<-locate_python_executable();preview<-format_command_for_display(python,init_args);cat(preview,"\n")
stopifnot(file.exists(init_script),!dir.exists(init_script),dir.exists(fastq_dir),length(init_args)==11L,init_args[[1]]==init_script,init_args[[5]]==fastq_dir,init_args[[11]]==repo_root,!any(init_args%in%c("2>&1",">","|","&&",";")))
if(execute){target<-file.path(repo_root,"output/projects",project_id);stopifnot(!file.exists(target));result<-run_process_sync(python,init_args,wd=repo_root);if(result$status!=0L)stop(result$stderr,call.=FALSE);config<-file.path(target,"config/project.yaml");stopifnot(validate_project_config_path(config,repo_root)$valid);cat("Created:",config,"\n")}
cat("PASS\n")
