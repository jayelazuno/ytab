#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript scripts/local/ytab_smoke_app_shell.R <project.yaml>",call.=FALSE)
all_args<-commandArgs(trailingOnly=FALSE);file_arg<-grep("^--file=",all_args,value=TRUE)
script<-normalizePath(sub("^--file=","",file_arg[[1]]),winslash="/",mustWork=TRUE)
repo_root<-normalizePath(file.path(dirname(script),"../.."),winslash="/",mustWork=TRUE)
`%||%`<-function(left,right)if(is.null(left)||!length(left))right else left
source(file.path(repo_root,"app/shiny/R/project_discovery.R"),local=TRUE)
failures<-character();check<-function(ok,label){if(!isTRUE(ok))failures<<-c(failures,label)}
config<-ytab_resolve_path(args[[1]],repo_root);valid<-validate_project_config_path(config,repo_root)
projects<-discover_ytab_projects(repo_root);info<-if(valid$valid)read_project_summary(config,repo_root)else NULL
check(identical(repo_root,normalizePath(file.path(dirname(script),"../.."),winslash="/")),"repository root resolves")
check(nrow(projects)>0,"existing projects discovered");check(any(grepl("smoke|demo",projects$project_id,ignore.case=TRUE)),"smoke project discovered")
check(config%in%projects$project_config,"exact project.yaml returned");check(!validate_project_config_path(repo_root,repo_root)$valid,"repository directory rejected")
check(!validate_project_config_path(dirname(dirname(config)),repo_root)$valid,"project directory rejected");check(valid$valid,"valid project.yaml accepted")
check(nzchar(info$project_id),"project ID loaded");check(nzchar(info$species),"species loaded");check(nzchar(info$fastq_directory),"FASTQ directory loaded")
check(file.exists(info$sample_sheet),"sample sheet resolves");check(info$included_count>0,"included samples present");check(is.logical(info$analysis_ready),"analysis readiness calculated")
check(!any(projects$project_config==repo_root),"repository root is never a config")
# the same project under a repository path containing a space resolves to that exact path
spaced<-file.path(tempfile("ytab smoke "),"repo with space");spaced_config<-file.path(spaced,"output/projects",basename(dirname(dirname(config))),"config")
dir.create(spaced_config,recursive=TRUE);on.exit(unlink(dirname(spaced),recursive=TRUE),add=TRUE)
invisible(file.copy(c(config,info$sample_sheet),spaced_config));spaced_root<-normalizePath(spaced,winslash="/")
spaced_valid<-validate_project_config_path(file.path(spaced_config,"project.yaml"),spaced_root)
check(isTRUE(spaced_valid$valid)&&grepl(" ",spaced_valid$path,fixed=TRUE)&&file.exists(spaced_valid$sample_sheet),"paths containing spaces remain intact")
unlink(dirname(spaced),recursive=TRUE)
if(length(failures)){cat("FAIL\n",paste("-",failures,collapse="\n"),"\n",sep="");quit(status=1)}
cat("PASS\n")
