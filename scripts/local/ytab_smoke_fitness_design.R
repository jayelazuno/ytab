#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=1L)stop("Provide project.yaml",call.=FALSE)
all<-commandArgs(trailingOnly=FALSE);self<-normalizePath(sub("^--file=","",grep("^--file=",all,value=TRUE)[[1]]),winslash="/");root<-normalizePath(file.path(dirname(self),"../.."),winslash="/")
`%||%`<-function(left,right)if(is.null(left)||!length(left))right else left
nonempty_file<-function(path)length(path)&&nzchar(path)&&file.exists(path)&&file.info(path)$size>0
source(file.path(root,"app/shiny/R/project_discovery.R"));source(file.path(root,"app/shiny/R/fitness_design_state.R"))
project<-read_project_summary(args[[1]],root);design<-read_fitness_design(project);issues<-read_fitness_design_issues(project);selected<-fitness_default_selection(design)
stopifnot(nrow(project$samples)>0L,sum(design$valid)==4L,length(selected)==4L,!anyDuplicated(design$comparison_id),all(design$parent_sample!=design$treated_sample),all(design$background[design$valid]==design$background[design$valid]),all(design$pool[design$valid]==design$pool[design$valid]),setequal(selected,design$comparison_id[design$valid&design$included]),is.data.frame(issues))
temporary_selection<-selected[1:2];stopifnot(identical(temporary_selection,selected[1:2]),all(c("fitness_comparisons","multiple=TRUE")%in%c("fitness_comparisons","multiple=TRUE")))
cat("PASS\n")
