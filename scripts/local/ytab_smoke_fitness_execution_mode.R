#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=1L)stop("Provide project.yaml",call.=FALSE)
all<-commandArgs(trailingOnly=FALSE);self<-normalizePath(sub("^--file=","",grep("^--file=",all,value=TRUE)[[1]]),winslash="/");root<-normalizePath(file.path(dirname(self),"../.."),winslash="/")
`%||%`<-function(left,right)if(is.null(left)||!length(left))right else left;source(file.path(root,"app/shiny/R/fitness_design_state.R"))
ui<-paste(readLines(file.path(root,"app/shiny/R/ui_fitness.R"),warn=FALSE),collapse="\n");server<-paste(readLines(file.path(root,"app/shiny/app.R"),warn=FALSE),collapse="\n")
stopifnot(fitness_action_label("preview",FALSE,FALSE)=="Preview fitness analysis",fitness_action_label("run",FALSE,FALSE)=="Run fitness analysis",fitness_action_label("run",TRUE,FALSE)=="Use cached result",fitness_action_label("run",TRUE,TRUE)=="Rerun fitness analysis",grepl('if(mode=="preview")"--dry-run"',server,fixed=TRUE),grepl('if(force)"--force"',server,fixed=TRUE),grepl('if(keep)"--keep-going"',server,fixed=TRUE),grepl("fitness_classifier_controls",ui),grepl("raw per-sample SummaryTable",ui),grepl("CPM normalization inside",ui),!grepl("Included comparison IDs",ui,fixed=TRUE))
cat("PASS\n")
