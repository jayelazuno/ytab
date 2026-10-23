#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=1L)stop("Provide project.yaml",call.=FALSE)
all<-commandArgs(trailingOnly=FALSE);self<-normalizePath(sub("^--file=","",grep("^--file=",all,value=TRUE)[[1]]),winslash="/");root<-normalizePath(file.path(dirname(self),"../.."),winslash="/")
`%||%`<-function(left,right)if(is.null(left)||!length(left))right else left;source(file.path(root,"app/shiny/R/fitness_design_state.R"))
ui<-paste(readLines(file.path(root,"app/shiny/R/ui_fitness.R"),warn=FALSE),collapse="\n");server<-paste(readLines(file.path(root,"app/shiny/app.R"),warn=FALSE),collapse="\n")
stopifnot(fitness_action_label("preview",FALSE,FALSE)=="Preview z-score analysis",fitness_action_label("run",FALSE,FALSE)=="Run z-score analysis",fitness_action_label("run",TRUE,FALSE)=="Run z-score analysis",fitness_action_label("run",TRUE,TRUE)=="Rerun z-score analysis",
  grepl('actionButton("fitness_open_match","Open existing result"',server,fixed=TRUE),grepl('actionButton("fitness_open_existing","Open selected result")',server,fixed=TRUE),!grepl("CPM inside R",server,fixed=TRUE),grepl('if(mode=="preview")"--dry-run"',server,fixed=TRUE),grepl('if(force)"--force"',server,fixed=TRUE),grepl('if(keep)"--keep-going"',server,fixed=TRUE),grepl("fitness_classifier_controls",ui),grepl("Total reads (CPM) needs no normalization run",ui,fixed=TRUE),grepl("Total reads (CPM)",ui,fixed=TRUE),!grepl("inside the R analysis",ui,fixed=TRUE),!grepl("Included comparison IDs",ui,fixed=TRUE))
cat("PASS\n")
