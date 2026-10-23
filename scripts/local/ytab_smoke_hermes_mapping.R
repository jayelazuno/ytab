#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
root<-normalizePath(args[[1L]],mustWork=TRUE)
library(shiny);library(bslib)
`%||%`<-function(x,y)if(is.null(x)||!length(x))y else x
for(name in c("ui_helpers.R","ui_components.R","sample_selector.R","job_progress.R","ui_preprocessing.R"))
  source(file.path(root,"app/shiny/R",name))
html<-as.character(run_module_ui())
stopifnot(grepl('mapping_junction_mate',html,fixed=TRUE),grepl('mapping_mate_agreement',html,fixed=TRUE),
          grepl('Read 1',html,fixed=TRUE),grepl('Read 2',html,fixed=TRUE))
stopifnot(length(mapping_read_arguments())==0L,
  identical(mapping_read_arguments("R2","no"),c("--junction-mate","R2","--no-require-mate-alignment")),
  identical(mapping_read_arguments("R1","yes"),c("--junction-mate","R1","--require-mate-alignment")))
invisible(parse(file=file.path(root,"app/shiny/app.R")))
cat("PASS: Hermes mapping UI, project defaults, explicit mate/concordance arguments, app syntax\n")
