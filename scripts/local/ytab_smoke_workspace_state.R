#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=1L)stop("Provide project.yaml",call.=FALSE)
all<-commandArgs(trailingOnly=FALSE);self<-normalizePath(sub("^--file=","",grep("^--file=",all,value=TRUE)[[1]]),winslash="/");root<-normalizePath(file.path(dirname(self),"../.."),winslash="/")
app_env<-new.env(parent=.GlobalEnv);source(file.path(root,"app/shiny/app.R"),local=app_env)
shiny::testServer(app_env$server,{session$setInputs(project_selector=normalizePath(args[[1]],winslash="/"));session$setInputs(open_project=1);session$flushReact();stopifnot(identical(state$values$view,"workspace"));stopifnot(nzchar(isolate(state$values$workspace_tab)))})
cat("PASS\n")
