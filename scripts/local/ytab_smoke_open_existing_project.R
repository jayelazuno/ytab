#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
all<-commandArgs(trailingOnly=FALSE);self<-normalizePath(sub("^--file=","",grep("^--file=",all,value=TRUE)[[1L]]),winslash="/");root<-normalizePath(file.path(dirname(self),"../.."),winslash="/");config<-normalizePath(args[[1L]],winslash="/")
options(shiny.fullstacktrace=TRUE);app_env<-new.env(parent=.GlobalEnv);source(file.path(root,"app/shiny/app.R"),local=app_env)
shiny::testServer(app_env$server,{
  stopifnot(state$values$view=="landing",is.null(state$values$project))
  session$setInputs(project_selector=config,open_project=1);session$flushReact();stopifnot(state$values$view=="workspace",identical(state$values$path,config),!is.null(state$values$project))
  for(top in c("overview","preprocessing","quality_control","essentiality","fitness","results_exports")){session$setInputs(workspace_tabs=top);session$flushReact();stopifnot(navigation_location()$active_top_tab==top)}
  nested<-list(preprocessing_tabs="hit_files",qc_tabs="diagnostic_files",essentiality_tabs="combined_summary",fitness_tabs="fitness_run")
  for(id in names(nested)){do.call(session$setInputs,setNames(list(nested[[id]]),id));session$flushReact();stopifnot(identical(navigation_location()[[app_env$navigation_nested_fields[[id]]]],nested[[id]]))}
  session$setInputs(switch_project=1);session$flushReact();stopifnot(state$values$view=="landing",is.null(state$values$path),is.null(state$values$project),is.null(jobs$current_progress()))
  session$setInputs(project_selector=config,open_project=2);session$flushReact();stopifnot(state$values$view=="workspace",identical(state$values$path,config))
})
cat("PASS\n")
