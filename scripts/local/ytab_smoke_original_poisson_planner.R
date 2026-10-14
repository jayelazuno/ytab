#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
repo_root<-normalizePath(args[[1L]],mustWork=TRUE)
library(shiny);library(bslib)
`%||%`<-function(x,y)if(is.null(x)||!length(x))y else x
for(name in c("ui_helpers.R","ui_components.R","plot_customization_helpers.R","plot_display_helpers.R",
              "project_discovery.R","screen_coverage.R","screen_coverage_active.R","mutant_screen.R","navigation.R"))
  source(file.path(repo_root,"app/shiny/R",name))
html<-as.character(screen_coverage_ui())
stopifnot(grepl("Screen Planner",html,fixed=TRUE),grepl("Methods &amp; References",html,fixed=TRUE),
  !grepl("Monte Carlo|Pool Contribution|sc_model_mode|sc_run|sc_allocation_mode|sc_gene_count",html),
  length(regmatches(html,gregexpr('class="tab-pane',html,fixed=TRUE))[[1L]])==3L)
project<-read_project_summary(file.path(repo_root,"output/projects/H2O2_screen_v1/config/project.yaml"),repo_root)
active<-function()project;active_project_path<-function()project$project_config
python_bin<-function()stop("Poisson must not request Python")
test_server<-screen_coverage_server;formals(test_server)<-alist(input=,output=,session=)
check_curve<-function(value){
  traces<-jsonlite::fromJSON(value,simplifyVector=FALSE)$x$data
  curves<-Filter(function(t)identical(t$mode,"lines")&&length(t$x)>=100L&&length(t$y)>=100L,traces)
  stopifnot(length(curves)>=1L,diff(range(as.numeric(unlist(curves[[1L]]$y))))>90)
}
withCallingHandlers(shiny::testServer(test_server,{
  session$setInputs(sc_target_preset="0.95")
  stopifnot(request()$k==13119,request()$model_mode=="poisson")
  check_curve(output$sc_plot_primary)
  summary<-as.character(output$sc_primary_summary$html)
  stopifnot(grepl("4,379",summary,fixed=TRUE),grepl("5,220",summary,fixed=TRUE),
    grepl("169,099",summary,fixed=TRUE),grepl("49,034",summary,fixed=TRUE))
  session$setInputs(sc_target_preset="0.99")
  stopifnot(request()$k==20167);check_curve(output$sc_plot_primary)
  session$setInputs(sc_target_preset="0.999")
  stopifnot(request()$k==30250);check_curve(output$sc_plot_primary)
  session$setInputs(sc_target_preset="custom",sc_target=97.5)
  stopifnot(request()$k==ceiling(-4379*log(.025)))
  session$elapse(1000);session$flushReact()
  stopifnot(is.null(job$process),!session$isClosed(),length(output$sc_poisson_target)>0)
}),warning=function(e)stop(conditionMessage(e),call.=FALSE))
cat("PASS: original parent-union summary, visible Poisson curve, dynamic targets and hidden advanced controls\n")
