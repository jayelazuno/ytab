args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
root <- normalizePath(args[1]); fixture <- normalizePath(args[2])
library(shiny)
library(bslib)
`%||%` <- function(x,y)if(is.null(x)||!length(x))y else x
for(name in c("ui_helpers.R","ui_components.R","plot_customization_helpers.R","plot_display_helpers.R","project_discovery.R","screen_coverage.R")) source(file.path(root,"app/shiny/R",name))
project <- read_project_summary(file.path(root,"output/projects/Zn_toxicity_screen/config/project.yaml"),root)
stopifnot(grepl("Screen Coverage",as.character(screen_coverage_ui()),fixed=TRUE))
test_server <- screen_coverage_server
formals(test_server)<-alist(input=,output=,session=)
repo_root<-root;active<-function()project
active_project_path<-function()file.path(project$project_root,"config/project.yaml")
python_bin<-function()"/Users/joshuaayelazuno/miniforge3/envs/ytab-local/bin/python"
shiny::testServer(test_server,{
  catalogue<-sc_parent_catalogue(root,active_project_path())
  stopifnot(nrow(catalogue)==8L,!any(grepl("treated",catalogue$sample)))
  session$setInputs(sc_target_preset="0.95",sc_r_preset="1",sc_weights="1",sc_target=95,sc_r=1,sc_k=10000,
    sc_max_k=100000,sc_iterations=1000,sc_seed=2026,sc_basis="sites",sc_conditioning="physical",sc_universe="informative",
    sc_samples=catalogue$key,sc_classifier=file.path(root,"output/projects/H2O2_screen_v1/classifier/T102400/essentiality_predictions.T102400.csv"),sc_simulate=TRUE,
    sc_pmf_mode="planned",sc_depth=1,sc_diagnostic="lorenz",sc_per_plate=96,sc_plates=100,sc_desired=10000)
  session$flushReact()
  stopifnot(length(request()$selected_sources)==8L,is.null(request()$classifier))
  result(jsonlite::fromJSON(fixture));submitted(request());session$flushReact()
  for(id in c("theory","redundancy","pmf","comparison","simulation","rarefaction","diagnostic","alleles")) {
    value<-output[[paste0("sc_plot_",id)]]
    stopifnot(length(value)>0)
  }
  for(id in c("summary","inventory","status","plate_result","planned_interval","limitations"))stopifnot(length(output[[paste0("sc_",id)]])>0)
  for(id in c("targets","comparison_targets","rarefaction_targets","genes","pools","planned","simulation_table","alleles_table","planned_alleles"))stopifnot(length(output[[paste0("sc_",id)]])>0)
  saved<-result()
  session$setInputs(sc_grid="show",sc_plot_style="presentation",sc_text_size="large",sc_diagnostic="length")
  session$flushReact();stopifnot(identical(saved,result()),is.null(job$process))
  session$setInputs(sc_classifier="a different classifier that must have no effect")
  session$elapse(1000);session$flushReact();stopifnot(identical(saved,result()),is.null(job$process))
  stopifnot(length(output$sc_plot_diagnostic)>0)
  session$setInputs(sc_diagnostic="abundance");session$flushReact();stopifnot(length(output$sc_plot_diagnostic)>0)
  session$setInputs(sc_diagnostic="probability");session$flushReact();stopifnot(length(output$sc_plot_diagnostic)>0)
  # Exercise the real asynchronous worker, followed by automatic analytical
  # recomputation. Cosmetic changes above must never start that worker.
  session$setInputs(sc_simulate=FALSE,sc_run=1)
  session$flushReact()
  for(i in seq_len(600)) {
    Sys.sleep(.1);session$elapse(500);session$flushReact()
    if(is.null(job$process))break
  }
  stopifnot(is.null(job$process),grepl("complete",message()),result()$settings$k==10000,result()$n==4379,result()$library$annotated_genes==5220)
  session$setInputs(sc_k=12000)
  session$elapse(1000);session$flushReact()
  for(i in seq_len(600)) {
    Sys.sleep(.1);session$elapse(500);session$flushReact()
    if(is.null(job$process)&&result()$settings$k==12000)break
  }
  stopifnot(is.null(job$process),result()$settings$k==12000,!isTRUE(result()$settings$simulate))
})
cat("PASS: Screen Coverage UI, plots, tables, planner, asynchronous worker, automatic analytical updates and cosmetic isolation\n")
