#!/usr/bin/env Rscript
# Screen Planning -> Coverage Planner (Poisson / equal-gene benchmark) on the H2O2 project's parent libraries.
# The planner's summary must match ytab_screen_coverage_project_summary.py, and each target's colony count must be
# K_q = ceil(-N ln(1 - q)) with N the represented genes. Also checks the curve, the downloads and that the retired
# Monte Carlo pool planner is gone.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
repo_root<-normalizePath(args[[1L]],mustWork=TRUE)
library(shiny);library(bslib)
`%||%`<-function(x,y)if(is.null(x)||!length(x))y else x
for(name in c("job_manager.R","ui_helpers.R","ui_components.R","plot_customization_helpers.R","plot_display_helpers.R",
              "project_discovery.R","analysis_runs.R","screen_planning.R","mutant_screen.R","navigation.R"))
  source(file.path(repo_root,"app/shiny/R",name))
html<-as.character(screen_coverage_ui())
stopifnot(grepl("Coverage Planner",html,fixed=TRUE),grepl("Mutant Identification",html,fixed=TRUE),grepl("Knockout Sudoku",html,fixed=TRUE),
  grepl("sc_download_curve",html,fixed=TRUE),!grepl("Methods &amp; References",html,fixed=TRUE),
  !grepl("Monte Carlo|Pool Contribution|sc_model_mode|sc_run|sc_allocation_mode|sc_gene_count",html),
  length(regmatches(html,gregexpr('class="tab-pane',html,fixed=TRUE))[[1L]])==2L)
config<-file.path(repo_root,"output/projects/H2O2_screen_v1/config/project.yaml")
project<-read_project_summary(config,repo_root)
python<-unname(Sys.which("python"));expected_file<-tempfile(fileext=".json")
stopifnot(system2(python,c(file.path(repo_root,"scripts/local/ytab_screen_coverage_project_summary.py"),"--project-config",config,"--output",expected_file),stdout=FALSE,stderr=FALSE)==0L)
expected<-jsonlite::fromJSON(expected_file);n<-expected$represented
fmt<-function(x)format(round(x),big.mark=",",scientific=FALSE,trim=TRUE)
server<-function(input,output,session)screen_coverage_server(input,output,session,repo_root,function()project,function()project$project_config,function()python)
shiny::testServer(server,{
  session$setInputs(sc_target_preset="0.95");deadline<-Sys.time()+120
  repeat{session$elapse(350);Sys.sleep(.05);session$flushReact()
    summary<-tryCatch(output$sc_primary_summary,shiny.silent.error=function(e)NULL);if(!is.null(summary))break
    if(Sys.time()>deadline)stop("Timed out loading the project's parent-library summary")}
  summary<-paste(unlist(summary),collapse=" ")
  stopifnot(grepl(fmt(expected$annotated),summary,fixed=TRUE),grepl(fmt(n),summary,fixed=TRUE),
            grepl(fmt(expected$unique_sites),summary,fixed=TRUE),grepl(fmt(expected$coding_sites),summary,fixed=TRUE))
  for(q in c(.95,.99,.999)){
    session$setInputs(sc_target_preset=as.character(q))
    stopifnot(grepl(paste(fmt(ceiling(-n*log1p(-q))),"colonies required"),paste(unlist(output$sc_poisson_target),collapse=" "),fixed=TRUE))
    traces<-jsonlite::fromJSON(output$sc_plot_primary,simplifyVector=FALSE)$x$data
    curves<-Filter(function(t)identical(t$mode,"lines")&&length(t$x)>=100L,traces);stopifnot(length(curves)>=1L)
  }
  session$setInputs(sc_target_preset="custom",sc_target=97.5)
  stopifnot(grepl(paste(fmt(ceiling(-n*log(.025))),"colonies required"),paste(unlist(output$sc_poisson_target),collapse=" "),fixed=TRUE))
  for(id in c("sc_download_png","sc_download_pdf","sc_download_curve"))stopifnot(file.size(output[[id]])>100)
  curve<-utils::read.delim(output$sc_download_curve);stopifnot(nrow(curve)==201L,all(curve$represented_genes==n))
})
stopifnot(navigation_target_from_hash("#tab=screen_coverage&subtab=methods")$nested_value=="methods",
          normalize_navigation_state(list(active_screen_coverage_subtab="methods"))$active_screen_coverage_subtab=="planner")
cat("PASS: Coverage Planner summary from the project's parents, K_q = ceil(-N ln(1-q)), curve, downloads; legacy planner retired\n")
