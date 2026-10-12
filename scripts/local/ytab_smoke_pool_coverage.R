args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);fixture<-normalizePath(args[2])
library(shiny);library(bslib)
`%||%`<-function(x,y)if(is.null(x)||!length(x))y else x
for(name in c("ui_helpers.R","ui_components.R","plot_customization_helpers.R","plot_display_helpers.R","project_discovery.R","screen_coverage.R"))source(file.path(root,"app/shiny/R",name))
project<-read_project_summary(file.path(root,"output/projects/Zn_toxicity_screen/config/project.yaml"),root)
html<-as.character(screen_coverage_ui())
stopifnot(all(vapply(c("Screen Planner","Pool Contribution / QC","Library QC","Parent Sequencing","Statistical Details","Math","Methods"),function(x)grepl(x,html,fixed=TRUE),logical(1))))
stopifnot(grepl("Coverage model",html,fixed=TRUE),!grepl("Total planned colonies K",html,fixed=TRUE))
test_server<-screen_coverage_server;formals(test_server)<-alist(input=,output=,session=)
repo_root<-root;active<-function()project;active_project_path<-function()file.path(project$project_root,"config/project.yaml")
python_bin<-function()"/Users/joshuaayelazuno/miniforge3/bin/python"
screen_runner<-sc_python_command(python_bin,root)
stopifnot(grepl("conda$",screen_runner$command),identical(tail(screen_runner$prefix,1),"python"))
shiny::testServer(test_server,{
  catalogue<-sc_parent_catalogue(root,active_project_path())
  stopifnot(nrow(catalogue)==8,length(unique(catalogue$physical_pool_id))==4)
  saved_fixture<-jsonlite::fromJSON(fixture)
  stopifnot(setequal(catalogue$physical_pool_id,saved_fixture$pool_screen$pools$pool_id))
  # Let the UI's selected theoretical default initialize before simulating a
  # user change to the empirical model; the app intentionally ignores init.
  session$flushReact()
  session$setInputs(sc_model_mode="equal",sc_k=10000,sc_target=95,sc_r=1,sc_max_k=100000,sc_iterations=2,sc_seed=2026,
    sc_basis="sites",sc_conditioning="physical",sc_universe="informative",sc_samples=catalogue$key,
    sc_allocation_mode="complexity",sc_technical_simulations=FALSE,
    sc_pmf_mode="planned",sc_depth=1,sc_diagnostic="length")
  session$flushReact()
  if(!grepl("Calculating coverage|Calculation complete",message()))stop(paste("Unexpected status after input initialization:",message()))
  stopifnot(length(output$sc_plot_primary)>0)
  # The default preview starts asynchronously. Cancel it here so the smoke
  # test can use the deterministic fixture for broad renderer assertions;
  # worker completion is exercised below with a focused one-pool request.
  stop_worker()
  result(saved_fixture);submitted(request());session$flushReact()
  stopifnot(request()$allocation_mode=="complexity",result()$pool_screen$settings$primary_strategy=="complexity")
  selected<-selected_allocations(result());stopifnot(all(selected$strategy=="complexity"),sum(selected$colonies)==selected_target_k(result()))
  predictions<-subset(result()$pool_screen$practical,k==10000)
  stopifnot(all(predictions$represented_denominator==4379),all(predictions$reference_denominator==5220))
  for(id in c("primary","allocation","allocation_curve","allocations_by_k","sampling_depth","union_reference","allele_curves","pool_gene_saturation","pool_allele_saturation","pool_site_saturation","marginal_curves","discovery","theory","redundancy","pmf","diagnostic")){
    base::message("Plot: ",id);stopifnot(length(output[[paste0("sc_plot_",id)]])>0)
  }
  for(id in c("primary_summary","primary_targets","primary_prediction","primary_allocation","primary_alleles","pool_inventory","pool_overlap","gene_access","parent_observations","pool_baselines","pool_probabilities","pool_abundance","pool_marginal","complexity_weights","pool_limitations")){
    base::message("Output: ",id);stopifnot(length(output[[paste0("sc_",id)]])>0)
  }
  before<-request();saved<-result()
  session$setInputs(sc_grid="show",sc_plot_style="presentation",sc_text_size="large")
  session$elapse(1000);session$flushReact()
  stopifnot(identical(before,request()),identical(saved,result()),is.null(job$process))
  session$setInputs(sc_allocation_mode="manual",sc_manual_1=1000,sc_manual_2=2000,sc_manual_3=3000,sc_manual_4=4000)
  session$flushReact();stopifnot(sum(unlist(request()$manual_allocation))==10000)
  # Worker integration with simulations off (the seeded 1000-iteration results
  # above already test every primary renderer). No styling can trigger sampling.
  # Exercise the real worker on a single pool; the all-pool deterministic
  # fixture above already covers rendering and allocation of all four pools.
  session$setInputs(sc_model_mode="equal",sc_samples=catalogue$key[catalogue$pool_label=="Pool 1"],
    sc_allocation_mode="manual",sc_manual_1=10000,sc_iterations=2,sc_run=1);session$flushReact()
  for(i in seq_len(600)){Sys.sleep(.1);session$elapse(500);session$flushReact();if(is.null(job$process))break}
  stopifnot(is.null(job$process),grepl("complete",message()),result()$pool_screen$settings$k==10000)
  x<-subset(result()$pool_screen$allocations,k==10000&strategy=="manual")
  stopifnot(nrow(x)==1,all(as.integer(x$colonies)==10000L))
})
cat("PASS: pool-aware UI, three-model plots, allocation, manual worker, diagnostics and cosmetic isolation\n")
