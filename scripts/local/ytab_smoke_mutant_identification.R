#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2L)
repo_root<-normalizePath(args[1]);fixture<-normalizePath(args[2])
library(shiny);library(bslib)
`%||%`<-function(x,y)if(is.null(x)||!length(x))y else x
for(name in c("ui_helpers.R","ui_components.R","table_display_helpers.R","project_discovery.R",
              "process_helpers.R","fitness_result_state.R","glabrata_annotation_lookup.R","qc_plot_utils.R","mutant_screen.R"))source(file.path(repo_root,"app/shiny/R",name))
html<-as.character(mutant_screen_ui())
stopifnot(grepl("gel_image.png",html,fixed=TRUE),grepl("HR-2",html,fixed=TRUE),grepl("2005.02.011",html,fixed=TRUE))
stopifnot(!grepl("Joshua",html,fixed=TRUE),identical(ms_natural_order(c("10-HR2","2-HR2","1-HR2")),c(3L,2L,1L)))
stopifnot(!grepl("Gene names use S. cerevisiae common names",html,fixed=TRUE))
for(id in c("ms_inspect_mutant","ms_gene_information")){
  opening<-regmatches(html,regexpr(paste0('<details[^>]*id="',id,'"[^>]*>'),html))
  stopifnot(nzchar(opening),!grepl("\\bopen(=|[ >])",opening,perl=TRUE))
}
project<-read_project_summary(file.path(repo_root,"output/projects/H2O2_screen_v1/config/project.yaml"),repo_root)
stopifnot(identical(ms_reference_label(project),"glabrata — GCA_014217725.1"),
  nchar(ms_reference_label(list(species="other",reference=list(fasta=paste(rep("x",100),collapse="")))))<40L)
fixture_data<-jsonlite::fromJSON(fixture,simplifyVector=FALSE)
lookup<-ytab_load_glabrata_annotation_lookup(repo_root)
fallback<-lookup[nzchar(lookup$gwk60_id_clean)&nzchar(lookup$cagl_id)&!nzchar(lookup$scer_gene_name),,drop=FALSE]
stopifnot(nrow(fallback)>0L)
fallback_display<-ms_gene_display(fallback$gwk60_id_clean[[1]],repo_root)
stopifnot(identical(fallback_display[["Gene name"]],fallback_display[["Gene ID"]]))
display<-ms_identification_table(fixture_data$results,repo_root,project)
first_read<-Filter(function(r)r$mutant=="1-HR2",fixture_data$results)[[1]]
pre_qc_example<-first_read;pre_qc_example$hermes_overlap_before_qc_bases<-12;pre_qc_example$hermes_overlap_bases<-0
pre_qc_display<-ms_identification_table(list(pre_qc_example),repo_root,project)
stopifnot(pre_qc_display[["Hermes overlap (bp)"]][1]==12,
  pre_qc_display[["Aligned bases after QC (bp)"]][1]==first_read$aligned_bases)
visible<-function(x)names(x)[!startsWith(names(x),".ytab_")]
stopifnot(identical(visible(display),c("Mutant","Gene ID","Gene name","Chromosome","Aligned bases after QC (bp)","Hermes overlap (bp)","Identity (%)")),
  nrow(display)==length(fixture_data$results),identical(display$Mutant,paste0(1:10,"-HR2")),
  display[["Aligned bases after QC (bp)"]][1]==first_read$aligned_bases,
  is.na(display[["Hermes overlap (bp)"]][1]),
  display[["Identity (%)"]][1]==round(100*first_read$identity,2),
  grepl(";",display[["Gene ID"]][display$Mutant=="10-HR2"],fixed=TRUE),
  "CAGL0H05885g"%in%display[["Gene ID"]],"FMP30"%in%display[["Gene name"]],
  !any(grepl("^CP",display$Chromosome)),"Chr H"%in%display$Chromosome)
candidate<-Filter(function(r)any(vapply(r$candidates,function(g)g$gene_id=="GWK60_H05819",FALSE)),fixture_data$results)[[1]]
stopifnot(length(visible(ms_candidate_table(candidate,repo_root,project)))==5L,
  length(visible(ms_linked_table(candidate,repo_root)))==5L,
  length(visible(ms_linked_table(candidate,repo_root,"Raw sample summary")))==5L)
widget<-ms_identification_widget(display)
stopifnot(inherits(ytab_gene_details_datatable(data.frame()),"datatables"))
stopifnot(grepl("ytab-gene-detail-toggle",paste(widget$x$options$columnDefs,collapse=" "),fixed=TRUE),
  identical(widget$x$data[[".ms_mutant_order"]],seq_len(10)),
  any(vapply(widget$x$options$columnDefs,function(d)identical(d$orderData,10L)&&identical(d$targets,0),FALSE)))
intergenic<-first_read;intergenic$candidates<-list();intergenic$gene_assignment<-list(rule="no_cds_overlap")
stopifnot(identical(ms_identification_table(list(intergenic),repo_root,project)[["Gene name"]],"Intergenic"))
test_root<-tempfile("ytab-mutant-ui-");dir.create(test_root);project$project_root<-test_root
storage_root<-file.path(test_root,"universal_mutant_screen")
active<-function()project;active_project_path<-function()project$project_config
python_bin<-function()normalizePath(file.path(R.home(),"../..","bin","python"))
go_to<-function(...)invisible(NULL)
test_server<-mutant_screen_server;formals(test_server)<-alist(input=,output=,session=)
withCallingHandlers(shiny::testServer(test_server,{
  session$flushReact();results(jsonlite::fromJSON(fixture,simplifyVector=FALSE))
  r<-Filter(function(x)length(x$candidates)>0,reads())[[1]]
  session$setInputs(ms_selected=r$read_id)
  stopifnot(length(output$ms_results)>0,length(output$ms_genes)>0,length(output$ms_screen_results)>0,
    length(output$ms_details)>0,length(output$ms_explorer_control)>0,is.null(job$process))
  session$setInputs(ms_input_mode="path",ms_paths="data/mutant_screen/pfasta-Oct13-9-18-48.tar",ms_run=1,
    ms_min_length=80,ms_identity=90,ms_mapq=20,ms_quality=20,ms_overlap=12,ms_error=.15)
  session$flushReact();stopifnot(!is.null(job$process))
  for(i in seq_len(600)){Sys.sleep(.05);session$elapse(500);session$flushReact();if(is.null(job$process))break}
  stopifnot(is.null(job$process),identical(status(),"Identification complete."),length(reads())==10L,
    file.exists(file.path(results()$output_dir,"mutant_identifications.csv")),!session$isClosed())
  saved<-results();stopifnot(all(vapply(reads(),function(r)!is.null(r$hermes_overlap_bases)&&!is.null(r$hermes_overlap_before_qc_bases)&&nzchar(r$junction_review_reason),FALSE)))
  stopifnot(identical(saved$project_id,"universal"),is.null(saved$linked_project_id),
    all(vapply(reads(),function(r)!length(r$screen_results),FALSE)),
    identical(saved$output_dir,normalizePath(file.path(storage_root,"current"))))
  previous_dirs<-list.dirs(storage_root,recursive=FALSE)
  session$setInputs(ms_run=2);session$flushReact()
  for(i in seq_len(600)){Sys.sleep(.05);session$elapse(500);session$flushReact();if(is.null(job$process))break}
  stopifnot(identical(status(),"Previous matching identification reused."),identical(saved,results()),
    identical(previous_dirs,list.dirs(storage_root,recursive=FALSE)))
  session$setInputs(ms_saved_run=file.path(saved$output_dir,"results.json"))
  stopifnot(identical(status(),"Saved identification run loaded."),identical(saved,results()))
  bad<-file.path(storage_root,"archive/bad_run");dir.create(bad,recursive=TRUE)
  writeLines("not JSON",file.path(bad,"results.json"));session$setInputs(ms_saved_run=file.path(bad,"results.json"))
  stopifnot(grepl("Could not load",status()),identical(saved,results()),!session$isClosed())
  project$project_id<-"another_open_project";session$flushReact()
  stopifnot(identical(saved,results()),!grepl("another_open_project",output$ms_link_context,fixed=TRUE))
  session$setInputs(ms_input_mode="upload",
    ms_files=data.frame(name="uploaded.seq",datapath=file.path(repo_root,"data/mutant_screen/pfasta-Oct13-9-18-48.tar"),size=1,type="text/plain"),ms_run=3)
  session$flushReact();stopifnot(!is.null(job$process));session$setInputs(ms_cancel=1)
  stopifnot(is.null(job$process),identical(saved,results()))
}),warning=function(e)stop(conditionMessage(e),call.=FALSE))
cat("PASS: method/gel, tables, real async archive run, review flags, saved-run loading, upload launch and cancellation\n")
