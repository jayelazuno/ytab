#!/usr/bin/env Rscript
# Quality Control: `ytab qc` builds the flagged library QC table for a real project, and the QC module shows it. The
# Overview table carries coloured cells with reasons and normalization hints; every themed subtab shows one plot at a
# time from its menu (each table metric plus the theme's curves), with the shared display options and downloads;
# Library Agreement compares libraries within the project or with other projects. TRANSIT is named only in the
# methodology references.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
repo_root<-normalizePath(args[[1L]],mustWork=TRUE)
suppressMessages({library(shiny);library(bslib)})
`%||%`<-function(x,y)if(is.null(x)||!length(x))y else x
for(f in c("process_helpers.R","project_discovery.R","project_state.R","ui_helpers.R","ui_components.R","plot_customization_helpers.R","plot_display_helpers.R",
           "job_progress.R","analysis_runs.R","navigation.R","qc_plot_utils.R","qc_mapping_stats_plot.R","qc_summary_library_plots.R",
           "qc_library_diagnostics_plots.R","ui_qc.R","qc_library_correlation.R","qc_module.R"))source(file.path(repo_root,"app/shiny/R",f))
html<-as.character(quality_control_ui())
values<-unique(sub('data-value="([a-z_]+)"',"\\1",regmatches(html,gregexpr('data-value="[a-z_]+"',html))[[1]]))
stopifnot(setequal(values,c("overview","reads_mapping","insertions","saturation_bias","agreement")),grepl("qc_table_run",html,fixed=TRUE),
          !grepl("transit|himar",gsub('(?s)<p class="ytab-reference">.*?</p>',"",html,perl=TRUE),ignore.case=TRUE),grepl('class="ytab-reference"',html,fixed=TRUE),grepl("qcp_insertions_plot_height",html,fixed=TRUE),grepl("qcp_agreement_heatmap_palette",html,fixed=TRUE),
          !grepl("mapping_qc_plot_choice|summary_qc_plot_choice|library_diagnostics_plot_choice|qc_library_plots",html),grepl("libcor_scope",html,fixed=TRUE))
for(old in names(navigation_qc_legacy))stopifnot(navigation_translate_legacy("quality_control",old)$nested_value==navigation_qc_legacy[[old]])
config<-file.path(repo_root,"output/projects/Zn_toxicity_screen/config/project.yaml")
project<-read_project_summary(config,repo_root)
temp<-file.path(tempfile(),"Zn");run<-file.path(temp,"qc","qc_smoke");dir.create(dirname(run),recursive=TRUE)
# The temporary project links to the real project's outputs; only its qc and library_correlation folders are new.
for(sub in c("mapfastq","create_hit_file","summary","library_diagnostics","config"))file.symlink(file.path(project$project_root,sub),file.path(temp,sub))
python<-unname(Sys.which("python"))
stopifnot(system2(python,c(file.path(repo_root,"scripts/local/ytab_qc.py"),"--project-config",config,"--output",run),stdout=FALSE,stderr=FALSE)==0L)
outputs<-c(readLines(file.path(run,"qc_flags.tsv")),readLines(file.path(run,"metrics.tsv")))
stopifnot(!any(grepl("transit|himar",outputs,ignore.case=TRUE)))
# A small correlation run with gene-level values, in the temporary project.
libcor<-file.path(temp,"library_correlation","libcor_smoke");dir.create(libcor,recursive=TRUE)
src<-"/private/tmp/claude-501/-Users-joshuaayelazuno-Desktop-ytab/bdf99455-d24f-41ee-a9a3-b73738c349e7/scratchpad/libcor_gale_zn"
if(dir.exists(src))invisible(file.copy(list.files(src,full.names=TRUE),libcor)) else
  stopifnot(system2(python,c(file.path(repo_root,"scripts/local/ytab_library_correlation.py"),"--project","Zn_toxicity_screen","--bootstrap","50","--output",libcor),stdout=FALSE,stderr=FALSE)==0L)
active<-reactive(modifyList(project,list(project_root=temp,project_id="Zn_toxicity_screen")))
started<-NULL;jobs<-list(job_is_running=function()FALSE,start_job=function(command,args,...)started<<-args)
testServer(function(input,output,session)qc<-qc_table_server(input,output,session,repo_root,active,function()config,jobs,function()python),{
  session$elapse(3100);session$flushReact()
  stopifnot(nrow(qc$runs())==1L,nrow(qc$table())==8L)
  d<-qc_display_table(qc$table(),qc$flags(),qc$metrics());stopifnot(nrow(d)==8L,ncol(d)==2L+nrow(qc$metrics()),!any(grepl("TRANSIT",unlist(d))))
  # Every QC-table metric is either in a menu or in the hover of the plot of the metric it rescales; menus have no repeats.
  listed<-sub("^metric:","",grep("^metric:",unlist(qc_plot_menus),value=TRUE));companions<-unlist(qc_metric_companions)
  stopifnot(setequal(union(listed,companions),qc$metrics()$metric),setequal(setdiff(companions,listed),c("insertions_per_kb","density","mean_ct","max_ct","mapped_reads")),!anyDuplicated(unlist(qc_plot_menus)),
            !any(c("role_summary","jackpot")%in%unlist(qc_plot_menus)))
  shown<-0L;empty<-character()
  for(id in c("reads_mapping","insertions","saturation_bias","agreement")){
    choices<-unname(unlist(qc_plot_menus[[id]]))
    for(choice in choices){
      if(choice=="pair_scatter")session$setInputs(qcp_agreement_lib_a="Gale_2020::SRR11561506",qcp_agreement_lib_b="Gale_2020::SRR11561507")
      do.call(session$setInputs,stats::setNames(list(choice),paste0("qcp_",id,"_plot")))
      out<-tryCatch(output[[paste0("qcp_",id,"_plotly")]],error=function(e)if(inherits(e,"shiny.silent.error")||inherits(e,"validation")||grepl("^No ",conditionMessage(e)))paste("EMPTY",conditionMessage(e))else paste("ERROR",conditionMessage(e)))
      if(startsWith(out,"EMPTY")){empty<-c(empty,sprintf("%s/%s: %s",id,choice,sub("^EMPTY ","",out)));next}
      if(!grepl("plotly",out,fixed=TRUE))stop(sprintf("%s / %s did not draw: %s",id,choice,substr(out,1,200)))
      shown<-shown+1L
    }
  }
  # Display options change the plot; downloads work for the plot shown.
  session$setInputs(qcp_insertions_plot="metric:skewness",qcp_insertions_bar_orientation="horizontal",qcp_insertions_show_value_labels=TRUE,qcp_insertions_text_size="large")
  b<-qc$built$insertions();stopifnot(inherits(b$plot$coordinates,"CoordFlip"),any(vapply(b$plot$layers,function(l)inherits(l$geom,"GeomText"),TRUE)))
  for(id in c("qcp_insertions_download_png","qcp_insertions_download_pdf","qcp_insertions_download_data","qc_table_download_xlsx","qc_table_download_tsv","qc_table_download_flags"))stopifnot(file.size(output[[id]])>100)
  session$setInputs(qcp_insertions_plot="metric:unique_sites");b<-qc$built$insertions()
  stopifnot(any(grepl("Insertions per kb of genome",b$plot$data$tip,fixed=TRUE)),any(grepl("Catalogue density",b$plot$data$tip,fixed=TRUE)))
  session$setInputs(qcp_insertions_plot="gene_hit_bins");stopifnot(identical(levels(qc$built$insertions()$plot$data$bin),c("0","1\u20133","4\u20137","8\u20139","10\u201319","\u226520")))
  session$setInputs(qcp_saturation_bias_plot="midlc",qcp_saturation_bias_colour="role");stopifnot(grepl(".role",rlang::quo_text(qc$built$saturation_bias()$plot$mapping$colour),fixed=TRUE))
  # Role filter on the library-correlation heatmaps: only the ticked roles are drawn.
  r<-qc$libcor_results();roles<-libcor_roles(rownames(r$matrix),dirname(temp),r$samples);one<-names(which.max(table(roles)))
  session$setInputs(qcp_agreement_plot="cor_heatmap",qcp_agreement_roles=one)
  stopifnot(length(unique(qc$built$agreement()$data$a))==sum(roles==one),sum(roles==one)<length(roles))
  session$setInputs(qc_table_run=1);stopifnot(any(grepl("ytab_qc.py",started)))
  cat(sprintf("plots drawn one at a time: %d; shown as an empty-state message: %d\n",shown,length(empty)));if(length(empty))cat(" ",empty,sep="\n  ")
  stopifnot(all(grepl("input_reads|retained_fraction|adapter_fraction",empty)))  # Zn was mapped elsewhere: no cleanup records
})
cat("PASS: QC table and one-plot-at-a-time QC subtabs (metrics, curves, agreement), display options, downloads; TRANSIT named only in references\n")
