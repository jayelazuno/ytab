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
           "qc_library_diagnostics_plots.R","ui_qc.R","reference_projects.R","qc_library_correlation.R","qc_module.R"))source(file.path(repo_root,"app/shiny/R",f))
html<-as.character(quality_control_ui())
values<-unique(sub('data-value="([a-z_]+)"',"\\1",regmatches(html,gregexpr('data-value="[a-z_]+"',html))[[1]]))
stopifnot(setequal(values,c("overview","reads_mapping","insertions","saturation_bias","agreement")),grepl("qc_table_run",html,fixed=TRUE),
          !grepl("transit|himar",gsub('(?s)<p class="ytab-reference">.*?</p>',"",html,perl=TRUE),ignore.case=TRUE),grepl('class="ytab-reference"',html,fixed=TRUE),grepl("qcp_insertions_plot_height",html,fixed=TRUE),grepl("qcp_agreement_heatmap_palette",html,fixed=TRUE),
          !grepl("mapping_qc_plot_choice|summary_qc_plot_choice|library_diagnostics_plot_choice|qc_library_plots",html),grepl("libcor_scope",html,fixed=TRUE),grepl("qcp_reads_mapping_reference_selector",html,fixed=TRUE))
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
reference<-ytab_reference_state()
others<-discover_ytab_projects(repo_root,include_tests=FALSE)
testServer(function(input,output,session)qc<-qc_table_server(input,output,session,repo_root,active,function()config,jobs,function()python,reactive(others),reference),{
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
  b<-qc$built$insertions();stopifnot(rlang::quo_text(b$plot$mapping$y)==".key",any(vapply(b$plot$layers,function(l)inherits(l$geom,"GeomText"),TRUE)))
  for(id in c("qcp_insertions_download_png","qcp_insertions_download_pdf","qcp_insertions_download_data","qc_table_download_xlsx","qc_table_download_tsv","qc_table_download_flags"))stopifnot(file.size(output[[id]])>100)
  session$setInputs(qcp_insertions_plot="metric:unique_sites");b<-qc$built$insertions()
  stopifnot(any(grepl("Insertions per kb of genome",b$plot$data$tip,fixed=TRUE)),any(grepl("Catalogue density",b$plot$data$tip,fixed=TRUE)))
  session$setInputs(qcp_insertions_plot="gene_hit_bins");stopifnot(identical(levels(qc$built$insertions()$plot$data$bin),c("0","1\u20133","4\u20137","8\u20139","10\u201319","\u226520")))
  session$setInputs(qcp_saturation_bias_plot="midlc",qcp_saturation_bias_colour="role");stopifnot(grepl(".role",rlang::quo_text(qc$built$saturation_bias()$plot$mapping$colour),fixed=TRUE))
  # Role filter on the library-correlation heatmaps: only the ticked roles are drawn.
  r<-qc$libcor_results();roles<-libcor_roles(rownames(r$matrix),dirname(temp),r$samples);one<-names(which.max(table(roles)))
  session$setInputs(qcp_agreement_plot="cor_heatmap",qcp_agreement_roles=one)
  stopifnot(length(unique(qc$built$agreement()$data$a))==sum(roles==one),sum(roles==one)<length(roles))
  # Reference projects: opt-in, so nothing is drawn until one is chosen; then the plot gains its band and median line
  # and the plotted-data download carries the reference libraries themselves.
  session$setInputs(qcp_reads_mapping_plot="metric:percent_mapped")
  stopifnot(is.null(qc$built$reads_mapping()$reference))
  gale<-others$project_config[others$project_id=="Gale_2020"]
  if(length(gale)){
    reference(gale[[1L]]);session$flushReact()
    b<-qc$built$reads_mapping();s0<-b$reference$summary
    # Bars: the reference libraries join the plot as bars of their own, filled by project rather than by QC flag.
    own<-nrow(qc$table());stopifnot(nrow(s0)==1L,s0$libraries>0L,is.finite(s0$median),
              nrow(b$plot$data)==own+s0$libraries,s0$.project%in%levels(b$plot$data$.fill),
              grepl("Gale_2020",ytab_reference_caption(s0,project$species),fixed=TRUE))
    # Library names repeat between projects, so each bar is keyed by project as well as by name.
    stopifnot(!anyDuplicated(as.character(b$plot$data$.key)))
    # Panels and lines are the other two ways of drawing the same reference.
    # Panels: facet_grid with free space, so a panel's height follows the number of libraries it holds.
    session$setInputs(qcp_reads_mapping_reference_mode="panels")
    pf<-qc$built$reads_mapping()$plot$facet
    stopifnot(inherits(pf,"FacetGrid"),isTRUE(pf$params$space_free$y)||isTRUE(pf$params$space_free$x))
    session$setInputs(qcp_reads_mapping_reference_mode="lines")
    stopifnot(sum(vapply(qc$built$reads_mapping()$plot$layers,function(l)class(l$geom)[[1L]],"")=="GeomHline")==2L)
    session$setInputs(qcp_reads_mapping_reference_mode="auto")
    # Automatic switches to panels once the bars no longer fit on one axis.
    stopifnot(identical(ytab_reference_mode("auto",8L,5L),"bars"),identical(ytab_reference_mode("auto",8L,60L),"panels"))
    d<-ytab_reference_download(b$data,b$reference,"Zn_toxicity_screen")
    stopifnot(length(unique(d$project))==2L,nrow(d)>nrow(b$data))
    # Every themed subtab draws the same reference: a metric as bars, a distribution and a curve as further series.
    session$setInputs(qcp_insertions_plot="metric:skewness")
    stopifnot(nrow(qc$built$insertions()$reference$summary)==1L)
    session$setInputs(qcp_insertions_plot="dist_histogram")
    hist<-qc$built$insertions()$data
    stopifnot(any(grepl("^Gale_2020 \u00b7 ",unique(hist$.library))),any(!grepl("^Gale_2020 \u00b7 ",unique(hist$.library))))
    session$setInputs(qcp_saturation_bias_plot="metric:depth_ratio")
    stopifnot(nrow(qc$built$saturation_bias()$reference$summary)==1L)
    # Library Agreement compares projects through its own correlation run, so it is left out.
    session$setInputs(qcp_agreement_plot="cor_heatmap");stopifnot(is.null(qc$built$agreement()$reference))
    # Every plot still draws with a reference chosen, and the ones that can show it do.
    drawn<-0L;overlaid<-character();quiet<-character()
    for(id in c("reads_mapping","insertions","saturation_bias")){
      for(choice in unname(unlist(qc_plot_menus[[id]]))){
        do.call(session$setInputs,stats::setNames(list(choice),paste0("qcp_",id,"_plot")))
        out<-tryCatch(output[[paste0("qcp_",id,"_plotly")]],error=function(e)
          if(inherits(e,"shiny.silent.error")||inherits(e,"validation")||grepl("^No ",conditionMessage(e)))"EMPTY" else paste("ERROR",conditionMessage(e)))
        if(startsWith(out,"ERROR"))stop(sprintf("%s / %s broke with a reference project: %s",id,choice,substr(out,1,200)))
        if(identical(out,"EMPTY"))next
        drawn<-drawn+1L
        b3<-qc$built[[id]]()
        text<-unlist(lapply(c(".library",".project",".group","display_group"),function(k)as.character(b3$data[[k]] %||% character())))
        shows<-nrow(b3$reference$summary %||% data.frame())>0L||any(grepl("Gale_2020",text,fixed=TRUE))
        if(shows)overlaid<-c(overlaid,choice) else quiet<-c(quiet,paste0(id,"/",choice))
      }
    }
    cat(sprintf("with a reference chosen: %d plots drew, %d of them show it\n",drawn,length(overlaid)))
    if(length(quiet))cat("  no reference overlay:",paste(quiet,collapse=", "),"\n")
    stopifnot(drawn>0L,length(overlaid)>=8L)
    # Unticking removes the project: its bars leave the plot and the plot is the project's own again.
    reference(character());session$flushReact()
    session$setInputs(qcp_reads_mapping_plot="metric:percent_mapped")
    b2<-qc$built$reads_mapping()
    stopifnot(nrow(b2$reference$summary %||% data.frame())==0L,nrow(b2$plot$data)==own,
              setequal(levels(b2$plot$data$.fill),names(YTAB_FLAG_COLOURS)))
  }
  session$setInputs(qc_table_run=1);stopifnot(any(grepl("ytab_qc.py",started)))
  cat(sprintf("plots drawn one at a time: %d; shown as an empty-state message: %d\n",shown,length(empty)));if(length(empty))cat(" ",empty,sep="\n  ")
  stopifnot(all(grepl("input_reads|retained_fraction|adapter_fraction",empty)))  # Zn was mapped elsewhere: no cleanup records
})
cat("PASS: QC table and one-plot-at-a-time QC subtabs (metrics, curves, agreement), display options, downloads, opt-in reference projects; TRANSIT named only in references\n")
