relative_project_path <- function(path, project_root) {
  path <- normalizePath(path, winslash="/", mustWork=FALSE); root <- normalizePath(project_root,winslash="/",mustWork=FALSE)
  if (startsWith(path,paste0(root,"/"))) substring(path,nchar(root)+2L) else basename(path)
}

read_qc_csvs <- function(root, pattern) {
  files <- if(dir.exists(root)) list.files(root,pattern=pattern,recursive=TRUE,full.names=TRUE) else character()
  rows <- lapply(files,function(file)tryCatch(transform(read.csv(file,stringsAsFactors=FALSE,check.names=FALSE),.source=file),error=function(e)NULL))
  rows <- Filter(Negate(is.null),rows); if(length(rows))do.call(rbind,rows)else data.frame()
}

format_qc_metrics <- function(data, count_columns=character(), percentage_columns=character(), mean_columns=character()) {
  for(name in intersect(count_columns,names(data))) data[[name]] <- ifelse(is.na(data[[name]]),"",format(round(as.numeric(data[[name]])),big.mark=",",scientific=FALSE,trim=TRUE))
  for(name in intersect(percentage_columns,names(data))) { value<-suppressWarnings(as.numeric(sub("%$","",as.character(data[[name]])))); data[[name]]<-ifelse(is.na(value),"",sprintf("%.1f%%",value)) }
  for(name in intersect(mean_columns,names(data))) data[[name]] <- ifelse(is.na(data[[name]]),"",sprintf("%.2f",as.numeric(data[[name]])))
  data
}

mapping_qc_data <- function(project) {
  raw<-read_qc_csvs(file.path(project$project_root,"mapfastq"),"mapping_stats.*\\.csv$"); if(!nrow(raw))return(raw)
  sample<-raw$sample %||% basename(dirname(raw$.source)); bam<-file.path(project$project_root,"mapfastq",sample,paste0(sample,".sorted.bam"))
  result<-data.frame(Sample=sample,check.names=FALSE,stringsAsFactors=FALSE)
  available<-list(`Reads processed`="total_records",`Aligned reads`="primary_mapped",`Alignment percentage`="percent_mapped")
  for(label in names(available))if(available[[label]]%in%names(raw))result[[label]]<-raw[[available[[label]]]]
  result$BAM<-ifelse(file.exists(bam),"Available","Missing");result$`Mapping statistics`<-"Available"
  attr(result,"details")<-data.frame(Sample=sample,BAM=vapply(bam,relative_project_path,"",project_root=project$project_root),Statistics=vapply(raw$.source,relative_project_path,"",project_root=project$project_root),stringsAsFactors=FALSE)
  format_qc_metrics(result,c("Reads processed","Aligned reads"),"Alignment percentage")
}

summary_qc_data <- function(project) {
  raw<-read_qc_csvs(file.path(project$project_root,"summary"),"summary_stats.*\\.csv$");if(!nrow(raw))return(raw)
  wanted<-c("sample","total_reads","total_hits","percent_hits_in_features","percent_intergenic_hits","percent_features_hit","mean_hits_per_feature"); visible<-raw[intersect(wanted,names(raw))]
  names(visible)<-c(sample="Sample",total_reads="Total reads",total_hits="Total hits",percent_hits_in_features="Hits in features",percent_intergenic_hits="Intergenic hits",percent_features_hit="Features hit",mean_hits_per_feature="Mean hits / feature")[names(visible)]
  detail_names<-setdiff(names(raw),c(wanted,".source"));details<-raw[,c("sample",detail_names),drop=FALSE];details$source<-vapply(raw$.source,relative_project_path,"",project_root=project$project_root);attr(visible,"details")<-details
  format_qc_metrics(visible,c("Total reads","Total hits"),c("Hits in features","Intergenic hits","Features hit"),"Mean hits / feature")
}

summary_qc_cards <- function(project) {
  raw<-read_qc_csvs(file.path(project$project_root,"summary"),"summary_stats.*\\.csv$"); num<-function(name){if(!name%in%names(raw))numeric()else suppressWarnings(as.numeric(sub("%$","",raw[[name]])))}; safe<-function(x,fun,digits=0){x<-x[is.finite(x)];if(!length(x))"Unavailable"else format(round(fun(x),digits),big.mark=",",trim=TRUE)}
  tags$div(class="ytab-stat-grid",tags$div(tags$b(nrow(raw)),"Samples"),tags$div(tags$b(safe(num("total_hits"),median)),"Median total hits"),tags$div(tags$b(paste0(safe(num("percent_features_hit"),median,1),"%")),"Median features hit"),tags$div(tags$b(paste0(safe(num("percent_features_hit"),min,1),"%")),"Minimum features hit"),tags$div(tags$b(paste0(safe(num("percent_features_hit"),max,1),"%")),"Maximum features hit"))
}

compact_qc_table <- function(data) DT::datatable(data,rownames=FALSE,selection="none",options=list(pageLength=10,lengthMenu=c(10,25,50),autoWidth=FALSE,ordering=TRUE,searching=TRUE,scrollX=TRUE,columnDefs=list(list(className="dt-right",targets=which(vapply(data,is.numeric,FALSE))-1L),list(className="ytab-nowrap",targets=0))))

qc_sample_eligibility <- function(project,status) {
  samples<-project$samples; included<-if("include"%in%names(samples))tolower(as.character(samples$include))%in%c("true","1","yes")else rep(TRUE,nrow(samples));idx<-match(samples$sample,status$sample);hit_status<-status$hit_file[idx];path<-status$hits_path[idx];valid<-included & hit_status=="success" & vapply(path,nonempty_file,FALSE)
  data.frame(Sample=as.character(samples$sample),Condition=as.character(samples$guessed_condition%||%""),Background=as.character(samples$guessed_background%||%""),Pool=as.character(samples$guessed_pool%||%""),`Hit-file status`=ifelse(valid,"Available",ifelse(included,"Unavailable","Excluded")),Eligible=ifelse(valid,"Yes","No"),Reason=ifelse(valid,"",ifelse(!included,"Not included in project","Valid successful hit file unavailable")),hit_path=path,stringsAsFactors=FALSE,check.names=FALSE)
}

diagnostic_cache_manifest <- function(project,run_id) { path<-file.path(project$project_root,"manifests","library_diagnostics","runs",run_id,"manifest.json");if(!file.exists(path))return(NULL);tryCatch(jsonlite::fromJSON(path,simplifyVector=FALSE),error=function(e)NULL) }
diagnostic_action_label <- function(mode="preview",cache_match=FALSE,force=FALSE) if(mode=="preview")"Preview diagnostics"else if(force)"Rerun diagnostics"else if(cache_match)"Use cached diagnostics"else"Run diagnostics"
diagnostic_result_heading <- function(status) switch(status,preview="Preview complete",cached="Cached result reused",success="Library Diagnostics complete",failed="Library Diagnostics failed","Library Diagnostics")

human_file_size <- function(bytes) {bytes<-as.numeric(bytes);vapply(bytes,function(x)if(is.na(x))"Unknown"else if(x<1024)sprintf("%d B",round(x))else if(x<1024^2)sprintf("%.0f KB",x/1024)else if(x<1024^3)sprintf("%.1f MB",x/1024^2)else sprintf("%.1f GB",x/1024^3),"")}
diagnostic_plot_type <- function(filename) {value<-tolower(filename);if(grepl("centromere_bias",value))"Centromere bias"else if(grepl("metaplots",value))"Feature metaplots"else if(grepl("tss_metaplot",value))"TSS metaplot"else if(grepl("trna_metaplot",value))"tRNA metaplot"else if(grepl("seqbias",value))"Sequence bias"else if(grepl("midlc",value))"MidLC"else"Diagnostic plot"}
diagnostic_plot_title <- function(filename) diagnostic_plot_type(filename)

diagnostic_resource_roots <- function(project) {
  roots <- c(
    `ytab-diagnostics-project` = file.path(project$project_root, "library_diagnostics"),
    `ytab-diagnostics-export` = file.path(project$export_root, "qc", "diagnostics")
  )
  roots[dir.exists(roots)]
}

diagnostic_path_within <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  startsWith(path, paste0(root, "/"))
}

diagnostic_encode_relative_path <- function(path) {
  parts <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1L]]
  if (!length(parts) || any(parts %in% c("", ".", ".."))) return("")
  paste(vapply(parts, URLencode, character(1), reserved = FALSE), collapse = "/")
}

diagnostic_served_file <- function(path, project) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  display <- relative_project_path(path, project$project_root)
  result <- list(served_url = "", preview_available = FALSE,
                 relative_path = display,
                 preview_message = "Preview unavailable; file exists but could not be served.")
  if (!file.exists(path) || dir.exists(path)) {
    result$preview_message <- "Preview unavailable; file is missing."
    return(result)
  }
  roots <- diagnostic_resource_roots(project)
  for (route in names(roots)) {
    root <- roots[[route]]
    if (!diagnostic_path_within(path, root)) next
    relative <- substring(path, nchar(normalizePath(root, winslash = "/", mustWork = FALSE)) + 2L)
    encoded <- diagnostic_encode_relative_path(relative)
    if (!nzchar(encoded)) return(result)
    result$served_url <- paste0(route, "/", encoded)
    result$preview_available <- TRUE
    result$relative_path <- if (identical(route, "ytab-diagnostics-project"))
      file.path("library_diagnostics", relative) else
      file.path("exports", "qc", "diagnostics", relative)
    return(result)
  }
  result
}

build_diagnostic_file_inventory <- function(project) {
  roots<-c(file.path(project$project_root,"library_diagnostics"),file.path(project$export_root,"qc","diagnostics"),file.path(project$project_root,"manifests","library_diagnostics"),file.path(project$project_root,"logs","library_diagnostics"));files<-unique(unlist(lapply(roots,function(root)if(dir.exists(root))list.files(root,recursive=TRUE,full.names=TRUE)else character())));files<-files[file.exists(files)&!dir.exists(files)];if(!length(files))return(data.frame())
  served<-lapply(files,diagnostic_served_file,project=project);filename<-basename(files);ext<-tolower(tools::file_ext(files));is_plot<-ext%in%c("png","jpg","jpeg","svg");type<-ifelse(is_plot,"Plot",ifelse(ext%in%c("csv","tsv"),"Table",ifelse(ext=="json","Manifest",ifelse(ext=="html","HTML",ifelse(ext%in%c("txt","log"),"Text","Log")))));relative<-vapply(served,`[[`,"","relative_path");parts<-strsplit(relative,"/");sample<-vapply(parts,function(x)if(length(x)>2&&x[[1]]%in%c("library_diagnostics","exports")&&!x[[length(x)-1]]%in%c("runs","library_diagnostics","diagnostics"))x[[length(x)-1]]else"Project","");run<-vapply(parts,function(x){index<-match("runs",x);if(!is.na(index)&&length(x)>index)x[[index+1L]]else"current_legacy"},"");info<-file.info(files)
  data.frame(file=files,filename=filename,extension=ext,display_type=type,plot_type=vapply(filename,diagnostic_plot_type,""),display_title=vapply(filename,diagnostic_plot_title,""),sample=sample,sample_set=ifelse(run=="current_legacy","Legacy diagnostics result",tools::toTitleCase(gsub("_"," ",run))),run_id=run,size_bytes=as.numeric(info$size),size_display=human_file_size(info$size),modified=format(info$mtime,"%Y-%m-%d %H:%M"),relative_path=relative,served_url=vapply(served,`[[`,"","served_url"),preview_available=vapply(served,`[[`,FALSE,"preview_available"),preview_message=vapply(served,`[[`,"","preview_message"),viewable=ext%in%c("png","jpg","jpeg","svg","csv","tsv","txt","html"),is_plot=is_plot,stringsAsFactors=FALSE)
}

filter_diagnostic_inventory <- function(data,sample="All",plot_type="All",sample_set="All",run_id="All",type="All") {
  if(!nrow(data))return(data);keep<-rep(TRUE,nrow(data));if(sample!="All")keep<-keep&data$sample==sample;if(plot_type!="All")keep<-keep&data$plot_type==plot_type;if(sample_set!="All")keep<-keep&data$sample_set==sample_set;if(run_id!="All")keep<-keep&data$run_id==run_id;if(type!="All")keep<-keep&data$display_type==type;data[keep,,drop=FALSE]
}
paginate_diagnostic_inventory <- function(data,page=1L,page_size=8L) {pages<-max(1L,ceiling(nrow(data)/page_size));page<-max(1L,min(as.integer(page),pages));start<-(page-1L)*page_size+1L;rows<-if(nrow(data)&&start<=nrow(data))seq.int(start,min(start+page_size-1L,nrow(data)))else integer();list(data=data[rows,,drop=FALSE],page=page,pages=pages,total=nrow(data))}

qc_mapping_ui <- function() panel_card(
  "Mapping QC",
  ytab_two_column_layout(
    controls = ytab_control_panel(
      "Mapping display",
      uiOutput("mapping_qc_readiness"),
      selectInput("mapping_qc_plot_choice", "Visualization",
                  choices = c("Total and mapped reads" = "read_counts",
                              "Percent mapped reads" = "percent_mapped",
                              "Mapping quality score" = "mapq"),
                  selected = "read_counts"),
      ytab_plot_customization_controls("mapping_qc", include_bars = TRUE,
                                       default_height = "medium",
                                       default_width = "wide",
                                       default_bar_orientation = "horizontal",
                                       default_show_value_labels = FALSE),
      tags$details(
        class = "ytab-more-options",
        tags$summary("Tables / downloads"),
        tags$div(class = "ytab-actions",
                 downloadButton("download_mapping_qc_plot", "Download plot"),
                 downloadButton("download_mapping_qc_plotted_data", "Download plotted data"),
                 downloadButton("download_mapping_qc_table", "Download mapping summary"),
                 downloadButton("download_mapping_qc_details", "Download file details")),
        DT::DTOutput("mapping_qc_table"),
        tags$details(tags$summary("File details"), DT::DTOutput("mapping_qc_details"))
      )
    ),
    main = ytab_plot_card(
      "Mapping summary",
      tagList(
        uiOutput("mapping_qc_selected_plot"),
        tags$details(class = "ytab-more-options",
                     tags$summary("Plot provenance"),
                     tags$p(class = "text-muted",
                            "Rendered from existing mapping statistics tables."))
      )
    )
  )
)

qc_summary_ui <- function() panel_card(
  "Summary QC",
  ytab_two_column_layout(
    controls = ytab_control_panel(
      "Summary display",
      uiOutput("summary_qc_readiness"),
      selectInput(
        "summary_qc_plot_choice", "Visualization",
        choices = c("Library complexity" = "complexity",
                    "Features hit" = "features",
                    "Combined features hit" = "combined_features",
                    "Reads per insertion site" = "reads_per_hit",
                    "Feature vs intergenic hits" = "feature_intergenic",
                    "Genome-wide binned signal" = "genome_bins",
                    "Library concordance" = "pairwise")
      ),
      uiOutput("summary_qc_combined_group_selector"),
      ytab_plot_customization_controls("summary_qc", include_bars = TRUE,
                                       default_height = "medium",
                                       default_bar_orientation = "horizontal",
                                       default_show_value_labels = FALSE),
      tags$details(
        class = "ytab-more-options",
        tags$summary("Tables / downloads"),
        uiOutput("summary_qc_cards"),
        tags$div(class = "ytab-actions",
                 downloadButton("download_summary_qc_plot", "Download plot"),
                 downloadButton("download_summary_qc_plotted_data", "Download plotted data"),
                 downloadButton("download_summary_qc_table", "Download summary table"),
                 downloadButton("download_summary_qc_details", "Download detailed metrics")),
        DT::DTOutput("summary_qc_table"),
        tags$details(tags$summary("Detailed metrics"), DT::DTOutput("summary_qc_details"))
      )
    ),
    main = ytab_plot_card(
      "Summary QC",
      tagList(
        uiOutput("summary_qc_selected_plot"),
        tags$details(class = "ytab-more-options",
                     tags$summary("Plot provenance"),
                     tags$p(class = "text-muted",
                            "Rendered from existing SummaryTable outputs."))
      )
    )
  )
)

qc_library_diagnostics_ui <- function() panel_card(
  "Library Diagnostics",
  tags$section(
    class = "ytab-qc-section",
    tags$h4("1. Select samples"),
    uiOutput("library_diagnostics_readiness"),
    uiOutput("qc_current_selection"),
    selectizeInput("qc_samples", "Samples to diagnose", choices = character(),
                   multiple = TRUE, options = list(placeholder = "Select samples…")),
    uiOutput("qc_sample_presets"),
    uiOutput("qc_selection_message"),
    tags$details(tags$summary("Review eligible samples"), DT::DTOutput("qc_sample_preview"))
  ),
  tags$section(
    class = "ytab-qc-section",
    tags$h4("2. Configure and run"),
    radioButtons("library_diagnostics_mode", "Execution mode",
                 choices = list("Preview command" = "preview", "Run diagnostics" = "run"),
                 selected = "preview", inline = TRUE),
    uiOutput("qc_execution_summary"),
    uiOutput("library_diagnostics_cache_ui"),
    uiOutput("library_diagnostics_action")
  ),
  tags$section(
    class = "ytab-qc-section",
    tags$h4("3. Results"),
    ytab_two_column_layout(
      controls = ytab_control_panel(
        "Diagnostics display",
        selectInput("library_diagnostics_plot_choice", "Visualization",
                    choices = c("MidLC saturation" = "midlc",
                                "Jackpots and library depth" = "jackpot",
                                "Centromere bias" = "centromere",
                                "Feature metaplots" = "metaplots",
                                "Sequence bias" = "sequence_bias")),
        uiOutput("library_diagnostics_group_selector"),
        uiOutput("library_diagnostics_metaplot_selector"),
        selectInput("library_diagnostics_color_by", "Color by",
                    choices = qc_library_color_choices(), selected = "group"),
        ytab_plot_customization_controls("library_diagnostics", include_bars = TRUE,
                                         include_heatmap = TRUE,
                                         default_height = "medium",
                                         default_bar_orientation = "horizontal",
                                         default_show_value_labels = FALSE),
        tags$details(
          class = "ytab-more-options",
          tags$summary("Downloads"),
          tags$div(class = "ytab-actions",
                   downloadButton("download_library_diagnostics_plot", "Download plot"),
                   downloadButton("download_library_diagnostics_plotted_data", "Download plotted data"))
        ),
        tags$details(class = "ytab-technical-details", tags$summary("Technical details"),
                     tags$div(class = "ytab-technical-console",
                              verbatimTextOutput("library_diagnostics_technical")))
      ),
      main = ytab_plot_card("Library Diagnostics",
                            uiOutput("library_diagnostics_selected_plot"))
    )
  )
)
qc_files_ui <- function() panel_card(
  "Diagnostic Files",
  ytab_two_column_layout(
    controls = ytab_control_panel(
      "Diagnostic file display",
      actionButton("refresh_diagnostic_files", "Refresh diagnostic files", class = "btn-secondary"),
      uiOutput("diagnostic_result_selector"),
      uiOutput("diagnostic_file_filters"),
      checkboxInput("diagnostic_show_archived_static",
                    "Show archived/static generated image files", FALSE)
    ),
    main = ytab_plot_card(
      "Diagnostic files",
      tagList(
        uiOutput("diagnostic_files_empty"),
        DT::DTOutput("diagnostic_files_table")
      )
    )
  )
)
quality_control_ui <- function() navset_tab(id="qc_tabs",nav_panel("Mapping QC",value="mapping_qc",qc_mapping_ui()),nav_panel("Summary QC",value="summary_qc",qc_summary_ui()),nav_panel("Library Diagnostics",value="library_diagnostics",qc_library_diagnostics_ui()),nav_panel("Diagnostic Files",value="diagnostic_files",qc_files_ui()))
