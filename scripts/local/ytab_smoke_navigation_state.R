#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
all<-commandArgs(trailingOnly=FALSE);self<-normalizePath(sub("^--file=","",grep("^--file=",all,value=TRUE)[[1L]]),winslash="/");root<-normalizePath(file.path(dirname(self),"../.."),winslash="/")
`%||%`<-function(left,right)if(is.null(left)||!length(left))right else left;source(file.path(root,"app/shiny/R/navigation.R"))
state<-navigation_state_defaults("preprocessing");state<-record_navigation_selection(state,"preprocessing_tabs","hit_files");stopifnot(state$active_preprocessing_subtab=="hit_files");unchanged<-state;stopifnot(identical(unchanged,state))
state<-record_navigation_selection(state,"workspace_tabs","quality_control");state<-record_navigation_selection(state,"qc_tabs","diagnostic_files");stopifnot(state$active_qc_subtab=="diagnostic_files")
state<-record_navigation_selection(state,"workspace_tabs","essentiality");state<-record_navigation_selection(state,"essentiality_tabs","combined_summary");stopifnot(state$active_essentiality_subtab=="combined_summary");restored<-normalize_navigation_state(state);stopifnot(restored$active_essentiality_subtab=="combined_summary")
state<-record_navigation_selection(state,"workspace_tabs","fitness");state<-record_navigation_selection(state,"fitness_tabs","fitness_run");stopifnot(state$active_fitness_subtab=="fitness_run");state<-record_navigation_selection(state,"essentiality_results_view","predictions");stopifnot(state$active_results_view=="predictions")
state<-record_navigation_selection(state,"workspace_tabs","go_analysis");state<-record_navigation_selection(state,"go_analysis_tabs","semantic_reduction");stopifnot(state$active_go_analysis_subtab=="semantic_reduction")
app<-paste(readLines(file.path(root,"app/shiny/app.R"),warn=FALSE),collapse="\n")
nav_js<-paste(readLines(file.path(root,"app/shiny/www/ytab_nav_menus.js"),warn=FALSE),collapse="\n")
stopifnot(
  grepl("restore_navigation_state",app,fixed=TRUE),
  grepl("record_navigation_selection",app,fixed=TRUE),
  grepl("navigation_history_restore <- reactiveVal(FALSE)",app,fixed=TRUE),
  grepl("!isTRUE(navigation_history_restore())",app,fixed=TRUE),
  grepl("sync_navigation_url(session, location, replace = TRUE",app,fixed=TRUE),
  grepl("ytab_nav_menus.js?v=20260909_1",app,fixed=TRUE),
  grepl("function restoreFromBrowserHistory()",nav_js,fixed=TRUE),
  grepl('window.addEventListener("popstate", restoreFromBrowserHistory)',nav_js,fixed=TRUE),
  grepl('window.addEventListener("hashchange", restoreFromBrowserHistory)',nav_js,fixed=TRUE),
  grepl('sendShinyValue("ytab_nav_hash", target);',nav_js,fixed=TRUE),
  !grepl('clientNavigate(target.top, target.nested_value);\n    sendShinyValue("ytab_nav_hash", target);',nav_js,fixed=TRUE),
  !grepl("state$load(active_project_path(),repo_root);go_to",app,fixed=TRUE)
)
landing<-navigation_state_defaults("overview");stopifnot(landing$active_top_tab=="overview",landing$active_preprocessing_subtab=="samples")
cat("PASS\n")
