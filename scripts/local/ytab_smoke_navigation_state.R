#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
all<-commandArgs(trailingOnly=FALSE);self<-normalizePath(sub("^--file=","",grep("^--file=",all,value=TRUE)[[1L]]),winslash="/");root<-normalizePath(file.path(dirname(self),"../.."),winslash="/")
`%||%`<-function(left,right)if(is.null(left)||!length(left))right else left;source(file.path(root,"app/shiny/R/navigation.R"))
state<-navigation_state_defaults("run");state<-record_navigation_selection(state,"run_tabs","hit_files");stopifnot(state$active_run_subtab=="hit_files");unchanged<-state;stopifnot(identical(unchanged,state))
state<-record_navigation_selection(state,"workspace_tabs","quality_control");state<-record_navigation_selection(state,"qc_tabs","library_correlation");stopifnot(state$active_qc_subtab=="agreement")
for(old in names(navigation_qc_legacy))stopifnot(navigation_target_from_hash(paste0("#tab=quality_control&subtab=",old))$nested_value==navigation_qc_legacy[[old]])
stopifnot(navigation_state_defaults("quality_control")$active_qc_subtab=="overview")
state<-record_navigation_selection(state,"workspace_tabs","essentiality");state<-record_navigation_selection(state,"essentiality_tabs","combined_summary");stopifnot(state$active_essentiality_subtab=="run");restored<-normalize_navigation_state(state);stopifnot(restored$active_essentiality_subtab=="run")
state<-record_navigation_selection(state,"essentiality_tabs","results");stopifnot(state$active_essentiality_subtab=="results",navigation_state_defaults()$active_essentiality_subtab=="run")
legacy<-navigation_target_from_hash("#tab=essentiality&subtab=classifier");stopifnot(legacy$top=="essentiality",legacy$nested_value=="run")
stopifnot(navigation_target_from_hash("#tab=essentiality&subtab=hybrid_inputs")$nested_value=="run")
state<-record_navigation_selection(state,"workspace_tabs","fitness");state<-record_navigation_selection(state,"fitness_tabs","fitness_run");stopifnot(state$active_fitness_subtab=="fitness_run");state<-record_navigation_selection(state,"essentiality_results_view","predictions");stopifnot(state$active_results_view=="predictions")
state<-record_navigation_selection(state,"workspace_tabs","go_analysis");state<-record_navigation_selection(state,"go_analysis_tabs","semantic_reduction");stopifnot(state$active_top_tab=="gene_explorer",state$active_go_analysis_subtab=="semantic_reduction")
# GO and Across Species keep inner tabs inside Gene Explorer; old links open the inner tab.
legacy<-navigation_target_from_hash("#tab=go_analysis&subtab=enrichment");stopifnot(legacy$top=="gene_explorer",legacy$nested_value=="go",legacy$inner_id=="go_analysis_tabs",legacy$inner_value=="enrichment")
legacy<-navigation_state_from_url(list(top="comparative_species",nested_value="gene_group"));stopifnot(legacy$active_gene_explorer_subtab=="comparative",legacy$active_comparative_subtab=="gene_group")
legacy<-navigation_target_from_hash("#tab=genome_browser");stopifnot(legacy$top=="gene_explorer",legacy$nested_value=="browser")
stopifnot(navigation_state_defaults()$active_gene_explorer_subtab=="browser","go_analysis_tabs"%in%names(navigation_nested_fields))
app<-paste(readLines(file.path(root,"app/shiny/app.R"),warn=FALSE),collapse="\n")
nav_js<-paste(readLines(file.path(root,"app/shiny/www/ytab_nav_menus.js"),warn=FALSE),collapse="\n")
stopifnot(
  grepl("restore_navigation_state",app,fixed=TRUE),
  grepl("record_navigation_selection",app,fixed=TRUE),
  grepl("navigation_history_restore <- reactiveVal(FALSE)",app,fixed=TRUE),
  grepl("!isTRUE(navigation_history_restore())",app,fixed=TRUE),
  grepl("sync_navigation_url(session, location, replace = TRUE",app,fixed=TRUE),
  grepl("ytab_nav_menus.js?v=20261022_2",app,fixed=TRUE),
  grepl("function restoreFromBrowserHistory()",nav_js,fixed=TRUE),
  !grepl('classList.add("dropdown"',nav_js,fixed=TRUE),
  grepl('window.addEventListener("popstate", restoreFromBrowserHistory)',nav_js,fixed=TRUE),
  grepl('window.addEventListener("hashchange", restoreFromBrowserHistory)',nav_js,fixed=TRUE),
  grepl('sendShinyValue("ytab_nav_hash", target);',nav_js,fixed=TRUE),
  grepl('window.setTimeout(function () {\n      clientNavigate(target.top, target.nested_value);\n    }, 180);',nav_js,fixed=TRUE),
  !grepl('clientNavigate(target.top, target.nested_value);\n    sendShinyValue("ytab_nav_hash", target);',nav_js,fixed=TRUE),
  grepl("initial_navigation_probe <- reactiveVal(0L)",app,fixed=TRUE),
  grepl("invalidateLater(100, session)",app,fixed=TRUE),
  !grepl("state$load(active_project_path(),repo_root);go_to",app,fixed=TRUE)
)
landing<-navigation_state_defaults("project");stopifnot(landing$active_top_tab=="project",landing$active_project_subtab=="samples",landing$active_run_subtab=="samples")
# Retired tabs (old links and bookmarks) open their new module.
legacy<-navigation_target_from_hash("#tab=preprocessing&subtab=hit_files");stopifnot(legacy$top=="run",legacy$nested_id=="run_tabs",legacy$nested_value=="mapping")
for(old in c("hit_files","summary_tables","status"))stopifnot(navigation_target_from_hash(paste0("#tab=run&subtab=",old))$nested_value=="mapping")
legacy<-navigation_target_from_hash("#tab=logs_provenance");stopifnot(legacy$top=="project",legacy$nested_value=="provenance")
legacy<-navigation_target_from_hash("#tab=project_data&subtab=reference");stopifnot(legacy$top=="run",legacy$nested_value=="reference")
legacy<-navigation_state_from_url(list(top="pipeline"));stopifnot(legacy$active_top_tab=="run",legacy$active_run_subtab=="profiles")
stopifnot(normalize_navigation_state(list(active_top_tab="results_exports"))$active_top_tab=="project")
cat("PASS\n")
