#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
all<-commandArgs(trailingOnly=FALSE);self<-normalizePath(sub("^--file=","",grep("^--file=",all,value=TRUE)[[1L]]),winslash="/");root<-normalizePath(file.path(dirname(self),"../.."),winslash="/");config<-normalizePath(args[[1L]],winslash="/")
options(shiny.fullstacktrace=TRUE);library(shiny);library(bslib)
app_env<-new.env(parent=.GlobalEnv);source(file.path(root,"app/shiny/app.R"),local=app_env)

landing<-app_env$landing_ui(character(),2L,FALSE);workspace<-app_env$workspace_ui(2L)
landing_html<-htmltools::renderTags(landing)$html;workspace_html<-htmltools::renderTags(workspace)$html
stopifnot(grepl("Open an existing YTAB project",landing_html,fixed=TRUE),grepl("workspace_tabs",workspace_html,fixed=TRUE))
for(label in c("Overview","Preprocessing","Quality Control","Essentiality","Fitness Screen","Results &amp; Exports"))stopifnot(grepl(label,workspace_html,fixed=TRUE))
static_code<-paste(deparse(body(app_env$workspace_ui)),deparse(body(app_env$landing_ui)),collapse="\n")
stopifnot(!grepl("active\\(\\)|active_project\\(\\)|state\\$values|input\\$",static_code))

server_text<-paste(readLines(file.path(root,"app/shiny/app.R"),warn=FALSE),collapse="\n")
stopifnot(!grepl("onFlushed\\(function\\(\\) navigation_location",server_text),!grepl("onFlushed\\(function\\(\\)[^{]*\\{[^}]*navigation_location\\(",server_text,perl=TRUE))

testServer(function(input,output,session){
  project_state<-app_env$new_project_state();location<-reactiveVal(app_env$navigation_state_defaults("overview"))
  observeEvent(input$open,{project_state$load(config,root);project_state$enter();snapshot<-app_env$navigation_state_defaults(isolate(project_state$values$workspace_tab));location(snapshot);session$onFlushed(function()app_env$restore_navigation_state(session,snapshot),once=TRUE)})
  observeEvent(input$switch,{project_state$clear();location(app_env$navigation_state_defaults("overview"))})
  session$userData$project_state<-project_state;session$userData$location<-location
},{
  session$setInputs(open=1);session$flushReact();stopifnot(session$userData$project_state$values$view=="workspace",!is.null(session$userData$project_state$values$project))
  session$setInputs(switch=1);session$flushReact();stopifnot(session$userData$project_state$values$view=="landing",is.null(session$userData$project_state$values$project),session$userData$location()$active_top_tab=="overview")
})
cat("PASS\n")
