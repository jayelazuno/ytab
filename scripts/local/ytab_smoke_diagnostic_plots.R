#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
all<-commandArgs(trailingOnly=FALSE);self<-normalizePath(sub("^--file=","",grep("^--file=",all,value=TRUE)[[1L]]),winslash="/");root<-normalizePath(file.path(dirname(self),"../.."),winslash="/")
library(shiny);`%||%`<-function(left,right)if(is.null(left)||!length(left))right else left
source(file.path(root,"app/shiny/R/project_discovery.R"));source(file.path(root,"app/shiny/R/ui_qc.R"))
project<-read_project_summary(args[[1L]],root);inventory<-build_diagnostic_file_inventory(project);plots<-inventory[inventory$extension=="png",,drop=FALSE]
stopifnot(nrow(plots)>0L,all(nzchar(plots$relative_path)),all(!grepl("^/|file://",plots$served_url)),all(plots$preview_available),all(grepl("^ytab-diagnostics-(project|export)/",plots$served_url)))
routes<-diagnostic_resource_roots(project);for(route in names(routes)){suppressWarnings(try(removeResourcePath(route),silent=TRUE));addResourcePath(route,routes[[route]])};stopifnot(all(names(routes)%in%names(resourcePaths())))
space_root<-file.path(tempdir(),"YTAB diagnostics with spaces");dir.create(file.path(space_root,"library_diagnostics","runs","all samples"),recursive=TRUE,showWarnings=FALSE);space_project<-project;space_project$project_root<-space_root;space_project$export_root<-file.path(space_root,"exports");png<-file.path(space_root,"library_diagnostics","runs","all samples","plot with spaces.png");file.create(png);served<-diagnostic_served_file(png,space_project);stopifnot(served$preview_available,grepl("%20",served$served_url,fixed=TRUE),!grepl("file://",served$served_url,fixed=TRUE))
missing<-diagnostic_served_file(file.path(space_root,"library_diagnostics","missing.png"),space_project);stopifnot(!missing$preview_available,grepl("Preview unavailable",missing$preview_message,fixed=TRUE))
stopifnot(all(c("file","is_plot","served_url","relative_path")%in%names(inventory)))
cat("PASS\n")
