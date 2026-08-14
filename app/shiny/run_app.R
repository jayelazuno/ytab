#!/usr/bin/env Rscript
if (!requireNamespace("shiny", quietly = TRUE)) stop("The R package 'shiny' is required. Activate ytab-local.", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE); opts <- list(host="127.0.0.1",port="3838",project_config="",browser=TRUE)
i <- 1L
while (i <= length(args)) {
  if (args[[i]] == "--no-browser") { opts$browser <- FALSE; i <- i+1L; next }
  if (i == length(args) || !args[[i]] %in% c("--host","--port","--project-config")) stop("Unknown or incomplete argument: ",args[[i]])
  key <- gsub("-","_",substring(args[[i]],3)); opts[[key]] <- args[[i+1L]]; i <- i+2L
}
all_args <- commandArgs(trailingOnly=FALSE); file_arg <- grep("^--file=",all_args,value=TRUE)
launcher <- if(length(file_arg)) sub("^--file=","",file_arg[[1]]) else "app/shiny/run_app.R"
# Rscript encodes spaces in --file on some platforms.
launcher <- gsub("~\\+~", " ", launcher)
app_dir <- normalizePath(dirname(launcher),mustWork=TRUE); repo_root <- normalizePath(file.path(app_dir,"../.."),mustWork=TRUE)
if(nzchar(opts$project_config)) {
  path <- opts$project_config; if(!grepl("^(/|[A-Za-z]:)",path)) path <- file.path(repo_root,path)
  candidate <- normalizePath(path, winslash="/", mustWork=FALSE)
  warning_text <- if(!file.exists(candidate)) paste("Project configuration does not exist:",candidate) else if(dir.exists(candidate)) paste("Expected a project.yaml file, but received a directory:",candidate) else if(basename(candidate)!="project.yaml") paste("Expected a file named project.yaml, received:",basename(candidate)) else ""
  if(nzchar(warning_text)){warning(warning_text,call.=FALSE);Sys.setenv(YTAB_PROJECT_CONFIG_INVALID=warning_text);Sys.unsetenv("YTAB_PROJECT_CONFIG")} else Sys.setenv(YTAB_PROJECT_CONFIG=normalizePath(candidate,winslash="/",mustWork=TRUE))
}
url <- sprintf("http://%s:%d",opts$host,as.integer(opts$port)); cat("YTAB app URL:",url,"\n")
shiny::runApp(appDir=app_dir,host=opts$host,port=as.integer(opts$port),launch.browser=isTRUE(opts$browser))
