#!/usr/bin/env Rscript
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]]))
root <- normalizePath(file.path(dirname(self), "../.."))
library(shiny)
`%||%` <- function(x,y) if(is.null(x)||!length(x)) y else x
source(file.path(root,"app/shiny/R/project_discovery.R"))
source(file.path(root,"app/shiny/R/preprocessing_status.R"))
source(file.path(root,"app/shiny/R/sample_selector.R"))
samples <- data.frame(sample=c("read-A","read-B"),display_name=c("Control A","Drug B"),
  guessed_condition=c("parent","treated"),include=TRUE)
status <- data.frame(sample=samples$sample,included=TRUE,mapping=c("success","failed"),
  hit_file="success",summary="success")
registered <- new.env(parent=emptyenv())
registrations <- 0L
server <- function(input,output,session) {
  session$registerDataObj <- function(name,data,filterFunc) {
    if (!identical(filterFunc, getFromNamespace("selectizeJSON", "shiny")))
      return(paste0("session/test/dataobj/",name))
    registrations <<- registrations + 1L
    registered[[name]] <- list(data=data,filter=filterFunc)
    paste0("session/test/dataobj/",name)
  }
  session$userData$selected <- sample_selector_server("pick",reactive(samples),status_data=reactive(status))
}
check_requests <- function() {
  stopifnot(length(ls(registered))>0L)
  for (name in ls(registered)) {
    entry <- registered[[name]]
    stopifnot(all(c("label","value") %in% names(entry$data)), nrow(entry$data)==2L)
    for (query in c("", "Drug")) {
      request <- list(QUERY_STRING=paste0("query=",query,"&field=%5B%22label%22%5D&maxop=100&value=value&conju=and"))
      # Exercise the actual HTTP data-object filter used by Shiny's handleRequest.
      response <- entry$filter(entry$data,request)
      stopifnot(response$status==200L)
      decoded <- jsonlite::fromJSON(response$content)
      stopifnot(all(decoded$value %in% samples$sample))
    }
  }
}
shiny::testServer(server, {
  session$flushReact()
  check_requests()
  initial_registrations <- registrations
  for (button in c("none","all","parents","treated","incomplete","failed")) {
    do.call(session$setInputs,setNames(list(1),paste0("pick-",button)))
    session$flushReact()
    check_requests()
    stopifnot(registrations==initial_registrations)
  }
  session$setInputs(`pick-samples`="read-A")
  session$flushReact()
  check_requests()
  stopifnot(identical(session$userData$selected(),"read-A"))
})
cat("PASS: actual Selectize request callbacks after all selection presets and manual selection\n")
