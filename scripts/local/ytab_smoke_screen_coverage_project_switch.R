#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
self <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[[1]]))
root <- normalizePath(file.path(dirname(self),"../.."))
library(shiny);library(bslib)
`%||%` <- function(x,y) if(is.null(x)||!length(x)) y else x
for (name in c("project_discovery.R","ui_helpers.R","ui_components.R","plot_customization_helpers.R",
  "plot_display_helpers.R","analysis_runs.R","screen_planning.R")) source(file.path(root,"app/shiny/R",name))
paths <- vapply(args,normalizePath,"")
python <- Sys.getenv("YTAB_PYTHON",unset="/Users/joshuaayelazuno/miniforge3/envs/ytab-local/bin/python")
server <- function(input,output,session) {
  path <- reactive(input$test_project %||% paths[[1]])
  project <- reactive(read_project_summary(path(),root))
  screen_coverage_server(input,output,session,root,project,path,function()python)
}
shiny::testServer(server, {
  for (path in paths) {
    session$setInputs(test_project=path,sc_target_preset="0.95")
    session$flushReact()
    expected <- yaml::read_yaml(path)$project_id
    deadline <- Sys.time()+90
    repeat {
      session$elapse(350);Sys.sleep(.05);session$flushReact()
      rendered <- tryCatch(output$sc_primary_summary,shiny.silent.error=function(e)NULL)
      if (!is.null(rendered)) break
      if (Sys.time()>deadline) stop("Timed out loading project inventory: ",expected)
    }
    html <- paste(unlist(rendered),collapse=" ")
    stopifnot(grepl(paste("Project:",expected),html,fixed=TRUE),length(output$sc_plot_primary)>0,
      length(output$sc_poisson_target)>0)
    other <- setdiff(paths,path)
    for (other_path in other) stopifnot(!grepl(paste("Project:",yaml::read_yaml(other_path)$project_id),html,fixed=TRUE))
    cat("PASS: project-scoped summary and analytical plot:",expected,"\n")
  }
})
