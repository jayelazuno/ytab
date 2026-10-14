#!/usr/bin/env Rscript
# A missing optional catalogue must not disconnect the whole Shiny session.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
repo_root <- normalizePath(args[[1L]], mustWork=TRUE)
library(shiny)
library(bslib)
`%||%` <- function(x,y) if(is.null(x)||!length(x)) y else x
for(name in c("ui_helpers.R","ui_components.R","plot_customization_helpers.R",
              "plot_display_helpers.R","screen_coverage.R"))
  source(file.path(repo_root,"app/shiny/R",name))
missing_audit <- tempfile("ytab-missing-audit-")
active <- function() list()
active_project_path <- function() file.path(repo_root,"output/projects/H2O2_screen_v1/config/project.yaml")
python_bin <- function() "python"
test_server <- screen_coverage_server
formals(test_server) <- alist(input=,output=,session=)
expect_validation <- function(expr) {
  value <- tryCatch(force(expr), shiny.silent.error=function(e)e)
  stopifnot(inherits(value,"shiny.silent.error"))
}
withCallingHandlers(shiny::testServer(test_server, {
  session$setInputs(sc_audit_dir=missing_audit,sc_model_mode="poisson",
                    sc_target_preset="0.95",sc_allocation_mode="manual")
  session$elapse(1000)
  session$flushReact()
  expect_validation(request())
  expect_validation(output$sc_manual_controls)
  expect_validation(output$sc_poisson_target)
  expect_validation(output$sc_plot_primary)
  stopifnot(is.null(job$process),!session$isClosed())
  # The debounced observer remains safe when controls change after startup.
  session$setInputs(sc_target_preset="0.99")
  session$elapse(1000)
  session$flushReact()
  stopifnot(is.null(job$process),!session$isClosed())
}), warning=function(e)stop(conditionMessage(e),call.=FALSE))
cat("PASS: missing Screen Coverage audit leaves the Shiny session available\n")
