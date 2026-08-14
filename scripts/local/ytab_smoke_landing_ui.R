#!/usr/bin/env Rscript
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1]]), winslash = "/")
repo_root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
library(shiny)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(repo_root, "app/shiny/R/ui_helpers.R"), local = TRUE)
source(file.path(repo_root, "app/shiny/R/ui_landing.R"), local = TRUE)
source(file.path(repo_root, "app/shiny/R/project_discovery.R"), local = TRUE)

html <- as.character(landing_ui(c("glabrata"), 2L, FALSE))
app_text <- paste(readLines(file.path(repo_root, "app/shiny/app.R"), warn = FALSE), collapse = "\n")
projects <- discover_ytab_projects(repo_root, include_tests = TRUE)
demo <- select_demo_project(projects, read_launcher_config(file.path(repo_root, "app/shiny"))$demo_project_id %||% "")
failures <- character()
check <- function(ok, label) if (!isTRUE(ok)) failures <<- c(failures, label)

check(length(gregexpr(">YTAB<", html, fixed = TRUE)[[1]]) == 1L, "one YTAB hero title")
check(length(gregexpr("Yeast Transposon Analysis Browser", html, fixed = TRUE)[[1]]) == 1L, "one browser subtitle")
check(!grepl("visually-hidden.*Create new project", app_text), "no duplicate top-level plain text")
check(grepl("Create new project", html, fixed = TRUE), "new-project action")
check(grepl("Select a project…", html, fixed = TRUE), "selector placeholder")
check(grepl("Open project", app_text, fixed = TRUE), "open-project action")
check(grepl("ytab-project-preview", app_text, fixed = TRUE), "compact metadata preview")
check(!is.null(demo) && !demo$project_type[[1]] %in% c("test", "temporary"), "demo is not internal test")
check(grepl("View setup guide", html, fixed = TRUE) && grepl("Show launch command", html, fixed = TRUE), "documentation actions")
check(grepl("btn-primary ytab-button", html, fixed = TRUE) && grepl("btn-secondary ytab-button", html, fixed = TRUE), "consistent button classes")
check(!grepl(repo_root, html, fixed = TRUE), "no repository paths in visible landing text")

if (length(failures)) { cat("FAIL\n", paste("-", failures, collapse = "\n"), "\n", sep = ""); quit(status = 1L) }
cat("PASS\n")
