#!/usr/bin/env Rscript
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1]]), winslash = "/")
repo_root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(repo_root, "app/shiny/R/project_discovery.R"), local = TRUE)

failures <- character()
check <- function(ok, label) if (!isTRUE(ok)) failures <<- c(failures, label)
visible <- discover_ytab_projects(repo_root)
all_projects <- discover_ytab_projects(repo_root, include_tests = TRUE)
labels <- if (nrow(all_projects)) vapply(seq_len(nrow(all_projects)), function(i) project_display_label(as.list(all_projects[i, ])), "") else character()
smoke_id <- "app_init_command_smoke_20260720114148"

check(nrow(visible) > 0L, "valid user projects are discovered")
check(all(labels %in% all_projects$display_name), "option labels contain display name only")
check(!any(grepl("glabrata — analysis ready", labels, fixed = TRUE)), "labels omit species and status")
check(all(file.exists(all_projects$project_config)) && all(basename(all_projects$project_config) == "project.yaml"), "exact project.yaml paths remain values")
check(!smoke_id %in% visible$project_id, "test projects are hidden by default")
check(smoke_id %in% all_projects$project_id, "test projects appear when requested")
configured <- select_demo_project(all_projects, "local_glabrata_smoke_v1")
check(!is.null(configured) && identical(configured$project_id[[1]], "local_glabrata_smoke_v1"), "configured demo is selected")
check(!is.null(configured) && !startsWith(configured$project_id[[1]], "app_init_command_smoke_"), "internal smoke project is not demo")
check(is.null(select_demo_project(all_projects, "missing_demo_project")), "missing demo has no random replacement")
check(!any(grepl(repo_root, labels, fixed = TRUE)), "labels contain no absolute paths")

current <- if (nrow(visible)) visible$project_config[[1]] else ""
refreshed <- discover_ytab_projects(repo_root)
preserved <- if (current %in% refreshed$project_config) current else ""
check(identical(preserved, current), "refresh preserves a valid selection")

if (length(failures)) { cat("FAIL\n", paste("-", failures, collapse = "\n"), "\n", sep = ""); quit(status = 1L) }
cat("PASS\n")
