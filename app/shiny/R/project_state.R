new_project_state <- function() {
  values <- reactiveValues(view = "landing", path = NULL, project = NULL, status = NULL, landing_warning = "", workspace_tab = "overview")
  load_project <- function(path, repo_root) {
    checked <- validate_project_config_path(path, repo_root)
    if (!checked$valid) stop(checked$message, call. = FALSE)
    project <- read_project_summary(checked$path, repo_root)
    values$path <- checked$path; values$project <- project; values$status <- project$status
    invisible(project)
  }
  list(
    values = values,
    load = load_project,
    clear = function() { values$path <- NULL; values$project <- NULL; values$status <- NULL; values$view <- "landing" },
    refresh = function(repo_root) { req(values$path); load_project(values$path, repo_root) },
    enter = function() { req(values$project); values$view <- "workspace"; values$workspace_tab <- if (isTRUE(values$project$analysis_ready)) "overview" else "preprocessing" },
    path = reactive(values$path), project = reactive(values$project), status = reactive(values$status), view = reactive(values$view)
  )
}
