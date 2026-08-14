ytab_resolve_path <- function(path, base) {
  if (is.null(path) || !length(path) || !nzchar(trimws(path[[1]]))) return("")
  path <- path[[1]]
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(base, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

validate_project_config_path <- function(path, repo_root) {
  fail <- function(message) list(valid = FALSE, path = NULL, data = NULL, message = message)
  if (is.null(path) || !length(path) || !nzchar(trimws(path[[1]]))) return(fail("Select a project.yaml file."))
  candidate <- ytab_resolve_path(path[[1]], repo_root)
  if (!file.exists(candidate)) return(fail(paste("Project configuration does not exist:", candidate)))
  if (dir.exists(candidate)) return(fail(paste("Expected a project.yaml file, but received a directory:", candidate)))
  info <- file.info(candidate)
  if (is.na(info$isdir) || info$isdir) return(fail(paste("Expected a regular project.yaml file:", candidate)))
  if (!identical(basename(candidate), "project.yaml")) return(fail(paste("Expected a file named project.yaml, received:", basename(candidate))))
  if (!requireNamespace("yaml", quietly = TRUE)) return(fail("The R package yaml is required to read project configuration."))
  config <- tryCatch(yaml::read_yaml(candidate), error = function(e) e)
  if (inherits(config, "error")) return(fail(paste("Could not read project YAML:", conditionMessage(config))))
  project_id <- trimws(as.character(config$project_id %||% ""))
  if (!nzchar(project_id)) return(fail("The project ID could not be resolved from project.yaml."))
  project_root <- normalizePath(dirname(dirname(candidate)), winslash = "/", mustWork = TRUE)
  sample_sheet <- ytab_resolve_path(config$sample_sheet %||% file.path("output", "projects", project_id, "config", "sample_sheet.csv"), repo_root)
  if (!file.exists(sample_sheet) || dir.exists(sample_sheet)) return(fail(paste("The referenced sample sheet could not be resolved:", sample_sheet)))
  expected_config <- normalizePath(file.path(project_root, "config"), winslash = "/", mustWork = TRUE)
  if (!startsWith(sample_sheet, paste0(expected_config, "/"))) return(fail("The sample sheet must remain inside the project's config directory."))
  list(valid = TRUE, path = normalizePath(candidate, winslash = "/", mustWork = TRUE), data = config,
       project_root = project_root, sample_sheet = sample_sheet, message = "")
}

read_project_status <- function(project_root) {
  path <- file.path(project_root, "manifests", "project_status.json")
  if (!file.exists(path) || !requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
}

stage_status <- function(status, stage) {
  if (is.null(status$stages)) return("not calculated")
  hit <- Filter(function(x) identical(x$stage %||% "", stage), status$stages)
  if (!length(hit)) "not calculated" else hit[[1]]$status %||% "not calculated"
}

project_is_analysis_ready <- function(info, status = info$status) {
  identical(stage_status(status, "summary"), "success")
}

read_project_summary <- function(project_config, repo_root) {
  checked <- validate_project_config_path(project_config, repo_root)
  if (!checked$valid) stop(checked$message, call. = FALSE)
  config <- checked$data
  samples <- tryCatch(read.csv(checked$sample_sheet, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  include <- if ("include" %in% names(samples)) tolower(as.character(samples$include)) %in% c("true", "1", "yes") else rep(TRUE, nrow(samples))
  condition <- if ("condition" %in% names(samples)) {
    tolower(samples$condition)
  } else if ("guessed_condition" %in% names(samples)) {
    tolower(samples$guessed_condition)
  } else rep("", nrow(samples))
  status <- read_project_status(checked$project_root)
  project_id <- as.character(config$project_id)
  fallback_test <- grepl("^(app_init_command_smoke_|progress_test_|app_smoke_|ytab_test_|Test_run|test_run)", project_id)
  project_type <- tolower(trimws(as.character(config$project_type %||% if (fallback_test) "test" else "user")))
  if (!project_type %in% c("user", "demo", "test", "temporary", "release")) project_type <- "user"
  show_in_launcher <- config$show_in_launcher %||% TRUE
  show_in_launcher <- isTRUE(show_in_launcher) || identical(tolower(as.character(show_in_launcher)), "true")
  display_name <- trimws(as.character(config$display_name %||% project_id))
  if (!nzchar(display_name)) display_name <- project_id
  info <- list(project_config = checked$path, project_id = project_id, display_name = display_name,
    project_type = project_type, show_in_launcher = show_in_launcher, species = as.character(config$species %||% "unknown"),
    start_stage = as.character(config$start_stage %||% "fastq"),
    fastq_directory = ytab_resolve_path(config$fastq_dir %||% "", repo_root), sample_sheet = checked$sample_sheet,
    reference = config$reference %||% list(), threads = as.integer(config$threads %||% 2L), repo_root = repo_root,
    project_root = checked$project_root,
    export_root = ytab_resolve_path(config$output_export_dir %||% file.path("output", "exports", config$project_id), repo_root),
    samples = samples, sample_count = nrow(samples), included_count = sum(include), parent_count = sum(include & condition == "parent"),
    treated_count = sum(include & condition == "treated"), status = status, modified = file.info(checked$path)$mtime)
  info$raw_summary_complete <- project_is_analysis_ready(info, status)
  info$analysis_ready <- info$raw_summary_complete
  info
}

project_display_label <- function(info) {
  label <- trimws(as.character(info$display_name %||% info$project_id %||% ""))
  if (nzchar(label)) label else as.character(info$project_id %||% "")
}

discover_ytab_projects <- function(repo_root, include_tests = FALSE) {
  base <- file.path(repo_root, "output", "projects")
  paths <- if (dir.exists(base)) list.files(base, pattern = "^project\\.yaml$", recursive = TRUE, full.names = TRUE) else character()
  paths <- paths[grepl("/config/project\\.yaml$", gsub("\\\\", "/", paths))]
  records <- lapply(paths, function(path) tryCatch(read_project_summary(path, repo_root), error = function(e) NULL))
  records <- Filter(Negate(is.null), records)
  records <- Filter(function(x) isTRUE(x$show_in_launcher) &&
    (isTRUE(include_tests) || !x$project_type %in% c("test", "temporary")), records)
  if (!length(records)) return(data.frame())
  data.frame(project_config = vapply(records, `[[`, "", "project_config"),
    project_id = vapply(records, `[[`, "", "project_id"), display_name = vapply(records, `[[`, "", "display_name"),
    project_type = vapply(records, `[[`, "", "project_type"), species = vapply(records, `[[`, "", "species"),
    sample_count = vapply(records, `[[`, 0L, "sample_count"), included_count = vapply(records, `[[`, 0L, "included_count"),
    parent_count = vapply(records, `[[`, 0L, "parent_count"), treated_count = vapply(records, `[[`, 0L, "treated_count"),
    analysis_ready = vapply(records, `[[`, FALSE, "analysis_ready"),
    status = vapply(records, function(x) if (is.null(x$status)) "status not yet calculated" else stage_status(x$status, "summary"), ""),
    modified = as.POSIXct(vapply(records, function(x) as.numeric(x$modified), 0)), stringsAsFactors = FALSE)
}

read_launcher_config <- function(app_dir) {
  path <- file.path(app_dir, "config.yaml")
  if (!file.exists(path) || !requireNamespace("yaml", quietly = TRUE)) return(list())
  tryCatch(yaml::read_yaml(path) %||% list(), error = function(e) list())
}

select_demo_project <- function(projects, configured_id = "") {
  if (!nrow(projects)) return(NULL)
  explicit <- which(projects$project_type == "demo")
  if (length(explicit)) return(projects[explicit[[1]], , drop = FALSE])
  configured_id <- trimws(as.character(configured_id %||% ""))
  configured <- which(nzchar(configured_id) & projects$project_id == configured_id &
    !projects$project_type %in% c("test", "temporary"))
  if (length(configured)) projects[configured[[1]], , drop = FALSE] else NULL
}

scan_fastq_directory <- function(path, limit = 5000L) {
  if (!nzchar(path) || !dir.exists(path)) return(data.frame())
  files <- list.files(path, pattern = "\\.(fastq|fq)(\\.gz)?$", ignore.case = TRUE, full.names = TRUE)
  if (!length(files)) return(data.frame())
  truncated <- length(files) > limit
  if (truncated) files <- files[seq_len(limit)]
  meta <- file.info(files)
  filename <- basename(files)
  data.frame(filename = filename,
    inferred_sample = sub("(_R?[12]|[._-][12])?(\\.fastq|\\.fq)(\\.gz)?$", "", filename, ignore.case = TRUE),
    extension = sub("^.*?(\\.(fastq|fq)(\\.gz)?)$", "\\1", filename, ignore.case = TRUE),
    compressed = ifelse(grepl("\\.gz$", filename, ignore.case = TRUE), "yes", "no"),
    size = format(meta$size, big.mark = ","), modified = format(meta$mtime, "%Y-%m-%d %H:%M"), include = "yes",
    stringsAsFactors = FALSE) |> structure(truncated = truncated)
}
