essentiality_stage_script <- function(stage) {
  scripts <- c(
    sample_normalization = "ytab_run_sample_normalization.py",
    summary_normalized = "ytab_run_summary_normalized.py",
    combined_hits = "ytab_run_combine_hits.py",
    summary_combined = "ytab_run_summary_combined.py",
    classifier = "ytab_run_classifier.py"
  )
  script <- unname(scripts[[stage]])
  if (is.null(script)) stop("Unknown Essentiality stage: ", stage, call. = FALSE)
  script
}

essentiality_optional_number_arg <- function(flag, value) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) == 1L && !is.na(value) && is.finite(value) && value > 0)
    c(flag, essentiality_numeric_text(value)) else character()
}

build_essentiality_command <- function(repo_root, project_config, stage,
                                       target = "", parents = character(),
                                       execution_mode = "preview", force = FALSE,
                                       threads = 2L, min_site_retention = .95,
                                       min_feature_retention = .95,
                                       auto_min = NA_real_, auto_max = NA_real_,
                                       auto_step = 5, seed = 0L,
                                       keep_going = FALSE) {
  stopifnot(execution_mode %in% c("preview", "run"))
  script <- file.path(repo_root, "scripts", "local", essentiality_stage_script(stage))
  base <- c(script, "--project-config", normalizePath(project_config, winslash = "/", mustWork = TRUE))
  selected <- if (length(parents)) c("--samples", paste(unique(parents), collapse = ",")) else character()
  args <- switch(stage,
    sample_normalization = c(
      "--targets", target, "--sample-mode", "parents", selected,
      "--threads", as.character(as.integer(threads)),
      "--min-site-retention", as.character(min_site_retention),
      if (identical(target, "auto")) essentiality_optional_number_arg("--auto-min-target", auto_min),
      if (identical(target, "auto")) essentiality_optional_number_arg("--auto-max-target", auto_max),
      if (identical(target, "auto")) c("--auto-step", essentiality_numeric_text(auto_step)),
      if (keep_going) "--keep-going"
    ),
    summary_normalized = c(
      "--targets", target, "--sample-mode", "parents", selected,
      "--threads", as.character(as.integer(threads)),
      "--min-feature-retention", as.character(min_feature_retention),
      if (keep_going) "--keep-going"
    ),
    combined_hits = c("--target", target, selected),
    summary_combined = c("--target", target, "--threads", as.character(as.integer(threads))),
    classifier = c("--target", target, "--seed", as.character(as.integer(seed)))
  )
  args <- c(base, args, if (identical(execution_mode, "preview")) "--dry-run",
            if (identical(execution_mode, "run") && isTRUE(force)) "--force")
  list(command = script, args = args[-1L], full = args, stage = stage,
       project_config = normalizePath(project_config, winslash = "/", mustWork = TRUE),
       execution_mode = execution_mode, force = identical(execution_mode, "run") && isTRUE(force),
       parents = unique(parents), target = target)
}
