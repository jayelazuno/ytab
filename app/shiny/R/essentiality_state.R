essentiality_truthy <- function(value) {
  tolower(trimws(as.character(value %||% ""))) %in% c("true", "1", "yes", "y", "parent")
}

essentiality_first_column <- function(data, candidates, default = "") {
  hit <- intersect(candidates, names(data))
  if (!length(hit)) rep(default, nrow(data)) else as.character(data[[hit[[1]]]])
}

essentiality_parent_inventory <- function(project) {
  data <- project$samples
  if (is.null(data) || !nrow(data)) {
    return(data.frame(Sample = character(), Background = character(), Pool = character(),
                      `Raw hit file` = character(), Eligible = character(), Reason = character(),
                      Explicit_treated = logical(), stringsAsFactors = FALSE, check.names = FALSE))
  }
  included <- if ("include" %in% names(data)) essentiality_truthy(data$include) else rep(TRUE, nrow(data))
  sample <- essentiality_first_column(data, c("sample", "Sample"))
  condition_columns <- intersect(c("condition", "guessed_condition", "treatment",
                                   "sample_role", "role"), names(data))
  parent_flag_columns <- intersect(c("parent", "is_parent", "parent_flag"), names(data))
  background <- essentiality_first_column(data, c("background", "guessed_background"), "")
  pool <- essentiality_first_column(data, c("pool", "guessed_pool"), "")
  rows <- lapply(seq_len(nrow(data)), function(index) {
    metadata <- tolower(trimws(unlist(lapply(condition_columns, function(column) as.character(data[[column]][[index]] %||% "")))))
    metadata <- metadata[nzchar(metadata)]
    explicit_treated <- any(grepl("treated|treatment|drug|stress", metadata)) &&
      !any(grepl("parent|control|untreated", metadata))
    flag_parent <- any(vapply(parent_flag_columns, function(column) essentiality_truthy(data[[column]][[index]]), logical(1)))
    metadata_parent <- any(grepl("parent", metadata))
    fallback_parent <- !length(metadata) && grepl("parent", sample[[index]], ignore.case = TRUE)
    detected <- included[[index]] && !explicit_treated && (flag_parent || metadata_parent || fallback_parent)
    reason <- if (!included[[index]]) "Excluded by project sample sheet" else if (explicit_treated) "Treated library excluded" else if (flag_parent) "Parent flag" else if (metadata_parent) "Parent metadata" else if (fallback_parent) "Parent name fallback" else "Parent role not identified"
    raw <- file.path("create_hit_file", sample[[index]], paste0(sample[[index]], "_hits.txt"))
    data.frame(
      Sample = sample[[index]], Background = background[[index]], Pool = pool[[index]],
      `Raw hit file` = raw,
      Eligible = if (detected && file.exists(file.path(project$project_root, raw))) "Yes" else if (detected) "No" else "No",
      Reason = if (detected && !file.exists(file.path(project$project_root, raw))) "Raw parent hit file missing" else reason,
      Detected_parent = detected, Explicit_treated = explicit_treated, Included = included[[index]],
      stringsAsFactors = FALSE, check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

detect_essentiality_parent_samples <- function(project, override = NULL) {
  inventory <- essentiality_parent_inventory(project)
  detected <- inventory$Sample[inventory$Detected_parent]
  if (is.null(override)) selected <- detected else {
    requested <- unique(as.character(override))
    allowed <- inventory$Sample[inventory$Included & !inventory$Explicit_treated]
    selected <- intersect(requested, allowed)
  }
  selected
}

essentiality_parent_review <- function(project, selected) {
  inventory <- essentiality_parent_inventory(project)
  review <- inventory[inventory$Sample %in% selected,
                      c("Sample", "Background", "Pool", "Raw hit file", "Eligible", "Reason"),
                      drop = FALSE]
  review[match(selected, review$Sample), , drop = FALSE]
}

essentiality_parent_override_choices <- function(project) {
  inventory <- essentiality_parent_inventory(project)
  choices <- inventory$Sample[inventory$Included & !inventory$Explicit_treated]
  as.list(setNames(unname(choices), unname(choices)))
}

essentiality_final_target <- function(project_root) {
  path <- file.path(project_root, "config", "final_classifier_target.txt")
  if (!file.exists(path) || dir.exists(path)) return("")
  value <- trimws(readLines(path, n = 1L, warn = FALSE))
  if (essentiality_valid_tag(value)) value else ""
}

essentiality_smoke_project <- function(project) {
  id <- tolower(as.character(project$project_id %||% ""))
  type <- tolower(as.character(project$project_type %||% ""))
  fastq <- tolower(as.character(project$fastq_directory %||% ""))
  identical(id, "h2o2-test2") || type %in% c("test", "temporary") ||
    grepl("smoke|test", id) || grepl("small_fastq_subset|smoke", fastq)
}

essentiality_normalized_inputs <- function(project_root, tag, parents) {
  paths <- file.path(project_root, "sample_normalization", tag, parents,
                     paste0(parents, "_normalized_hits.txt"))
  data.frame(Sample = parents, Path = paths,
             Available = file.exists(paths) & !dir.exists(paths) & file.info(paths)$size > 0,
             stringsAsFactors = FALSE)
}

essentiality_manifest <- function(project_root, stage, tag) {
  path <- switch(stage,
    normalize = file.path(project_root, "manifests", "sample_normalization",
                          paste0(tag, ".sample_normalization_manifest.json")),
    combine = file.path(project_root, "manifests", "combined_hits",
                        paste0(tag, ".combined_hits_manifest.json")),
    combined_summary = file.path(project_root, "manifests", "summary_combined",
                                 paste0(tag, ".summary_combined_manifest.json")),
    classifier = file.path(project_root, "manifests", "classifier",
                           paste0(tag, ".classifier_manifest.json")),
    ""
  )
  essentiality_read_json(path)
}

essentiality_classifier_resources <- function(project, repo_root, target_tag = "") {
  reference <- project$reference %||% list()
  resolve <- function(path) {
    path <- as.character(path %||% "")
    if (!nzchar(path)) return("")
    ytab_resolve_path(path, repo_root)
  }
  kornmann <- if (dir.exists(file.path(repo_root, "resources", "species", "Kornmann"))) {
    list.files(file.path(repo_root, "resources", "species", "Kornmann"),
               pattern = "WildType.*\\.wig$", full.names = TRUE)
  } else character()
  feature <- if (nzchar(target_tag)) file.path(project$project_root, "summary_combined", target_tag,
                                               paste0("combined_feature_table.", target_tag, ".txt")) else ""
  resources <- list(
    `BG2–S. cerevisiae orthology map` = resolve(reference$orthology_file),
    `Kornmann WildType WIG tracks` = kornmann,
    `SGD feature table` = file.path(repo_root, "resources", "species", "cerevisiae", "SGD_features.tab"),
    `Viable annotations` = file.path(repo_root, "resources", "species", "cerevisiae", "cerevisiae_viable_annotations.txt"),
    `Inviable annotations` = file.path(repo_root, "resources", "species", "cerevisiae", "cerevisiae_inviable_annotations.txt"),
    `S. cerevisiae paralog annotations` = file.path(repo_root, "resources", "species", "cerevisiae", "hasParalogs_sc.txt"),
    `Combined parent feature table` = feature
  )
  rows <- lapply(names(resources), function(name) {
    paths <- resources[[name]]
    available <- length(paths) > 0L && all(nzchar(paths)) && all(file.exists(paths)) &&
      all(!dir.exists(paths)) && all(file.info(paths)$size > 0)
    data.frame(Resource = name, Status = if (available) "Available" else "Missing",
               Requirement = "Required",
               Path = if (length(paths)) paste(paths, collapse = "; ") else "",
               stringsAsFactors = FALSE, check.names = FALSE)
  })
  result <- do.call(rbind, rows)
  attr(result, "supported_species") <- identical(as.character(project$species %||% ""), "glabrata")
  result
}

essentiality_stage_state <- function(project, parents, target_tag = "", job_progress = NULL,
                                     repo_root = dirname(dirname(project$project_root))) {
  root <- project$project_root
  available <- discover_essentiality_targets(root)
  normalized <- if (nzchar(target_tag) && length(parents)) essentiality_normalized_inputs(root, target_tag, parents) else data.frame()
  combined <- if (nzchar(target_tag)) file.path(root, "combined_hits", target_tag,
                                               paste0("combined_parent_hits.", target_tag, ".txt")) else ""
  feature <- if (nzchar(target_tag)) file.path(root, "summary_combined", target_tag,
                                              paste0("combined_feature_table.", target_tag, ".txt")) else ""
  prediction <- if (nzchar(target_tag)) file.path(root, "classifier", target_tag,
                                                 paste0("essentiality_predictions.", target_tag, ".csv")) else ""
  evaluation_state <- essentiality_feature_evaluation_state(root)
  stage_files <- list(
    normalize = nrow(available) > 0L,
    evaluate = identical(evaluation_state, "feature_evaluation_complete"),
    combine = nzchar(combined) && file.exists(combined) && file.info(combined)$size > 0,
    combined_summary = nzchar(feature) && file.exists(feature) && file.info(feature)$size > 0,
    classifier = nzchar(prediction) && file.exists(prediction) && file.info(prediction)$size > 0,
    results = length(list.files(file.path(root, "classifier"), pattern = "^essentiality_predictions\\..*\\.csv$",
                                recursive = TRUE, full.names = TRUE)) > 0L
  )
  states <- c(
    normalize = if (!length(parents)) "blocked" else if (stage_files$normalize) "complete" else "ready",
    evaluate = if (!stage_files$normalize) "blocked" else if (stage_files$evaluate) "complete" else "ready",
    combine = if (!nzchar(target_tag) || !nrow(normalized) || any(!normalized$Available)) "blocked" else if (stage_files$combine) "complete" else "ready",
    combined_summary = if (!stage_files$combine) "blocked" else if (stage_files$combined_summary) "complete" else "ready",
    classifier = if (!stage_files$combined_summary) "blocked" else if (stage_files$classifier) "complete" else "ready",
    results = if (stage_files$results) "complete" else "blocked"
  )
  for (stage in c("normalize", "combine", "combined_summary", "classifier")) {
    manifest <- if (nzchar(target_tag)) essentiality_manifest(root, stage, target_tag) else NULL
    if (!is.null(manifest) && identical(manifest$status %||% "", "failed") &&
        !isTRUE(stage_files[[stage]])) states[[stage]] <- "failed"
    if (!is.null(manifest) && identical(manifest$status %||% "", "skipped") && isTRUE(stage_files[[stage]])) states[[stage]] <- "cached"
  }
  if (!is.null(job_progress)) {
    map <- c(sample_normalization = "normalize", summary_normalized = "evaluate",
             combined_hits = "combine", summary_combined = "combined_summary",
             classifier = "classifier")
    progress_stage <- as.character(job_progress$stage %||% "")
    current <- if (progress_stage %in% names(map)) unname(map[[progress_stage]]) else ""
    progress_status <- as.character(job_progress$status %||% "")
    if (length(current) && nzchar(current)) {
      if (progress_status %in% c("queued", "starting", "running")) states[[current]] <- "running"
      else if (identical(progress_status, "failed") && !isTRUE(stage_files[[current]])) states[[current]] <- "failed"
      else if (progress_status %in% c("cached", "cached_success") && isTRUE(stage_files[[current]])) states[[current]] <- "cached"
    }
  }
  attr(states, "evaluation_state") <- evaluation_state
  states
}

essentiality_next_stage <- function(states) {
  order <- c("normalize", "evaluate", "combine", "combined_summary", "classifier", "results")
  labels <- c(normalize = "normalize", evaluate = "evaluate", combine = "combine",
              combined_summary = "combined_summary", classifier = "classifier", results = "results")
  for (stage in order) if (states[[stage]] %in% c("ready", "failed")) return(labels[[stage]])
  for (stage in order) if (states[[stage]] == "blocked") {
    previous <- match(stage, order) - 1L
    if (previous >= 1L) return(labels[[order[[previous]]]])
  }
  "results"
}
