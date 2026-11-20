# Parent-library detection for the Essentiality Run tab.
#
# The rest of this file served the retired classifier's stages (final target, stage state, manifests,
# normalized inputs) and went with that subsystem to archive/classifier_subsystem_2026-11-20/. Two of its
# functions called helpers that were removed with it, which would have failed the moment those paths ran.

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
    reason <- if (!included[[index]]) "Excluded by dataset sample sheet" else if (explicit_treated) "Treated library excluded" else if (flag_parent) "Parent flag" else if (metadata_parent) "Parent metadata" else if (fallback_parent) "Parent name fallback" else "Parent role not identified"
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

