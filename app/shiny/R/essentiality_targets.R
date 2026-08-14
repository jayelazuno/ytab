essentiality_numeric_text <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || is.na(value) || !is.finite(value) || value <= 0) return("")
  format(value, scientific = FALSE, trim = TRUE, digits = 15)
}

essentiality_target_tag <- function(value) {
  text <- essentiality_numeric_text(value)
  if (!nzchar(text)) stop("Target must be a finite positive number.", call. = FALSE)
  pieces <- strsplit(text, ".", fixed = TRUE)[[1]]
  whole <- pieces[[1]]
  fraction <- if (length(pieces) > 1L) sub("0+$", "", pieces[[2]]) else ""
  paste0("T", strrep("0", max(0L, 3L - nchar(whole))), whole,
         if (nzchar(fraction)) paste0("p", fraction) else "")
}

essentiality_tag_value <- function(tag) {
  tag <- trimws(as.character(tag %||% ""))
  if (!grepl("^T[0-9]+(?:p[0-9]+)?$", tag)) return(NA_real_)
  value <- suppressWarnings(as.numeric(sub("p", ".", substring(tag, 2L), fixed = TRUE)))
  if (is.na(value) || !is.finite(value) || value <= 0) return(NA_real_)
  if (!identical(essentiality_target_tag(value), tag)) return(NA_real_)
  value
}

parse_essentiality_targets <- function(text) {
  raw <- as.character(text %||% "")
  tokens <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  if (!nzchar(trimws(raw))) tokens <- ""
  seen <- character()
  rows <- lapply(tokens, function(token) {
    value <- suppressWarnings(as.numeric(token))
    valid <- nzchar(token) && !is.na(value) && is.finite(value) && value > 0
    tag <- if (valid) essentiality_target_tag(value) else ""
    reason <- if (!nzchar(token)) "Empty target" else if (!valid) "Enter a finite positive number" else if (tag %in% seen) "Duplicate target" else ""
    if (valid && !nzchar(reason)) seen <<- c(seen, tag)
    data.frame(Token = token, Value = if (valid) value else NA_real_, Tag = tag,
               Valid = if (nzchar(reason)) "No" else "Yes", Reason = reason,
               stringsAsFactors = FALSE, check.names = FALSE)
  })
  table <- do.call(rbind, rows)
  valid <- table$Valid == "Yes"
  list(
    table = table,
    valid = length(tokens) > 0L && all(valid),
    values = unname(table$Value[valid]),
    tags = unname(table$Tag[valid]),
    argument = paste(essentiality_numeric_text(table$Value[valid]), collapse = ","),
    errors = unname(table$Reason[!valid])
  )
}

essentiality_read_json <- function(path) {
  if (!file.exists(path) || !requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
}

essentiality_valid_tag <- function(tag) !is.na(essentiality_tag_value(tag))

format_target_label <- function(value = NA_real_, tag = "",
                                mode = c("user", "technical", "badge")) {
  mode <- match.arg(mode)
  target_tag <- as.character(tag %||% "")
  target_value <- essentiality_numeric_text(value)
  if (!essentiality_valid_tag(target_tag) && nzchar(target_value))
    target_tag <- essentiality_target_tag(value)
  if (!nzchar(target_value) && essentiality_valid_tag(target_tag))
    target_value <- essentiality_numeric_text(essentiality_tag_value(target_tag))
  if (!essentiality_valid_tag(target_tag)) return("Not available")
  switch(
    mode,
    user = paste("Normalization target:", target_value),
    technical = paste0("Target value: ", target_value,
                       "; technical tag: ", target_tag),
    badge = target_tag
  )
}

format_essentiality_target <- function(target_value = NA_real_, target_tag = "",
                                       context = c("badge", "value", "recommendation", "user", "technical")) {
  context <- match.arg(context)
  if (context %in% c("user", "technical", "badge"))
    return(format_target_label(target_value, target_tag, context))
  tag <- as.character(target_tag %||% "")
  value <- essentiality_numeric_text(target_value)
  if (!essentiality_valid_tag(tag) && nzchar(value)) tag <- essentiality_target_tag(target_value)
  if (!nzchar(value) && essentiality_valid_tag(tag)) value <- essentiality_numeric_text(essentiality_tag_value(tag))
  if (!essentiality_valid_tag(tag)) return("Not available")
  switch(context, badge = tag, value = value,
         recommendation = paste("Recommended", format_target_label(value, tag, "user")))
}

essentiality_target_details <- function(target_value = NA_real_, target_tag = "") {
  list(
    target_value = format_essentiality_target(target_value, target_tag, "value"),
    target_tag = format_essentiality_target(target_value, target_tag, "badge")
  )
}

essentiality_feature_evaluation_state <- function(project_root) {
  base <- file.path(project_root, "sample_normalization")
  evaluation <- file.path(base, "normalization_target_evaluation.csv")
  data <- if (file.exists(evaluation)) tryCatch(
    read.csv(evaluation, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) data.frame()
  ) else data.frame()
  target_column <- intersect(c("target", "target_value", "target_tag"), names(data))
  feature_column <- intersect(c("feature_retention", "min_feature_retention",
                                "min_feature_retention_fraction",
                                "passed_feature_threshold"), names(data))
  complete <- nrow(data) > 0L && length(target_column) > 0L && length(feature_column) > 0L
  normalized <- dir.exists(base) && length(list.files(
    base, pattern = "_normalized_hits\\.txt$", recursive = TRUE, full.names = TRUE
  )) > 0L
  recommendation <- essentiality_read_json(file.path(base, "normalization_recommendation.json"))
  if (complete) "feature_evaluation_complete" else if (normalized)
    "feature_evaluation_ready" else if (!is.null(recommendation))
      "preliminary_recommendation_only" else "blocked"
}

essentiality_named_list <- function(value) {
  if (is.null(value)) return(list())
  if (is.data.frame(value)) {
    if (!nrow(value)) return(list())
    return(as.list(value[1L, , drop = FALSE]))
  }
  if (is.list(value)) return(value)
  if (is.atomic(value) && length(names(value)) && any(nzchar(names(value))))
    return(as.list(value))
  list()
}

essentiality_first_value <- function(data, names, default = NULL) {
  data <- essentiality_named_list(data)
  for (name in names) {
    value <- data[[name]]
    if (!is.null(value) && length(value)) return(value[[1L]])
  }
  default
}

normalize_recommendation_state <- function(value, type = "", label = "", path = "",
                                           raw_data = NULL) {
  data <- essentiality_named_list(value)
  raw <- essentiality_named_list(raw_data %||% value)
  target <- suppressWarnings(as.numeric(essentiality_first_value(
    data, c("target", "recommended_target", "value"), NA_real_
  )))
  tag <- trimws(as.character(essentiality_first_value(
    data, c("target_tag", "recommended_target_tag", "tag"), ""
  )))
  if (!essentiality_valid_tag(tag) && length(target) == 1L &&
      is.finite(target) && target > 0)
    tag <- essentiality_target_tag(target)
  if ((!length(target) || !is.finite(target) || target <= 0) &&
      essentiality_valid_tag(tag))
    target <- essentiality_tag_value(tag)
  resolved_type <- trimws(as.character(essentiality_first_value(
    data, c("type", "recommendation_level"), type
  )))
  recommendation_type <- trimws(as.character(essentiality_first_value(
    data, c("recommendation_type", "label"), label
  )))
  if (!nzchar(recommendation_type)) {
    recommendation_type <- switch(
      resolved_type,
      feature = "Feature-evaluated recommendation",
      site = "Preliminary site-retention recommendation",
      "No recommendation"
    )
  }
  valid <- length(target) == 1L && is.finite(target) && target > 0 &&
    essentiality_valid_tag(tag) &&
    isTRUE(all.equal(as.numeric(target), essentiality_tag_value(tag)))
  site <- suppressWarnings(as.numeric(essentiality_first_value(
    data, c("site_retention", "min_site_retention_observed",
            "min_hit_site_retention_fraction"), NA_real_
  )))
  feature <- suppressWarnings(as.numeric(essentiality_first_value(
    data, c("feature_retention", "min_feature_retention",
            "min_feature_retention_fraction"), NA_real_
  )))
  parents <- suppressWarnings(as.integer(essentiality_first_value(
    data, c("parents_passing", "parent_libraries_passing", "sample_count"), NA_integer_
  )))
  reason <- as.character(essentiality_first_value(
    data, c("reason", "message"),
    if (valid) "" else "No valid target recommendation is available."
  ))
  list(
    available = valid,
    target = if (valid) as.numeric(target) else NA_real_,
    target_tag = if (valid) tag else "",
    recommendation_type = if (valid) recommendation_type else "No recommendation",
    site_retention = if (length(site)) site[[1L]] else NA_real_,
    feature_retention = if (length(feature)) feature[[1L]] else NA_real_,
    parents_passing = if (length(parents)) parents[[1L]] else NA_integer_,
    reason = if (length(reason)) reason[[1L]] else "",
    value = if (valid) as.numeric(target) else NA_real_,
    tag = if (valid) tag else "",
    type = if (valid) resolved_type else "none",
    label = if (valid) recommendation_type else "No recommendation",
    path = as.character(path %||% ""),
    data = raw,
    preliminary = valid && identical(resolved_type, "site")
  )
}

discover_essentiality_targets <- function(project_root) {
  stages <- c("sample_normalization", "summary_normalized", "combined_hits",
              "summary_combined", "classifier")
  found <- list()
  add <- function(tag, source) {
    tag <- as.character(tag %||% "")
    if (!essentiality_valid_tag(tag)) return()
    found[[tag]] <<- unique(c(found[[tag]] %||% character(), source))
  }
  for (stage in stages) {
    root <- file.path(project_root, stage)
    if (!dir.exists(root)) next
    paths <- c(list.dirs(root, recursive = TRUE, full.names = FALSE),
               list.files(root, recursive = TRUE, full.names = FALSE))
    hits <- unique(unlist(regmatches(paths, gregexpr("T[0-9]+(?:p[0-9]+)?", paths, perl = TRUE))))
    for (tag in hits[nzchar(hits)]) add(tag, stage)
  }
  normalization <- file.path(project_root, "sample_normalization")
  for (name in c("normalization_feature_recommendation.json", "normalization_recommendation.json")) {
    data <- essentiality_read_json(file.path(normalization, name))
    normalized <- normalize_recommendation_state(data)
    add(normalized$target_tag, name)
  }
  for (name in c("normalization_comparison.csv", "normalization_target_evaluation.csv",
                 "normalization_target_summary.csv")) {
    path <- file.path(normalization, name)
    if (!file.exists(path)) next
    data <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
    if ("target_tag" %in% names(data)) for (tag in unique(as.character(data$target_tag))) add(tag, name)
  }
  if (!length(found)) return(data.frame(Value = numeric(), Tag = character(), Sources = character(),
                                        stringsAsFactors = FALSE, check.names = FALSE))
  result <- data.frame(
    Value = vapply(names(found), essentiality_tag_value, numeric(1)),
    Tag = names(found),
    Sources = vapply(found, function(x) paste(sort(x), collapse = "; "), character(1)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  result[order(result$Value, result$Tag), , drop = FALSE]
}

essentiality_recommendation <- function(project_root) {
  base <- file.path(project_root, "sample_normalization")
  candidates <- list(
    list(file = "normalization_feature_recommendation.json",
         type = "feature", label = "Feature-evaluated recommendation"),
    list(file = "normalization_recommendation.json",
         type = "site", label = "Preliminary site-retention recommendation")
  )
  if (!identical(essentiality_feature_evaluation_state(project_root),
                 "feature_evaluation_complete")) candidates <- candidates[2L]
  for (candidate in candidates) {
    path <- file.path(base, candidate$file)
    data <- essentiality_read_json(path)
    normalized <- normalize_recommendation_state(
      data, type = candidate$type, label = candidate$label,
      path = path, raw_data = data
    )
    if (isTRUE(normalized$available)) return(normalized)
  }
  normalize_recommendation_state(NULL)
}

essentiality_target_choices <- function(targets, suffix = "") {
  if (is.null(targets) || !nrow(targets)) return(list())
  labels <- paste0(vapply(seq_len(nrow(targets)), function(index)
    format_target_label(targets$Value[[index]], targets$Tag[[index]], "user"),
    character(1)), suffix)
  as.list(setNames(unname(targets$Tag), unname(labels)))
}

essentiality_combined_feature_targets <- function(project_root) {
  all <- discover_essentiality_targets(project_root)
  if (!nrow(all)) return(all)
  keep <- vapply(all$Tag, function(tag) {
    path <- file.path(project_root, "summary_combined", tag,
                      paste0("combined_feature_table.", tag, ".txt"))
    file.exists(path) && !dir.exists(path) && file.info(path)$size > 0
  }, logical(1))
  all[keep, , drop = FALSE]
}
