essentiality_prediction_columns <- function(data) {
  names_lower <- tolower(names(data))
  first <- function(patterns) {
    for (pattern in patterns) {
      hit <- which(grepl(pattern, names_lower, perl = TRUE))
      if (length(hit)) return(names(data)[hit[[1]]])
    }
    ""
  }
  list(
    feature = first(c("^standard name$", "^feature.?id$", "^gene.?id$", "^id$")),
    gene = first(c("^common name$", "^gene$", "^standard.?name$")),
    label = first(c("ess\\. for fpr", "verdict$", "classifier.?label", "prediction.?label")),
    score = first(c("^rf - g4$", "probab", "prediction.?score", "^score$")),
    ortholog = first(c("^sc ortholog$", "s\\.? ?cerevisiae.*ortholog", "ortholog")),
    exclusion = first(c("exclusion.?reason", "excluded")),
    metric = first(c("^hits$", "^reads$", "insertion.?index", "freedom.?index", "neighborhood.?index"))
  )
}

essentiality_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE) || !file.exists(path) || dir.exists(path)) return("")
  tryCatch(digest::digest(file = path, algo = "sha256", serialize = FALSE),
           error = function(e) "")
}

essentiality_classifier_manifest_integrity <- function(manifest) {
  if (is.null(manifest) || !identical(as.character(manifest$status %||% ""), "success"))
    return(FALSE)
  required <- c("input_combined_feature_table", "input_sha256",
                "classifier_script", "classifier_script_sha256",
                "classifier_runner", "classifier_runner_sha256",
                "scientific_parameters", "classifier_resources", "seed",
                "species", "target", "target_tag")
  if (!all(required %in% names(manifest))) return(FALSE)
  if (!identical(essentiality_sha256(manifest$input_combined_feature_table),
                 as.character(manifest$input_sha256))) return(FALSE)
  if (!identical(essentiality_sha256(manifest$classifier_script),
                 as.character(manifest$classifier_script_sha256))) return(FALSE)
  if (!identical(essentiality_sha256(manifest$classifier_runner),
                 as.character(manifest$classifier_runner_sha256))) return(FALSE)
  resources <- manifest$classifier_resources
  if (!length(resources)) return(FALSE)
  for (entries in resources) {
    if (!length(entries)) return(FALSE)
    for (entry in entries) {
      if (!identical(essentiality_sha256(entry$path %||% ""),
                     as.character(entry$sha256 %||% ""))) return(FALSE)
    }
  }
  TRUE
}

essentiality_classifier_cache_matches <- function(project_root, tag, seed) {
  manifest <- essentiality_manifest(project_root, "classifier", tag)
  stable <- file.path(project_root, "classifier", tag,
                      paste0("essentiality_predictions.", tag, ".csv"))
  feature <- file.path(project_root, "summary_combined", tag,
                       paste0("combined_feature_table.", tag, ".txt"))
  essentiality_classifier_manifest_integrity(manifest) &&
    file.exists(stable) && file.info(stable)$size > 0 &&
    identical(as.character(manifest$target_tag), tag) &&
    isTRUE(all.equal(as.numeric(manifest$target), essentiality_tag_value(tag))) &&
    identical(as.character(manifest$input_combined_feature_table),
              normalizePath(feature, winslash = "/", mustWork = FALSE)) &&
    identical(as.integer(manifest$seed), as.integer(seed))
}

read_essentiality_predictions <- function(path) {
  if (!nzchar(path %||% "") || !file.exists(path) || dir.exists(path)) return(data.frame())
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) data.frame())
}

discover_classifier_results <- function(project_root, selected_tag = "", final_tag = "") {
  root <- file.path(project_root, "classifier")
  if (!dir.exists(root)) return(list())
  stable <- list.files(root, pattern = "^essentiality_predictions\\.T[0-9]+(?:p[0-9]+)?\\.csv$",
                       recursive = TRUE, full.names = TRUE)
  legacy <- list.files(root, pattern = "\\.(predictions\\.csv|tsv)$",
                       recursive = TRUE, full.names = TRUE)
  legacy <- setdiff(legacy, stable)
  build <- function(path, legacy_result = FALSE) {
    parent_tag <- basename(dirname(path))
    file_tag <- regmatches(basename(path), regexpr("T[0-9]+(?:p[0-9]+)?", basename(path), perl = TRUE))
    tag <- if (essentiality_valid_tag(parent_tag)) parent_tag else if (length(file_tag) && essentiality_valid_tag(file_tag)) file_tag else ""
    if (!nzchar(tag)) return(NULL)
    data <- read_essentiality_predictions(path)
    if (!nrow(data)) return(NULL)
    manifest_path <- file.path(project_root, "manifests", "classifier",
                               paste0(tag, ".classifier_manifest.json"))
    manifest <- essentiality_read_json(manifest_path)
    status <- as.character(manifest$status %||% if (legacy_result) "legacy" else "historical")
    valid_manifest <- essentiality_classifier_manifest_integrity(manifest)
    classification <- if (identical(tag, final_tag)) "final-target result" else if (identical(tag, selected_tag) && valid_manifest) "matching classifier result" else if (legacy_result || is.null(manifest)) "legacy result" else if (identical(status, "failed")) "failed result" else "historical classifier result"
    metadata_path <- file.path(dirname(path), paste0("classifier_run_metadata.", tag, ".json"))
    list(
      tag = tag, value = essentiality_tag_value(tag), path = path, manifest_path = manifest_path,
      manifest = manifest, metadata_path = metadata_path,
      metadata = essentiality_read_json(metadata_path), legacy = legacy_result || is.null(manifest),
      status = status, classification = classification, rows = nrow(data),
      modified = file.info(path)$mtime
    )
  }
  results <- Filter(Negate(is.null), c(lapply(stable, build), lapply(legacy, build, legacy_result = TRUE)))
  if (!length(results)) return(list())
  # Prefer the stable result when a legacy table exists for the same target.
  groups <- split(results, vapply(results, `[[`, "", "tag"))
  results <- lapply(groups, function(items) {
    stable_items <- Filter(function(x) !x$legacy, items)
    candidates <- if (length(stable_items)) stable_items else items
    candidates[[which.max(vapply(candidates, function(x) as.numeric(x$modified), numeric(1)))]]
  })
  results[order(vapply(results, `[[`, numeric(1), "value"))]
}

essentiality_result_labels <- function(results, recommendation_tag = "", final_tag = "") {
  if (!length(results)) return(list())
  values <- vapply(results, `[[`, "", "tag")
  labels <- vapply(results, function(result) {
    suffix <- c(if (identical(result$tag, recommendation_tag)) "recommended",
                if (identical(result$tag, final_tag)) "final target",
                if (result$legacy) "legacy")
    paste(result$tag, if (length(suffix)) paste0("— ", paste(suffix, collapse = ", ")) else "")
  }, character(1))
  as.list(setNames(unname(values), unname(labels)))
}

choose_essentiality_result <- function(results, selected = "", selected_target = "", final_target = "") {
  if (!length(results)) return(NULL)
  for (tag in c(selected, selected_target, final_target)) {
    hit <- Filter(function(x) identical(x$tag, tag), results)
    if (length(hit)) return(hit[[1]])
  }
  results[[length(results)]]
}

essentiality_result_summary <- function(result) {
  if (is.null(result)) return(list(total = 0L, counts = integer(), excluded = 0L,
                                   label_column = "", score_column = ""))
  data <- read_essentiality_predictions(result$path)
  columns <- essentiality_prediction_columns(data)
  counts <- if (nzchar(columns$label)) table(as.character(data[[columns$label]]), useNA = "ifany") else integer()
  excluded <- if (nzchar(columns$exclusion)) {
    values <- trimws(as.character(data[[columns$exclusion]]))
    sum(!is.na(values) & nzchar(values))
  } else {
    as.integer(result$manifest$excluded_feature_count %||%
                 max(0, as.integer(result$manifest$input_feature_count %||% nrow(data)) - nrow(data)))
  }
  list(total = nrow(data), counts = counts, excluded = excluded,
       label_column = columns$label, score_column = columns$score, columns = columns)
}

filter_essentiality_results <- function(data, search = "", label = "All", inclusion = "All") {
  if (!nrow(data)) return(data)
  columns <- essentiality_prediction_columns(data)
  keep <- rep(TRUE, nrow(data))
  query <- trimws(tolower(as.character(search %||% "")))
  if (nzchar(query)) {
    searchable <- unique(c(columns$feature, columns$gene, columns$ortholog))
    searchable <- searchable[nzchar(searchable)]
    if (length(searchable)) keep <- keep & apply(data[, searchable, drop = FALSE], 1L, function(row) {
      any(grepl(query, tolower(as.character(row)), fixed = TRUE))
    })
  }
  if (!identical(label, "All") && nzchar(columns$label))
    keep <- keep & as.character(data[[columns$label]]) == label
  if (!identical(inclusion, "All")) {
    excluded <- if (nzchar(columns$exclusion)) nzchar(trimws(as.character(data[[columns$exclusion]]))) else rep(FALSE, nrow(data))
    keep <- keep & if (identical(inclusion, "Excluded")) excluded else !excluded
  }
  data[keep, , drop = FALSE]
}

essentiality_visible_results <- function(data) {
  if (!nrow(data)) return(data.frame())
  columns <- essentiality_prediction_columns(data)
  chosen <- unique(c(columns$feature, columns$gene, columns$label, columns$score,
                     columns$ortholog, columns$metric, columns$exclusion))
  chosen <- chosen[nzchar(chosen) & chosen %in% names(data)]
  if (!length(chosen)) chosen <- names(data)[seq_len(min(8L, ncol(data)))]
  result <- data[, chosen, drop = FALSE]
  labels <- c(
    setNames("Feature ID", columns$feature), setNames("Gene", columns$gene),
    setNames("Classifier label", columns$label), setNames("Prediction score", columns$score),
    setNames("S. cerevisiae ortholog", columns$ortholog),
    setNames("Input feature metric", columns$metric),
    setNames("Exclusion reason", columns$exclusion)
  )
  labels <- labels[nzchar(names(labels))]
  names(result) <- ifelse(names(result) %in% names(labels), labels[names(result)], names(result))
  result
}

essentiality_target_summary_data <- function(project_root) {
  path <- file.path(project_root, "sample_normalization", "normalization_target_summary.csv")
  if (!file.exists(path)) return(data.frame())
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) data.frame())
}

essentiality_target_evaluation_data <- function(project_root) {
  path <- file.path(project_root, "sample_normalization", "normalization_target_evaluation.csv")
  if (!file.exists(path)) return(data.frame())
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) data.frame())
}

essentiality_target_summary_view <- function(data) {
  if (!nrow(data)) return(data.frame())
  aliases <- list(
    `Normalization target` = c("target"), `Technical tag` = c("target_tag"),
    `Parent libraries evaluated` = c("sample_count", "parent_libraries_evaluated"),
    `Minimum site retention` = c("min_hit_site_retention_fraction", "minimum_site_retention"),
    `Median site retention` = c("median_hit_site_retention_fraction", "median_site_retention"),
    `Minimum feature retention` = c("min_feature_retention_fraction", "minimum_feature_retention"),
    `Median feature retention` = c("median_feature_retention_fraction", "median_feature_retention"),
    `Passed site threshold` = c("passed_site_threshold"),
    `Passed feature threshold` = c("passed_feature_threshold"),
    Recommended = c("recommended_by_feature_retention", "recommended")
  )
  selected <- list()
  for (label in names(aliases)) {
    hit <- intersect(aliases[[label]], names(data))
    if (length(hit)) selected[[label]] <- data[[hit[[1]]]]
  }
  as.data.frame(selected, stringsAsFactors = FALSE, check.names = FALSE)
}

essentiality_target_visual_state <- function(data, recommendation = NULL) {
  if (is.null(data) || !is.data.frame(data)) data <- data.frame()
  count <- nrow(data)
  mode <- if (count == 0L) "empty" else if (count == 1L) "single" else
    if (count == 2L) "compact" else "trends"
  list(
    target_count = count,
    mode = mode,
    show_table = count >= 2L,
    show_site_plot = count >= 3L,
    show_feature_plot = count >= 3L,
    show_combined_plot = count == 2L || count >= 3L,
    recommendation = normalize_recommendation_state(recommendation),
    message = switch(
      mode,
      empty = "Run parent normalization before evaluating targets.",
      single = paste(
        "Only one target is available; a cross-target trend cannot be displayed.",
        "Additional targets are needed to compare retention trends."
      ),
      compact = "Two targets are available for a compact comparison.",
      trends = "Retention trends are available across three or more targets."
    )
  )
}

essentiality_results_view_contract <- function(view = "overview") {
  allowed <- c("overview", "predictions", "visualizations", "provenance", "downloads")
  view <- tolower(as.character(view %||% "overview")[[1L]])
  if (!view %in% allowed) view <- "overview"
  components <- switch(
    view,
    overview = c("summary", "label_distribution", "final_target_action"),
    predictions = c("search", "label_filter", "inclusion_filter", "predictions_table"),
    visualizations = c("label_distribution", "score_distribution", "feature_metric_plot"),
    provenance = c("target_provenance", "technical_provenance"),
    downloads = c("download_actions")
  )
  list(
    view = view,
    components = components,
    renders_predictions_table = identical(view, "predictions"),
    renders_visualizations = identical(view, "visualizations"),
    renders_provenance = identical(view, "provenance"),
    renders_downloads = identical(view, "downloads")
  )
}

essentiality_download_availability <- function(paths) {
  if (is.null(paths) || !length(paths)) return(logical())
  vapply(paths, function(path) {
    path <- as.character(path %||% "")
    length(path) == 1L && nzchar(path) && file.exists(path) && !dir.exists(path)
  }, logical(1))
}

essentiality_result_provenance <- function(project_root, result, parents) {
  if (is.null(result)) return(list())
  tag <- result$tag
  target_summary <- essentiality_target_summary_data(project_root)
  target_row <- if (nrow(target_summary) && "target_tag" %in% names(target_summary))
    target_summary[target_summary$target_tag == tag, , drop = FALSE] else data.frame()
  combine <- essentiality_manifest(project_root, "combine", tag)
  combined <- essentiality_manifest(project_root, "combined_summary", tag)
  manifest <- result$manifest %||% list()
  recommendation <- essentiality_recommendation(project_root)
  list(
    target_value = essentiality_tag_value(tag), target_tag = tag,
    recommendation_type = if (identical(recommendation$tag, tag)) recommendation$label else "Selected evaluated target",
    minimum_site_retention = if (nrow(target_row) && "min_hit_site_retention_fraction" %in% names(target_row)) target_row$min_hit_site_retention_fraction[[1]] else NA,
    minimum_feature_retention = if (nrow(target_row) && "min_feature_retention_fraction" %in% names(target_row)) target_row$min_feature_retention_fraction[[1]] else NA,
    parent_libraries = length(parents),
    combined_sites = combine$total_combined_sites %||% NA,
    combined_reads = combine$total_combined_reads %||% NA,
    combined_feature_count = manifest$input_feature_count %||% NA,
    random_seed = manifest$seed %||% result$metadata$seed %||% 0,
    resources = manifest$classifier_resources %||% result$metadata$classifier_resources %||% list(),
    final = identical(essentiality_final_target(project_root), tag)
  )
}

essentiality_download_paths <- function(project_root, result) {
  tag <- result$tag %||% ""
  list(
    predictions = result$path %||% "",
    combined_feature = if (nzchar(tag)) file.path(project_root, "summary_combined", tag, paste0("combined_feature_table.", tag, ".txt")) else "",
    target_evaluation = file.path(project_root, "sample_normalization", "normalization_target_evaluation.csv"),
    normalization_recommendation = essentiality_recommendation(project_root)$path,
    classifier_manifest = result$manifest_path %||% "",
    final_target = file.path(project_root, "config", "final_classifier_target.txt")
  )
}
