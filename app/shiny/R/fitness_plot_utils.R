fitness_ma_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
fitness_ma_lfc_support_threshold <- function() 0
fitness_ma_candidate_lfc_threshold <- function() 1
fitness_ma_candidate_cpm_threshold <- function() 1

fitness_ma_first_nonempty <- function(...) {
  values <- list(...)
  if (!length(values)) return(character())
  n <- max(vapply(values, length, integer(1)))
  out <- rep("", n)
  for (value in values) {
    value <- as.character(value %||% "")
    if (length(value) == 1L && n > 1L) value <- rep(value, n)
    value[is.na(value)] <- ""
    take <- !nzchar(out) & nzchar(trimws(value))
    out[take] <- value[take]
  }
  out
}

fitness_ma_apply_display_labels <- function(data, repo_root = NULL) {
  if (!is.data.frame(data) || !nrow(data) || !"feature_id" %in% names(data)) return(data)
  if (!is.null(repo_root) && exists("ytab_join_glabrata_display", mode = "function")) {
    data <- tryCatch(ytab_join_glabrata_display(data, repo_root, "feature_id"),
                     error = function(e) data)
  }
  data$label <- fitness_ma_first_nonempty(
    if ("gene_display_name" %in% names(data)) data$gene_display_name else "",
    if ("cagl_display_id" %in% names(data)) data$cagl_display_id else "",
    if ("label" %in% names(data)) data$label else "",
    data$feature_id
  )
  data
}

fitness_ma_pair_choices <- function(result) {
  if (is.null(result) || !nzchar(result$output_dir %||% "")) return(setNames(character(), character()))
  path <- file.path(result$output_dir, "tables", "treated_vs_parent.by_pool.log2fc_z.csv")
  if (!file.exists(path)) return(setNames(character(), character()))
  x <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (!nrow(x) || !"contrast" %in% names(x)) return(setNames(character(), character()))
  ids <- unique(as.character(x$contrast)); ids <- ids[nzchar(ids) & !is.na(ids)]
  labels <- vapply(ids, function(id) {
    row <- x[match(id, x$contrast), , drop = FALSE]
    pool <- if ("pool" %in% names(row)) paste0("pool", row$pool[[1]]) else id
    parent <- if ("parent_sample" %in% names(row)) row$parent_sample[[1]] else "control"
    treated <- if ("treated_sample" %in% names(row)) row$treated_sample[[1]] else "treated"
    paste(pool, parent, "vs", treated)
  }, "")
  setNames(ids, labels)
}

fitness_ma_available_pools <- function(result) {
  if (is.null(result) || !nzchar(result$output_dir %||% "")) return(character())
  path <- file.path(result$output_dir, "tables", "treated_vs_parent.by_pool.log2fc_z.csv")
  if (!file.exists(path)) return(character())
  x <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (!nrow(x) || !"contrast" %in% names(x)) return(character())
  keep <- !is.na(x$contrast) & nzchar(as.character(x$contrast))
  if (all(c("log2FC", "z") %in% names(x)))
    keep <- keep & is.finite(fitness_ma_numeric(x$log2FC)) & is.finite(fitness_ma_numeric(x$z))
  sort(unique(as.character(x$contrast[keep])))
}

fitness_ma_data <- function(result, mode = "combined", pair = "") {
  empty <- data.frame(feature_id = character(), label = character(), mean_abundance = numeric(),
                      log2fc = numeric(), call = character(), stringsAsFactors = FALSE)
  if (is.null(result) || !nzchar(result$output_dir %||% "")) return(empty)
  root <- file.path(result$output_dir, "tables")
  summary_path <- file.path(root, "treated_vs_parent.summary_by_feature.csv")
  pool_path <- file.path(root, "treated_vs_parent.by_pool.log2fc_z.csv")
  cpm_path <- file.path(root, "feature_reads_cpm.long.csv")
  if (!file.exists(summary_path)) return(empty)
  summary <- tryCatch(read.csv(summary_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (!nrow(summary) || !all(c("feature_id", "mean_log2FC") %in% names(summary))) return(empty)
  id <- as.character(summary$feature_id)
  labels <- if ("common_name" %in% names(summary)) as.character(summary$common_name) else rep("", nrow(summary))
  if ("standard_name" %in% names(summary)) labels[is.na(labels) | !nzchar(trimws(labels))] <- as.character(summary$standard_name[is.na(labels) | !nzchar(trimws(labels))])
  labels[is.na(labels) | !nzchar(trimws(labels))] <- id[is.na(labels) | !nzchar(trimws(labels))]
  if (identical(mode, "individual")) {
    if (!file.exists(pool_path)) return(empty)
    pool <- tryCatch(read.csv(pool_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
    if (!nrow(pool) || !"contrast" %in% names(pool)) return(empty)
    if (!nzchar(pair)) pair <- fitness_ma_pair_choices(result)[[1]] %||% ""
    pool <- pool[as.character(pool$contrast) == pair, , drop = FALSE]
    if (!nrow(pool) || !all(c("feature_id", "log2FC") %in% names(pool))) return(empty)
    out <- pool[, intersect(c("feature_id", "log2FC", "call", "standard_name", "common_name", "parent_cpm", "treated_cpm"), names(pool)), drop = FALSE]
    out$feature_id <- as.character(out$feature_id)
    out$log2fc <- fitness_ma_numeric(out$log2FC)
    out$mean_abundance <- if (all(c("parent_cpm", "treated_cpm") %in% names(out))) rowMeans(cbind(fitness_ma_numeric(out$parent_cpm), fitness_ma_numeric(out$treated_cpm)), na.rm = TRUE) else NA_real_
    out$call <- if ("call" %in% names(out)) as.character(out$call) else "none"
    out$label <- if ("common_name" %in% names(out)) as.character(out$common_name) else rep("", nrow(out))
    if ("standard_name" %in% names(out)) out$label[is.na(out$label) | !nzchar(trimws(out$label))] <- as.character(out$standard_name[is.na(out$label) | !nzchar(trimws(out$label))])
    out$label[is.na(out$label) | !nzchar(trimws(out$label))] <- out$feature_id[is.na(out$label) | !nzchar(trimws(out$label))]
    return(out[, c("feature_id", "label", "mean_abundance", "log2fc", "call"), drop = FALSE])
  }
  if (!file.exists(cpm_path)) return(empty)
  cpm <- tryCatch(read.csv(cpm_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (!nrow(cpm) || !all(c("feature_id", "cpm") %in% names(cpm))) return(empty)
  abundance <- aggregate(fitness_ma_numeric(cpm$cpm), list(feature_id = as.character(cpm$feature_id)), mean, na.rm = TRUE)
  names(abundance)[2] <- "mean_abundance"
  out <- data.frame(feature_id = id, label = labels, mean_abundance = abundance$mean_abundance[match(id, abundance$feature_id)],
                    log2fc = fitness_ma_numeric(summary$mean_log2FC),
                    call = if ("final_call" %in% names(summary)) as.character(summary$final_call) else "none",
                    stringsAsFactors = FALSE)
  out
}

fitness_ma_rank_data <- function(result, mode = "combined", pair = "", direction = "both", min_support = 0L) {
  data <- fitness_ma_data(result, mode, pair)
  if (!nrow(data)) return(data)
  min_support <- suppressWarnings(as.integer(min_support %||% 0L))
  if (is.na(min_support) || min_support < 0L) min_support <- 0L
  data$depleted_support_n <- 0L; data$enriched_support_n <- 0L; data$depleted_stored_call_n <- 0L; data$enriched_stored_call_n <- 0L; data$valid_pool_n <- 1L
  data$rank_z_strength <- 0; data$rank_lfc_strength <- abs(fitness_ma_numeric(data$log2fc)); data$rank_cpm_support <- fitness_ma_numeric(data$mean_abundance)
  data$rank_support_n <- 0L; data$ranking_basis <- "ranked by directional log2FC magnitude, CPM/read support, local z-score support, and stable feature ID"
  if (identical(mode, "combined")) {
    path <- file.path(result$output_dir, "tables", "treated_vs_parent.by_pool.log2fc_z.csv")
    by_pool <- if (file.exists(path)) tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame()) else data.frame()
    if (!nrow(by_pool) || !all(c("feature_id", "log2FC", "z", "call") %in% names(by_pool))) {
      data$ranking_warning <- "Stored local z-score evidence is unavailable for this comparison."
      data$support_pass <- FALSE
      data$selected_min_support_pools <- min_support
      return(data)
    }
    by_pool$feature_id <- as.character(by_pool$feature_id); by_pool$log2FC <- fitness_ma_numeric(by_pool$log2FC); by_pool$z <- fitness_ma_numeric(by_pool$z)
    by_pool$call <- tolower(as.character(by_pool$call))
    cpm_ok <- rep(TRUE, nrow(by_pool))
    if (all(c("parent_cpm", "treated_cpm") %in% names(by_pool))) {
      cpm_pair <- cbind(fitness_ma_numeric(by_pool$parent_cpm), fitness_ma_numeric(by_pool$treated_cpm))
      cpm_ok <- is.finite(rowMeans(cpm_pair, na.rm = TRUE)) & rowSums(is.finite(cpm_pair)) > 0L & rowSums(pmax(cpm_pair, 0), na.rm = TRUE) > 0
    }
    valid <- by_pool[is.finite(by_pool$log2FC) & is.finite(by_pool$z) & cpm_ok, , drop = FALSE]
    count_by <- function(x, name) { if (!nrow(x)) return(data.frame(feature_id = character(), value = integer())); out <- aggregate(x$feature_id, list(feature_id = x$feature_id), length); names(out)[2] <- name; out }
    lfc_threshold <- fitness_ma_lfc_support_threshold()
    dep <- valid[valid$z < 0 & valid$log2FC < -lfc_threshold, , drop = FALSE]; enr <- valid[valid$z > 0 & valid$log2FC > lfc_threshold, , drop = FALSE]
    dep_stored <- valid[valid$call == "depleted", , drop = FALSE]; enr_stored <- valid[valid$call == "enriched", , drop = FALSE]
    valid$pool_cpm <- if (all(c("parent_cpm", "treated_cpm") %in% names(valid))) rowMeans(cbind(fitness_ma_numeric(valid$parent_cpm), fitness_ma_numeric(valid$treated_cpm)), na.rm = TRUE) else NA_real_
    evidence <- merge(count_by(valid, "valid_pool_n"), count_by(dep, "depleted_support_n"), by = "feature_id", all = TRUE)
    evidence <- merge(evidence, count_by(enr, "enriched_support_n"), by = "feature_id", all = TRUE)
    evidence <- merge(evidence, count_by(dep_stored, "depleted_stored_call_n"), by = "feature_id", all = TRUE)
    evidence <- merge(evidence, count_by(enr_stored, "enriched_stored_call_n"), by = "feature_id", all = TRUE)
    combined_lfc <- aggregate(valid$log2FC, list(feature_id = valid$feature_id), mean, na.rm = TRUE); names(combined_lfc)[2] <- "combined_log2FC_valid"
    combined_cpm <- aggregate(valid$pool_cpm, list(feature_id = valid$feature_id), mean, na.rm = TRUE); names(combined_cpm)[2] <- "combined_CPM_valid"
    evidence <- merge(evidence, combined_lfc, by = "feature_id", all = TRUE)
    evidence <- merge(evidence, combined_cpm, by = "feature_id", all = TRUE)
    if (nrow(evidence)) { for (nm in c("valid_pool_n", "depleted_support_n", "enriched_support_n", "depleted_stored_call_n", "enriched_stored_call_n")) evidence[[nm]][is.na(evidence[[nm]])] <- 0L }
    evidence$feature_id <- as.character(evidence$feature_id)
    metric <- function(x, label, fun = max) { if (!nrow(x)) return(data.frame(feature_id = character(), value = numeric())); x$abs_z <- abs(x$z); x$abs_lfc <- abs(x$log2FC); x$cpm_support <- if (all(c("parent_cpm", "treated_cpm") %in% names(x))) rowMeans(cbind(fitness_ma_numeric(x$parent_cpm), fitness_ma_numeric(x$treated_cpm)), na.rm = TRUE) else NA_real_; x$value <- if (label == "z") x$abs_z else if (label == "lfc") x$abs_lfc else x$cpm_support; aggregate(value ~ feature_id, x, fun, na.rm = TRUE) }
    selected <- if (identical(direction, "depleted")) valid[valid$z < 0, , drop = FALSE] else if (identical(direction, "enriched")) valid[valid$z > 0, , drop = FALSE] else valid
    for (spec in list(c("z", "rank_z_strength"))) { m <- metric(selected, spec[[1]], max); names(m)[2] <- spec[[2]]; evidence <- merge(evidence, m, by = "feature_id", all = TRUE) }
    evidence$rank_support_n <- if (identical(direction, "depleted")) evidence$depleted_support_n else if (identical(direction, "enriched")) evidence$enriched_support_n else pmax(evidence$depleted_support_n, evidence$enriched_support_n)
    data <- data[, setdiff(names(data), setdiff(names(evidence), "feature_id")), drop = FALSE]
    data <- merge(data, evidence, by = "feature_id", all.x = TRUE, sort = FALSE)
    for (nm in c("depleted_support_n", "enriched_support_n", "depleted_stored_call_n", "enriched_stored_call_n", "valid_pool_n", "rank_support_n")) data[[nm]][is.na(data[[nm]])] <- 0L
    data$log2fc <- ifelse(is.finite(data$combined_log2FC_valid), data$combined_log2FC_valid, data$log2fc)
    data$mean_abundance <- ifelse(is.finite(data$combined_CPM_valid), data$combined_CPM_valid, data$mean_abundance)
    data$rank_lfc_strength <- abs(fitness_ma_numeric(data$log2fc))
    data$rank_cpm_support <- fitness_ma_numeric(data$mean_abundance)
    data$rank_z_strength <- ifelse(is.finite(data$rank_z_strength), data$rank_z_strength, 0)
  } else {
    path <- file.path(result$output_dir, "tables", "treated_vs_parent.by_pool.log2fc_z.csv")
    by_pool <- if (file.exists(path)) tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame()) else data.frame()
    row <- by_pool[as.character(by_pool$contrast) == pair & as.character(by_pool$feature_id) %in% as.character(data$feature_id), , drop = FALSE]
    if (nrow(row)) { row <- row[match(as.character(data$feature_id), as.character(row$feature_id)), , drop = FALSE]; z <- fitness_ma_numeric(row$z); lfc <- fitness_ma_numeric(row$log2FC); call <- tolower(as.character(row$call)); cpm <- if (all(c("parent_cpm", "treated_cpm") %in% names(row))) rowMeans(cbind(fitness_ma_numeric(row$parent_cpm), fitness_ma_numeric(row$treated_cpm)), na.rm = TRUE) else rep(NA_real_, length(z)); usable <- is.finite(lfc) & is.finite(cpm) & cpm > 0; z_usable <- usable & is.finite(z); lfc_threshold <- fitness_ma_lfc_support_threshold(); dep_evidence <- z_usable & z < 0 & lfc < -lfc_threshold; enr_evidence <- z_usable & z > 0 & lfc > lfc_threshold; data$depleted_support_n <- as.integer(dep_evidence); data$enriched_support_n <- as.integer(enr_evidence); data$depleted_stored_call_n <- as.integer(call == "depleted"); data$enriched_stored_call_n <- as.integer(call == "enriched"); data$valid_pool_n <- as.integer(usable); data$rank_z_strength <- if (direction == "enriched") ifelse(z > 0 & is.finite(z), z, 0) else if (direction == "depleted") ifelse(z < 0 & is.finite(z), abs(z), 0) else ifelse(is.finite(z), abs(z), 0); data$rank_lfc_strength <- ifelse(usable, abs(lfc), NA); if (all(c("parent_cpm", "treated_cpm") %in% names(row))) data$rank_cpm_support <- ifelse(usable, cpm, NA); data$rank_support_n <- pmax(data$depleted_support_n, data$enriched_support_n); data$ranking_basis <- "ranked by directional log2FC magnitude, CPM/read support, local z-score support, and stable feature ID" }
  }
  data$candidate_log2FC_threshold <- fitness_ma_candidate_lfc_threshold()
  data$candidate_CPM_threshold <- fitness_ma_candidate_cpm_threshold()
  lfc_value <- fitness_ma_numeric(data$log2fc)
  data$passes_log2FC_threshold <- if (identical(direction, "depleted")) {
    is.finite(lfc_value) & lfc_value <= -data$candidate_log2FC_threshold
  } else if (identical(direction, "enriched")) {
    is.finite(lfc_value) & lfc_value >= data$candidate_log2FC_threshold
  } else if (identical(direction, "both")) {
    is.finite(lfc_value) & abs(lfc_value) >= data$candidate_log2FC_threshold
  } else {
    FALSE
  }
  data$passes_CPM_threshold <- is.finite(data$rank_cpm_support) & data$rank_cpm_support >= data$candidate_CPM_threshold
  direction_support <- if (identical(direction, "depleted")) data$depleted_support_n else if (identical(direction, "enriched")) data$enriched_support_n else ifelse(lfc_value < 0, data$depleted_support_n, data$enriched_support_n)
  data$rank_support_n <- direction_support
  data$support_pass <- data$passes_log2FC_threshold & data$passes_CPM_threshold
  if (min_support > 0L) data$support_pass <- data$support_pass & direction_support >= min_support
  data$selected_min_support_pools <- min_support
  data$concordance_filter <- if (min_support > 0L) paste0(">=", min_support, " supporting pool", ifelse(min_support == 1L, "", "s")) else "All support classes"
  mixed <- data$depleted_support_n > 0L & data$enriched_support_n > 0L
  support_label <- ifelse(direction_support >= 2L, "multi-pool concordant", ifelse(direction_support == 1L, "single-pool directional", "no directional pool support"))
  support_label[mixed] <- paste(support_label[mixed], "mixed-direction", sep = "; ")
  none_call <- tolower(as.character(data$call %||% "")) %in% c("", "none", "not_called", "not called")
  support_label[none_call] <- paste(support_label[none_call], "stored-call none", sep = "; ")
  weak_z <- !is.finite(data$rank_z_strength) | data$rank_z_strength <= 0
  support_label[weak_z] <- paste(support_label[weak_z], "statistically weak / near-background", sep = "; ")
  data$support_class <- support_label
  data$candidate_rank_score <- NA_real_; data$candidate_rank_order <- NA_integer_; data$overall_candidate_rank <- NA_integer_
  data$rank_z_strength <- ifelse(is.finite(data$rank_z_strength), data$rank_z_strength, 0)
  eligible <- which(data$passes_log2FC_threshold & data$passes_CPM_threshold)
  eligible <- eligible[is.finite(data$rank_lfc_strength[eligible]) & is.finite(data$rank_cpm_support[eligible])]
  if (length(eligible)) { ord <- eligible[order(-data$rank_lfc_strength[eligible], -data$rank_cpm_support[eligible], -data$rank_z_strength[eligible], as.character(data$feature_id[eligible]))]; data$overall_candidate_rank[ord] <- seq_along(ord); data$candidate_rank_order[ord] <- data$overall_candidate_rank[ord]; data$candidate_rank_score[ord] <- seq_along(ord) }
  data
}

fitness_ma_highlight_rows <- function(data, direction = "both", n = 10L, min_support = 0L) {
  if (!nrow(data) || !n) return(integer())
  if ("candidate_rank_order" %in% names(data)) {
    keep <- which(!is.na(data$support_pass) & data$support_pass)
    keep <- keep[is.finite(data$candidate_rank_order[keep]) & data$passes_log2FC_threshold[keep] & data$passes_CPM_threshold[keep]]
    return(head(keep[order(data$candidate_rank_order[keep])], n))
  }
  call <- tolower(as.character(data$call %||% "")); y <- fitness_ma_numeric(data$log2fc)
  d <- which(grepl("depleted", call) & is.finite(y)); e <- which(grepl("enriched", call) & is.finite(y))
  d <- d[order(y[d], decreasing = FALSE)]; e <- e[order(y[e], decreasing = TRUE)]
  if (identical(direction, "depleted")) return(head(d, n))
  if (identical(direction, "enriched")) return(head(e, n))
  if (identical(direction, "both")) return(unique(c(head(d, n), head(e, n))))
  integer()
}

fitness_ma_annotation_rows <- function(data, highlighted, mode = "top", custom = "") {
  if (!nrow(data)) return(list(rows = integer(), unmatched = character()))
  queries <- trimws(unlist(strsplit(custom %||% "", "[,;[:space:][:cntrl:]]+"))); queries <- unique(queries[nzchar(queries)])
  custom_rows <- integer(); unmatched <- character()
  if (mode %in% c("custom", "top_custom") && length(queries)) for (q in queries) {
    fields <- intersect(c("feature_id", "original_feature_id", "label", "cagl_display_id",
                          "gene_display_name", "scer_gene_name", "scer_gene_id", "qng_id",
                          "cgla_common_name"), names(data))
    hit <- integer()
    if (length(fields)) {
      matches <- vapply(data[fields], function(x) tolower(as.character(x)) == tolower(q),
                        logical(nrow(data)))
      hit <- which(rowSums(as.matrix(matches), na.rm = TRUE) > 0)
    }
    if (length(hit)) custom_rows <- c(custom_rows, hit[[1]]) else unmatched <- c(unmatched, q)
  }
  top_rows <- if (mode %in% c("top", "top_custom")) highlighted else integer()
  list(rows = unique(c(top_rows, custom_rows)), unmatched = unmatched)
}

fitness_ma_label_table <- function(data, highlighted, annotation_mode = "top", custom = "", mode = "combined", pair = "", max_labels = 25L) {
  ann <- fitness_ma_annotation_rows(data, highlighted, annotation_mode, custom)
  rows <- ann$rows
  if (length(rows) > max_labels) rows <- rows[seq_len(max_labels)]
  if (!length(rows)) return(list(data = data.frame(), unmatched = ann$unmatched, limited = FALSE, requested = length(ann$rows)))
  x <- log10(data$mean_abundance[rows] + 1); y <- data$log2fc[rows]
  label_data <- data.frame(feature_id = data$feature_id[rows], display_label = data$label[rows], x = x, y = y,
                           x_label = x,
                           y_label = y,
                           side = 1,
                           annotation_source = ifelse(rows %in% highlighted, "top_hit", "custom"),
                           comparison_view = if (identical(mode, "individual")) "individual_pair" else "combined",
                           pair = if (nzchar(pair)) pair else NA_character_, stringsAsFactors = FALSE)
  list(data = label_data, unmatched = ann$unmatched, limited = length(ann$rows) > nrow(label_data), requested = length(ann$rows))
}

fitness_ma_pack_label_y <- function(target_y, labels, cex, y_min, y_max, yr) {
  n <- length(target_y)
  if (!n) return(target_y)
  if (n == 1L) return(min(max(target_y, y_min), y_max))
  label_heights <- vapply(labels, strheight, numeric(1), cex = cex, font = 2)
  label_height <- max(label_heights[is.finite(label_heights)], 0.025 * yr)
  min_gap <- max(label_height * 1.45, 0.035 * yr)
  available <- y_max - y_min
  if (!is.finite(available) || available <= 0) return(target_y)
  if ((n - 1L) * min_gap > available) {
    min_gap <- available / max(1L, n - 1L)
  }
  ord <- order(target_y)
  y <- target_y[ord]
  for (i in seq.int(2L, n)) {
    y[[i]] <- max(y[[i]], y[[i - 1L]] + min_gap)
  }
  if (y[[n]] > y_max) y <- y - (y[[n]] - y_max)
  for (i in seq.int(n - 1L, 1L)) {
    y[[i]] <- min(y[[i]], y[[i + 1L]] - min_gap)
  }
  if (y[[1L]] < y_min) y <- y + (y_min - y[[1L]])
  y <- pmin(pmax(y, y_min), y_max)
  out <- target_y
  out[ord] <- y
  out
}

fitness_ma_place_labels <- function(label_data, cex = 1) {
  if (!nrow(label_data)) return(label_data)
  usr <- par("usr")
  xr <- diff(usr[1:2]); yr <- diff(usr[3:4])
  if (!is.finite(xr) || xr <= 0) xr <- 1
  if (!is.finite(yr) || yr <= 0) yr <- 1
  edge_x <- 0.025 * xr
  edge_y <- 0.050 * yr
  stagger_pattern <- data.frame(
    side = c(1, -1, 1, -1, 1, -1, 1, -1, 1, -1),
    x_frac = c(0.120, 0.120, 0.155, 0.155, 0.190, 0.190, 0.135, 0.135, 0.170, 0.170),
    y_frac = c(0.000, 0.000, 0.018, -0.018, -0.032, 0.032, 0.050, -0.050, -0.065, 0.065),
    stringsAsFactors = FALSE
  )
  label_order <- order(label_data$x, label_data$y)
  for (rank in seq_along(label_order)) {
    row_ix <- label_order[[rank]]
    pattern_ix <- ((rank - 1L) %% nrow(stagger_pattern)) + 1L
    pattern <- stagger_pattern[pattern_ix, , drop = FALSE]
    side <- pattern$side[[1L]]
    label <- as.character(label_data$display_label[[row_ix]])
    text_width <- strwidth(label, cex = cex, font = 2)
    proposed_x <- label_data$x[[row_ix]] + side * pattern$x_frac[[1L]] * xr
    if (side > 0) {
      proposed_x <- min(proposed_x, usr[[2L]] - edge_x - text_width)
      proposed_x <- max(proposed_x, label_data$x[[row_ix]] + 0.085 * xr)
    } else {
      proposed_x <- max(proposed_x, usr[[1L]] + edge_x + text_width)
      proposed_x <- min(proposed_x, label_data$x[[row_ix]] - 0.085 * xr)
    }
    label_data$x_label[[row_ix]] <- proposed_x
    label_data$y_label[[row_ix]] <- label_data$y[[row_ix]] + pattern$y_frac[[1L]] * yr
    label_data$side[[row_ix]] <- side
  }
  for (side in c(-1, 1)) {
    idx <- which(label_data$side == side)
    if (length(idx)) {
      label_data$y_label[idx] <- fitness_ma_pack_label_y(
        target_y = label_data$y_label[idx],
        labels = as.character(label_data$display_label[idx]),
        cex = cex,
        y_min = usr[[3L]] + edge_y,
        y_max = usr[[4L]] - edge_y,
        yr = yr
      )
    }
  }
  label_data
}

fitness_ma_top_hits_table <- function(result, mode = "combined", pair = "", direction = "both", n = 10L,
                                     annotation_mode = "top", custom = "", min_support = 0L) {
  n <- suppressWarnings(as.integer(n %||% 0L))
  if (is.na(n) || n < 0L) n <- 0L
  min_support <- suppressWarnings(as.integer(min_support %||% 0L))
  if (is.na(min_support) || min_support < 0L) min_support <- 0L
  ranked <- fitness_ma_rank_data(result, mode, pair, direction, min_support)
  if (!nrow(ranked) || identical(direction, "none") || n <= 0L)
    return(ranked[0, , drop = FALSE])
  highlighted <- fitness_ma_highlight_rows(ranked, direction, n, min_support)
  labels <- fitness_ma_label_table(ranked, highlighted, annotation_mode, custom, mode, pair)
  if (!length(highlighted)) return(ranked[0, , drop = FALSE])
  out <- ranked[highlighted, , drop = FALSE]
  out$highlighted <- TRUE
  out$labeled <- seq_len(nrow(out)) %in% match(labels$data$feature_id, out$feature_id)
  out$annotation_source <- ifelse(out$labeled, "top_hit", "none")
  out$comparison_view <- if (identical(mode, "individual")) "individual_pair" else "combined"
  out$pair <- if (nzchar(pair)) pair else NA_character_
  out$hit_direction <- direction
  out$stored_call_or_final_call <- out$call
  if (identical(mode, "combined")) {
    pool_path <- file.path(result$output_dir, "tables", "treated_vs_parent.by_pool.log2fc_z.csv")
    pool_data <- if (file.exists(pool_path)) tryCatch(read.csv(pool_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame()) else data.frame()
    if (nrow(pool_data)) {
      pz <- fitness_ma_numeric(pool_data$z); pl <- fitness_ma_numeric(pool_data$log2FC)
      pkeep <- if (identical(direction, "depleted")) pz < 0 & pl < 0 else if (identical(direction, "enriched")) pz > 0 & pl > 0 else pz * pl > 0
      plabel <- if ("pool" %in% names(pool_data)) paste0("pool", pool_data$pool) else as.character(pool_data$contrast)
      pool_data$.support_label <- plabel; pool_data$.support_keep <- pkeep
      support_map <- tapply(pool_data$.support_label[pool_data$.support_keep], pool_data$feature_id[pool_data$.support_keep], function(v) paste(unique(v), collapse = ";"))
      out$supporting_pool_ids <- unname(support_map[as.character(out$feature_id)]); out$supporting_pool_ids[is.na(out$supporting_pool_ids)] <- ""
    } else out$supporting_pool_ids <- ""
  } else out$supporting_pool_ids <- "Current pair only"
  out$direction_specific_lfc_strength <- out$rank_lfc_strength
  out$direction_specific_cpm_support <- out$rank_cpm_support
  out$direction_specific_z_support <- out$rank_z_strength
  out$selected_min_support_pools <- min_support
  out$concordance_filter <- if (min_support > 0L) paste0(">=", min_support, " supporting pool", ifelse(min_support == 1L, "", "s")) else "All support classes"
  out$overall_candidate_rank <- ranked$overall_candidate_rank[highlighted]
  out$selected_rank <- seq_len(nrow(out))
  out$candidate_rank_order <- ranked$candidate_rank_order[highlighted]
  out
}

fitness_ma_supporting_pool_ids <- function(result, feature_id, direction = "both") {
  path <- file.path(result$output_dir, "tables", "treated_vs_parent.by_pool.log2fc_z.csv")
  if (!file.exists(path)) return("")
  x <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (!nrow(x)) return("")
  x <- x[as.character(x$feature_id) == as.character(feature_id), , drop = FALSE]
  if (!nrow(x)) return("")
  z <- fitness_ma_numeric(x$z); lfc <- fitness_ma_numeric(x$log2FC)
  keep <- if (identical(direction, "depleted")) z < 0 & lfc < 0 else if (identical(direction, "enriched")) z > 0 & lfc > 0 else z * lfc > 0
  label <- if ("pool" %in% names(x)) paste0("pool", x$pool) else as.character(x$contrast)
  paste(unique(label[keep & !is.na(label) & nzchar(label)]), collapse = ";")
}

plot_fitness_ma <- function(result, mode = "combined", pair = "", direction = "both", n = 10L,
                            annotation_mode = "top", custom = "", text_size = "medium", grid = TRUE,
                            point_size = 1.5, min_support = 1L, repo_root = NULL) {
  data <- fitness_ma_rank_data(result, mode, pair, direction, min_support)
  if (!nrow(data)) return(qc_plot_empty("No MA data are available for this comparison."))
  if (!identical(direction, "none") && any(nzchar(as.character(data$ranking_warning %||% ""))))
    return(qc_plot_empty(unique(as.character(data$ranking_warning))[1]))
  data <- fitness_ma_apply_display_labels(data, repo_root)
  data <- data[is.finite(data$mean_abundance) & is.finite(data$log2fc), , drop = FALSE]
  if (!nrow(data)) return(qc_plot_empty("No finite abundance/log2FC values are available for this comparison."))
  highlighted <- fitness_ma_highlight_rows(data, direction, as.integer(n %||% 10L), min_support)
  label_table <- fitness_ma_label_table(data, highlighted, annotation_mode, custom, mode, pair)
  scales <- qc_plot_text_sizes(); old <- qc_plot_par(mar = c(8.3, 6.2, 5.4, 2.2), mgp = c(3.2, 0.9, 0)); on.exit(par(old), add = TRUE)
  depleted_col <- "#2f6fb5"; enriched_col <- "#c83f3f"
  selected_hit_colors <- rep(NA_character_, nrow(data))
  if (length(highlighted)) {
    selected_hit_colors[highlighted] <- if (identical(direction, "depleted")) depleted_col else if (identical(direction, "enriched")) enriched_col else ifelse(data$log2fc[highlighted] < 0, depleted_col, enriched_col)
  }
  display_point_size <- qc_plot_point_cex(point_size)
  plot_x <- log10(data$mean_abundance + 1)
  xlim <- range(plot_x, finite = TRUE)
  if (diff(xlim) <= 0 || !all(is.finite(xlim))) xlim <- xlim + c(-0.5, 0.5)
  if (nrow(label_table$data)) xlim <- xlim + c(-0.16, 0.16) * diff(xlim)
  plot(plot_x, data$log2fc, pch = 16, cex = display_point_size, col = adjustcolor("grey70", alpha.f = 0.7),
       xlim = xlim,
       xlab = "M: mean abundance, log10(CPM + 1)", ylab = "A: fitness effect, log2FC treated/control",
       main = if (identical(mode, "individual")) "MA plot — individual treated-control pair" else "MA plot — combined across pools",
       cex.main = scales$title, cex.lab = scales$axis, cex.axis = scales$sample)
  if (isTRUE(grid)) grid(col = "grey85", lty = 1)
  if (length(highlighted)) points(log10(data$mean_abundance[highlighted] + 1), data$log2fc[highlighted], pch = 21, bg = selected_hit_colors[highlighted], col = "black", lwd = qc_plot_lwd, cex = display_point_size + qc_plot_point_cex(0.55))
  if (nrow(label_table$data)) {
    label_table$data <- fitness_ma_place_labels(label_table$data, cex = scales$key)
    segments(label_table$data$x, label_table$data$y, label_table$data$x_label, label_table$data$y_label,
             col = adjustcolor("#68727d", alpha.f = 0.65), lwd = qc_plot_lwd * 0.5)
    text(label_table$data$x_label, label_table$data$y_label, labels = label_table$data$display_label,
         cex = scales$key, font = 2, xpd = NA,
         adj = ifelse(label_table$data$side > 0, 0, 1))
  }
  present <- c(depleted = any(selected_hit_colors[highlighted] == depleted_col), enriched = any(selected_hit_colors[highlighted] == enriched_col))
  legend_labels <- c(depleted = "Depleted selected hit", enriched = "Enriched selected hit")
  legend_cols <- c(depleted = depleted_col, enriched = enriched_col)
  if (any(present)) legend("topright", legend = unname(legend_labels[present]), col = unname(legend_cols[present]), pch = 21, pt.bg = unname(legend_cols[present]), bty = "n", cex = scales$key)
  note_line <- 5.75
  if (length(highlighted)) {
    mtext("Depleted/enriched colors reflect stored calls from local-abundance z-scores.", side = 1, line = note_line, cex = scales$key, col = "#666666")
    note_line <- note_line + 0.9
  }
  if (length(highlighted) > 10L && nrow(label_table$data)) {
    mtext("More than 10 labels may overlap; reduce Number of hits for cleaner labels.", side = 1, line = note_line, cex = scales$key, col = "#666666")
    note_line <- note_line + 0.9
  }
  if (length(label_table$unmatched)) {
    mtext(paste("Unmatched:", paste(label_table$unmatched, collapse = ", ")), side = 1, line = note_line, cex = scales$key, col = "#666666")
    note_line <- note_line + 0.9
  }
  if (isTRUE(label_table$limited)) mtext(sprintf("Showing labels for the top %d of %d highlighted features to reduce overlap.", nrow(label_table$data), label_table$requested), side = 1, line = note_line, cex = scales$key, col = "#666666")
  invisible(data)
}

fitness_selected_hit_heatmap_data <- function(result, selected_hits, mode = "combined", pair = "") {
  empty <- list(data = data.frame(), matrix = matrix(numeric(), nrow = 0, ncol = 0), row_labels = character(), column_labels = character())
  if (is.null(result) || !nzchar(result$output_dir %||% "") || !nrow(selected_hits)) return(empty)
  path <- file.path(result$output_dir, "tables", "treated_vs_parent.by_pool.log2fc_z.csv")
  if (!file.exists(path)) return(empty)
  by_pool <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (!nrow(by_pool) || !all(c("feature_id", "contrast", "log2FC") %in% names(by_pool))) return(empty)
  by_pool$feature_id <- as.character(by_pool$feature_id)
  by_pool$contrast <- as.character(by_pool$contrast)
  by_pool$log2FC <- fitness_ma_numeric(by_pool$log2FC)
  if (identical(mode, "individual") && nzchar(pair)) by_pool <- by_pool[by_pool$contrast == pair, , drop = FALSE]
  if (!nrow(by_pool)) return(empty)
  selected_hits$feature_id <- as.character(selected_hits$feature_id)
  selected <- selected_hits[!duplicated(selected_hits$feature_id), , drop = FALSE]
  by_pool <- by_pool[by_pool$feature_id %in% selected$feature_id, , drop = FALSE]
  if (!nrow(by_pool)) return(empty)
  contrast_order <- unique(by_pool$contrast[order(if ("pool" %in% names(by_pool)) suppressWarnings(as.numeric(by_pool$pool)) else seq_len(nrow(by_pool)), by_pool$contrast)])
  feature_order <- selected$feature_id
  z <- matrix(NA_real_, nrow = length(feature_order), ncol = length(contrast_order), dimnames = list(feature_order, contrast_order))
  for (i in seq_len(nrow(by_pool))) z[by_pool$feature_id[[i]], by_pool$contrast[[i]]] <- by_pool$log2FC[[i]]
  row_labels <- fitness_ma_first_nonempty(
    if ("gene_display_name" %in% names(selected)) selected$gene_display_name else "",
    if ("cagl_display_id" %in% names(selected)) selected$cagl_display_id else "",
    if ("label" %in% names(selected)) selected$label else "",
    selected$feature_id
  )
  row_labels[is.na(row_labels) | !nzchar(trimws(row_labels))] <- selected$feature_id[is.na(row_labels) | !nzchar(trimws(row_labels))]
  column_labels <- vapply(contrast_order, function(id) {
    row <- by_pool[match(id, by_pool$contrast), , drop = FALSE]
    if ("pool" %in% names(row) && !is.na(row$pool[[1]])) paste0("pool", row$pool[[1]]) else id
  }, "")
  flat <- data.frame(
    gene = rep(row_labels, each = length(contrast_order)),
    rank = rep(selected$selected_rank %||% seq_len(nrow(selected)), each = length(contrast_order)),
    hit_direction = rep(selected$hit_direction %||% "", each = length(contrast_order)),
    support_class = rep(selected$support_class %||% "", each = length(contrast_order)),
    feature_id = rep(feature_order, each = length(contrast_order)),
    contrast = rep(contrast_order, times = length(feature_order)),
    log2FC = as.vector(t(z)),
    stringsAsFactors = FALSE
  )
  list(data = flat, matrix = z, row_labels = row_labels, column_labels = column_labels)
}

plot_fitness_selected_hit_heatmap <- function(result, selected_hits, mode = "combined", pair = "", text_size = "medium") {
  hm <- fitness_selected_hit_heatmap_data(result, selected_hits, mode, pair)
  if (!nrow(hm$data) || !length(hm$matrix)) return(qc_plot_empty("No selected top hits are available for the log2FC heatmap."))
  z <- hm$matrix
  rownames(z) <- hm$row_labels
  colnames(z) <- hm$column_labels
  scales <- qc_plot_text_sizes()
  max_row_label <- max(nchar(rownames(z)), na.rm = TRUE)
  row_cex <- min(scales$sample, if (max_row_label > 18) scales$sample * 18 / max_row_label else scales$sample)
  left_margin <- max(9, min(30, max_row_label * row_cex * 0.78 + 3))
  old <- qc_plot_par(mar = c(7.5, left_margin, 4.5, 2))
  on.exit(par(old), add = TRUE)
  finite <- z[is.finite(z)]
  max_abs <- if (length(finite)) max(abs(finite), na.rm = TRUE) else 1
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  palette <- grDevices::colorRampPalette(c("#2f6fb5", "white", "#c83f3f"))(100)
  image(t(z[nrow(z):1, , drop = FALSE]), col = palette, breaks = breaks, axes = FALSE,
        xlab = "", ylab = "", main = if (identical(mode, "individual") && nzchar(pair)) paste("Top selected hits log2FC heatmap —", pair) else "Top selected hits log2FC heatmap — combined across pools",
        cex.main = scales$title)
  axis(1, at = seq(0, 1, length.out = ncol(z)), labels = colnames(z), las = 2, cex.axis = scales$sample)
  axis(2, at = seq(0, 1, length.out = nrow(z)), labels = rev(rownames(z)), las = 1, cex.axis = row_cex)
  box()
  mtext("Rows = currently selected hits; columns = matched treated-control pool contrasts; values = log2FC treated/control.", side = 1, line = 5.4, cex = scales$key, col = "#666666")
  invisible(hm$data)
}
