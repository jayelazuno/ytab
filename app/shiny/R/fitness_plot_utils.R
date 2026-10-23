fitness_ma_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
fitness_ma_candidate_lfc_threshold <- function() 1
fitness_ma_candidate_cpm_threshold <- function() 1

# Results are read as the pipeline wrote them (ytab.fitness.contract); the app derives nothing. Gene values and
# abundances are in the result's value_unit (CPM, or the normalization's reads per million).
fitness_abundance_unit <- function(result) if (nzchar(result$value_unit %||% "")) result$value_unit else "CPM"
fitness_gene_values <- function(result) {
  path <- file.path(result$output_dir %||% "", "tables", "gene_values.long.csv")
  x <- if (!is.null(result) && file.exists(path)) tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame()) else data.frame()
  if (nrow(x) && "value" %in% names(x)) x$cpm <- fitness_ma_numeric(x$value)
  x
}

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
  if (!"mean_abundance" %in% names(summary)) return(empty)
  out <- data.frame(feature_id = id, label = labels, mean_abundance = fitness_ma_numeric(summary$mean_abundance),
                    log2fc = fitness_ma_numeric(if ("call_log2FC" %in% names(summary)) summary$call_log2FC else summary$mean_log2FC),
                    call = if ("final_call" %in% names(summary)) as.character(summary$final_call) else "none",
                    stringsAsFactors = FALSE)
  out
}

# Top hits as ranked by the pipeline (tables/top_hit_evidence*.csv, ytab.fitness.zscore.hit_evidence). The app only
# applies the viewer's choices: direction, minimum supporting pools and how many hits (fitness_ma_highlight_rows).
fitness_ma_rank_data <- function(result, mode = "combined", pair = "", direction = "both", min_support = 0L) {
  data <- fitness_ma_data(result, mode, pair)
  if (!nrow(data)) return(data)
  individual <- identical(mode, "individual")
  if (individual && !nzchar(pair)) pair <- unname(fitness_ma_pair_choices(result)[[1]] %||% "")
  path <- file.path(result$output_dir %||% "", "tables", if (individual) "top_hit_evidence.by_pair.csv" else "top_hit_evidence.csv")
  hits <- if (file.exists(path)) tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame()) else data.frame()
  if (individual && nrow(hits)) hits <- hits[as.character(hits$contrast) == pair, , drop = FALSE]
  min_support <- suppressWarnings(as.integer(min_support %||% 0L)); if (is.na(min_support) || min_support < 0L) min_support <- 0L
  if (!nrow(hits)) {
    data$ranking_warning <- "This result has no top-hit ranking; it was made by an earlier version of the pipeline. Run it again."
    data$support_pass <- FALSE; data$candidate_rank_order <- NA_integer_; data$selected_min_support_pools <- min_support
    return(data)
  }
  hits <- hits[match(as.character(data$feature_id), as.character(hits$feature_id)), , drop = FALSE]
  key <- if (direction %in% c("depleted", "enriched")) direction else "both"
  num <- function(x) fitness_ma_numeric(x)
  data$log2fc <- num(hits$log2FC); data$mean_abundance <- num(hits$mean_abundance)
  data$valid_pool_n <- as.integer(hits$n_pairs_valid); data$depleted_support_n <- as.integer(hits$n_pairs_depleted_support)
  data$enriched_support_n <- as.integer(hits$n_pairs_enriched_support); data$depleted_stored_call_n <- as.integer(hits$n_pairs_depleted_called)
  data$enriched_stored_call_n <- as.integer(hits$n_pairs_enriched_called)
  data$rank_lfc_strength <- abs(data$log2fc); data$rank_cpm_support <- data$mean_abundance
  data$rank_z_strength <- num(hits[[paste0("z_strength_", key)]])
  data$supporting_pool_ids <- ifelse(data$log2fc < 0, hits$pools_depleted_support, hits$pools_enriched_support)
  data$supporting_pool_ids[is.na(data$supporting_pool_ids)] <- ""
  data$ranking_basis <- "ranked by the pipeline: called genes with at least a 2-fold change, by |z|, then |log2FC|, then feature ID (Gale et al. 2020)"
  settings <- result$settings %||% list()
  data$candidate_log2FC_threshold <- as.numeric(settings$hit_min_log2FC %||% fitness_ma_candidate_lfc_threshold())
  data$candidate_CPM_threshold <- as.numeric(settings$hit_min_abundance %||% fitness_ma_candidate_cpm_threshold())
  lfc <- data$log2fc
  data$passes_log2FC_threshold <- if (identical(direction, "none")) FALSE else is.finite(lfc) &
    (if (key == "depleted") lfc <= -data$candidate_log2FC_threshold else if (key == "enriched") lfc >= data$candidate_log2FC_threshold else abs(lfc) >= data$candidate_log2FC_threshold)
  data$passes_CPM_threshold <- hits$passes_abundance %in% c(TRUE, "True", "TRUE")
  support <- if (key == "depleted") data$depleted_support_n else if (key == "enriched") data$enriched_support_n else ifelse(lfc < 0, data$depleted_support_n, data$enriched_support_n)
  support[is.na(support)] <- 0L; data$rank_support_n <- support
  rank <- num(hits[[paste0("rank_", key)]]); if (identical(direction, "none")) rank[] <- NA_real_
  data$overall_candidate_rank <- as.integer(rank); data$candidate_rank_order <- as.integer(rank); data$candidate_rank_score <- rank
  data$support_pass <- is.finite(rank) & (min_support == 0L | support >= min_support)
  data$selected_min_support_pools <- min_support
  data$concordance_filter <- if (min_support > 0L) paste0(">=", min_support, " supporting pool", ifelse(min_support == 1L, "", "s")) else "All support classes"
  # words for the pipeline's counts (display only)
  mixed <- data$depleted_support_n > 0L & data$enriched_support_n > 0L
  label <- ifelse(support >= 2L, "multi-pool concordant", ifelse(support == 1L, "single-pool directional", "no directional pool support"))
  label[mixed %in% TRUE] <- paste(label[mixed %in% TRUE], "mixed-direction", sep = "; ")
  none_call <- tolower(as.character(data$call %||% "")) %in% c("", "none", "not_called", "not called")
  label[none_call] <- paste(label[none_call], "stored-call none", sep = "; ")
  weak <- !is.finite(data$rank_z_strength) | data$rank_z_strength <= 0
  label[weak] <- paste(label[weak], "statistically weak / near-background", sep = "; ")
  data$support_class <- label
  data$rank_z_strength[!is.finite(data$rank_z_strength)] <- 0
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
  # Push each label only as far as it must: upward from the lower edge, then downward from the upper edge. (Shifting
  # the whole side at once coupled labels near the top with labels near the bottom and clamped some onto each other.)
  y[[1L]] <- max(y[[1L]], y_min)
  for (i in seq.int(2L, n)) y[[i]] <- max(y[[i]], y[[i - 1L]] + min_gap)
  y[[n]] <- min(y[[n]], y_max)
  for (i in seq.int(n - 1L, 1L)) y[[i]] <- min(y[[i]], y[[i + 1L]] - min_gap)
  out <- target_y
  out[ord] <- y
  out
}

# Separates labels whose text boxes overlap, whichever side of their point they are on (the side packing only sees
# one side): each overlapping pair moves apart vertically, kept inside [y_min, y_max]; repeats until none overlap.
fitness_ma_resolve_overlaps <- function(label_data, cex, y_min, y_max, yr, iterations = 300L) {
  n <- nrow(label_data); if (n < 2L) return(label_data)
  w <- vapply(as.character(label_data$display_label), strwidth, numeric(1), cex = cex, font = 2)
  h <- max(vapply(as.character(label_data$display_label), strheight, numeric(1), cex = cex, font = 2), 0.025 * yr) * 1.35
  left <- ifelse(label_data$side > 0, label_data$x_label, label_data$x_label - w); right <- left + w
  pad <- 0.01 * diff(range(c(left, right)))
  y <- label_data$y_label
  for (k in seq_len(iterations)) {
    moved <- FALSE
    for (i in seq_len(n - 1L)) for (j in seq.int(i + 1L, n)) {
      if (left[[i]] - pad >= right[[j]] || left[[j]] - pad >= right[[i]]) next
      gap <- abs(y[[i]] - y[[j]]); if (gap >= h) next
      push <- (h - gap) / 2 + 1e-9 * yr; up <- if (y[[i]] > y[[j]] || (y[[i]] == y[[j]] && i < j)) i else j; down <- if (up == i) j else i
      y[[up]] <- y[[up]] + push; y[[down]] <- y[[down]] - push
      if (y[[up]] > y_max) { y[[down]] <- y[[down]] - (y[[up]] - y_max); y[[up]] <- y_max }
      if (y[[down]] < y_min) { y[[up]] <- y[[up]] + (y_min - y[[down]]); y[[down]] <- y_min }
      moved <- TRUE
    }
    if (!moved) break
  }
  label_data$y_label <- pmin(pmax(y, y_min), y_max)
  label_data
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
  fitness_ma_resolve_overlaps(label_data, cex, usr[[3L]] + edge_y, usr[[4L]] - edge_y, yr)
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
  if (!identical(mode, "combined")) out$supporting_pool_ids <- "Current pair only"
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
  path <- file.path(result$output_dir %||% "", "tables", "top_hit_evidence.csv")
  if (!file.exists(path)) return("")
  x <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  x <- x[as.character(x$feature_id) == as.character(feature_id), , drop = FALSE]
  if (!nrow(x)) return("")
  ids <- c(if (direction %in% c("depleted", "both")) x$pools_depleted_support, if (direction %in% c("enriched", "both")) x$pools_enriched_support)
  paste(unique(unlist(strsplit(ids[!is.na(ids)], ";"))), collapse = ";")
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
       xlab = sprintf("A: mean abundance, log10(%s + 1)", fitness_abundance_unit(result)), ylab = "M: fitness effect, log2FC treated/control",
       main = if (identical(mode, "individual")) "MA plot — individual treated-control pair" else "MA plot — combined across pools",
       cex.main = scales$title, cex.lab = scales$axis, cex.axis = scales$sample)
  if (isTRUE(grid)) grid(col = "grey85", lty = 1)
  ytab_base_axes()
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
