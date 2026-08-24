fitness_condition_control_plot_data <- function(result, mode = "combined", pair = "") {
  if (is.null(result) || !nzchar(result$output_dir %||% "")) return(data.frame())
  path <- file.path(result$output_dir, "tables", "treated_vs_parent.by_pool.log2fc_z.csv")
  if (!file.exists(path)) return(data.frame())
  by_pool <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  required <- c("feature_id", "parent_cpm", "treated_cpm", "log2FC")
  if (!nrow(by_pool) || !all(required %in% names(by_pool))) return(data.frame())
  if (identical(mode, "individual")) {
    if (!nzchar(pair)) pair <- fitness_ma_pair_choices(result)[[1]] %||% ""
    by_pool <- by_pool[as.character(by_pool$contrast) == pair, , drop = FALSE]
  }
  if (!nrow(by_pool)) return(data.frame())
  by_pool$feature_id <- as.character(by_pool$feature_id)
  by_pool$control_abundance <- qc_plot_numeric(by_pool$parent_cpm)
  by_pool$treated_abundance <- qc_plot_numeric(by_pool$treated_cpm)
  by_pool$log2fc <- qc_plot_numeric(by_pool$log2FC)
  valid <- by_pool[is.finite(by_pool$control_abundance) & is.finite(by_pool$treated_abundance) &
                     is.finite(by_pool$log2fc), , drop = FALSE]
  if (!nrow(valid)) return(data.frame())
  if (identical(mode, "individual")) {
    out <- valid
    out$valid_pool_n <- 1L
  } else {
    control <- aggregate(valid$control_abundance, list(feature_id = valid$feature_id), mean, na.rm = TRUE)
    treated <- aggregate(valid$treated_abundance, list(feature_id = valid$feature_id), mean, na.rm = TRUE)
    valid_n <- aggregate(valid$feature_id, list(feature_id = valid$feature_id), length)
    names(control)[2] <- "control_abundance"; names(treated)[2] <- "treated_abundance"; names(valid_n)[2] <- "valid_pool_n"
    out <- merge(control, treated, by = "feature_id", all = TRUE)
    out <- merge(out, valid_n, by = "feature_id", all = TRUE)
    meta_cols <- intersect(c("feature_id", "call", "standard_name", "common_name"), names(valid))
    meta <- valid[!duplicated(valid$feature_id), meta_cols, drop = FALSE]
    out <- merge(out, meta, by = "feature_id", all.x = TRUE, sort = FALSE)
  }
  out$call <- as.character(out$call %||% "")
  out$label <- if ("common_name" %in% names(out)) as.character(out$common_name) else rep("", nrow(out))
  if ("standard_name" %in% names(out)) out$label[is.na(out$label) | !nzchar(trimws(out$label))] <- as.character(out$standard_name[is.na(out$label) | !nzchar(trimws(out$label))])
  out$label[is.na(out$label) | !nzchar(trimws(out$label))] <- out$feature_id[is.na(out$label) | !nzchar(trimws(out$label))]
  out
}

fitness_condition_control_highlights <- function(data, n = 10L) {
  if (!nrow(data) || !"mean_log2FC" %in% names(data)) return(integer())
  depleted <- which(grepl("depleted", data$final_call %||% "", ignore.case = TRUE) &
                      is.finite(data$mean_log2FC))
  enriched <- which(grepl("enriched", data$final_call %||% "", ignore.case = TRUE) &
                      is.finite(data$mean_log2FC))
  depleted <- depleted[order(data$mean_log2FC[depleted], decreasing = FALSE)]
  enriched <- enriched[order(data$mean_log2FC[enriched], decreasing = TRUE)]
  unique(c(head(depleted, n), head(enriched, n)))
}

plot_fitness_condition_control_scatter <- function(result, mode = "combined", pair = "", direction = "both", n = 10L,
                                                   annotation_mode = "top", custom = "", text_size = "medium", grid = TRUE,
                                                   point_size = 1.5, min_support = 0L, repo_root = NULL) {
  data <- fitness_condition_control_plot_data(result, mode, pair)
  if (!nrow(data)) return(qc_plot_empty("Condition-control log-log scatter data are not available for this selected fitness result."))
  ranked <- fitness_ma_rank_data(result, mode, pair, direction, min_support)
  if (!nrow(ranked)) return(qc_plot_empty("Selected-hit ranking data are not available for this selected fitness result."))
  ranked <- fitness_ma_apply_display_labels(ranked, repo_root)
  highlighted_ranked <- fitness_ma_highlight_rows(ranked, direction, as.integer(n %||% 10L), min_support)
  highlighted_features <- as.character(ranked$feature_id[highlighted_ranked])
  label_info <- fitness_ma_label_table(ranked, highlighted_ranked, annotation_mode, custom, mode, pair)
  label_features <- as.character(label_info$data$feature_id %||% character())

  control <- qc_plot_numeric(data$control_abundance)
  treated <- qc_plot_numeric(data$treated_abundance)
  keep <- is.finite(control) & is.finite(treated) & control >= 0 & treated >= 0
  if (!any(keep)) return(qc_plot_empty("Condition-versus-control read data are not available."))
  x <- log10(control + 1)
  y <- log10(treated + 1)
  highlight <- which(as.character(data$feature_id) %in% highlighted_features & keep)
  label_rows <- which(as.character(data$feature_id) %in% label_features & keep)
  label_data <- data.frame()
  if (length(label_rows)) {
    idx <- match(as.character(data$feature_id[label_rows]), as.character(label_info$data$feature_id))
    label_data <- data.frame(
      feature_id = as.character(data$feature_id[label_rows]),
      display_label = label_info$data$display_label[idx],
      x = x[label_rows], y = y[label_rows],
      x_label = x[label_rows], y_label = y[label_rows],
      side = 1,
      stringsAsFactors = FALSE
    )
  }

  scales <- qc_plot_text_sizes()
  old <- qc_plot_par(mar = c(8.3, 6.2, 5.4, 2.2), mgp = c(3.2, 0.9, 0))
  on.exit(par(old), add = TRUE)
  lim <- range(c(x[keep], y[keep]), na.rm = TRUE)
  if (diff(lim) <= 0 || !all(is.finite(lim))) lim <- lim + c(-0.5, 0.5)
  if (nrow(label_data)) lim <- lim + c(-0.16, 0.16) * diff(lim)
  display_point_size <- qc_plot_point_cex(point_size)
  plot(x[keep], y[keep], pch = 16, cex = display_point_size, col = adjustcolor("grey70", alpha.f = 0.7),
       xlim = lim, ylim = lim,
       xlab = "Control abundance, log10(CPM + 1)",
       ylab = "Treated abundance, log10(CPM + 1)",
       main = if (identical(mode, "individual") && nzchar(pair)) paste("Condition-control log-log scatter —", pair) else "Condition-control log-log scatter — combined across pools",
       cex.main = scales$title, cex.lab = scales$axis, cex.axis = scales$sample)
  if (isTRUE(grid)) grid(col = "grey85", lty = 1)
  abline(0, 1, lwd = qc_plot_lwd, lty = 2, col = adjustcolor("black", alpha.f = 0.45))
  if (length(highlight)) {
    depleted_col <- "#2f6fb5"; enriched_col <- "#c83f3f"
    hit_cols <- if (identical(direction, "depleted")) rep(depleted_col, length(highlight)) else if (identical(direction, "enriched")) rep(enriched_col, length(highlight)) else ifelse(qc_plot_numeric(ranked$log2fc[match(data$feature_id[highlight], ranked$feature_id)]) < 0, depleted_col, enriched_col)
    points(x[highlight], y[highlight], pch = 21, bg = hit_cols, col = "black", cex = display_point_size + qc_plot_point_cex(0.55), lwd = qc_plot_lwd)
  }
  if (nrow(label_data)) {
    label_data <- fitness_ma_place_labels(label_data, cex = scales$key)
    segments(label_data$x, label_data$y, label_data$x_label, label_data$y_label,
             col = adjustcolor("#68727d", alpha.f = 0.65), lwd = qc_plot_lwd * 0.5)
    text(label_data$x_label, label_data$y_label, labels = label_data$display_label,
         cex = scales$key, font = 2, xpd = NA,
         adj = ifelse(label_data$side > 0, 0, 1))
  }
  if (length(highlight)) {
    present <- c(depleted = any(ranked$log2fc[highlighted_ranked] < 0), enriched = any(ranked$log2fc[highlighted_ranked] > 0))
    legend("topleft", legend = c(depleted = "Depleted selected hit", enriched = "Enriched selected hit")[present],
           col = c(depleted = "#2f6fb5", enriched = "#c83f3f")[present], pch = 21,
           pt.bg = c(depleted = "#2f6fb5", enriched = "#c83f3f")[present], bty = "n", cex = scales$key)
  }
  note_line <- 5.75
  mtext("Below diagonal = depleted in treatment; above diagonal = enriched in treatment.", side = 1, line = note_line, cex = scales$key, col = "#666666")
  if (length(highlight) > 10L && nrow(label_data)) mtext("More than 10 labels may overlap; reduce Number of hits for cleaner labels.", side = 1, line = note_line + 0.9, cex = scales$key, col = "#666666")
}

fitness_condition_control_plot_file <- function(result) {
  if (is.null(result) || !nzchar(result$output_dir %||% "")) return("")
  file.path(result$output_dir, "plots", "condition_vs_control_loglog_scatter.png")
}

save_fitness_condition_control_scatter <- function(result, path = fitness_condition_control_plot_file(result),
                                                   width = 3000L, height = 2100L, res = 300L) {
  if (is.null(result) || !nzchar(path %||% "")) return("")
  data <- fitness_condition_control_plot_data(result)
  if (!nrow(data)) return("")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_fitness_condition_control_scatter(result)
  path
}
