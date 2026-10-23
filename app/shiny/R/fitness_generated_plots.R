fitness_generated_plot_inventory <- function(project) {
  root <- file.path(project$project_root, "treated_vs_parent")
  files <- if (dir.exists(root)) list.files(root, pattern = "\\.png$", recursive = TRUE, full.names = TRUE) else character()
  files <- files[grepl("/plots/", files) & file.exists(files)]
  # MA views are rendered from the selected result and current app controls.
  # Exclude the entire static MA family, regardless of strain/background name.
  # Keep the saved images on disk for provenance and ordinary file exports.
  files <- files[!grepl("^MA([_.-]|$)", basename(files), ignore.case = TRUE)]
  files <- files[!basename(files) %in% c(
    "top_features_log2FC_heatmap.png",
    "control_control_z_histogram.png",
    "library_sizes_feature_reads.png",
    "ranked_mean_log2FC.png",
    "mean_log2FC_distribution.png"
  )]
  if (!length(files)) return(data.frame())
  project_root <- normalizePath(project$project_root, winslash = "/", mustWork = FALSE)
  data.frame(
    file = files,
    filename = basename(files),
    title = tools::file_path_sans_ext(gsub("_", " ", basename(files))),
    served_url = vapply(files, function(path) {
      path <- normalizePath(path, winslash = "/", mustWork = FALSE)
      rel <- substring(path, nchar(project_root) + 2L)
      paste0("ytab-project-output/", diagnostic_encode_relative_path(rel), fitness_plot_cache_token(path))
    }, ""),
    stringsAsFactors = FALSE
  )
}

fitness_mean_log2fc_data <- function(result) {
  if (is.null(result) || !nzchar(result$table %||% "") || !file.exists(result$table)) return(data.frame())
  data <- tryCatch(read.csv(result$table, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (!nrow(data) || !"mean_log2FC" %in% names(data)) return(data.frame())
  if (!"feature_id" %in% names(data)) data$feature_id <- seq_len(nrow(data))
  data$feature_id <- as.character(data$feature_id)
  # The effect each call is based on (contract column call_log2FC); results made before the contract have only the pair mean.
  data$mean_log2FC <- qc_plot_numeric(if ("call_log2FC" %in% names(data)) data$call_log2FC else data$mean_log2FC)
  data <- data[is.finite(data$mean_log2FC), , drop = FALSE]
  if (!nrow(data)) return(data.frame())
  data <- data[order(data$mean_log2FC, data$feature_id), , drop = FALSE]
  data$rank <- seq_len(nrow(data))
  data
}

plot_fitness_ranked_mean_log2fc <- function(result) {
  data <- fitness_mean_log2fc_data(result)
  if (!nrow(data)) return(qc_plot_empty("log2FC values are unavailable for this fitness result."))
  scales <- qc_plot_text_sizes()
  old <- qc_plot_par(mar = c(5.6, 6.0, 4.3, 2.2), mgp = c(3.2, 0.9, 0))
  on.exit(par(old), add = TRUE)
  x <- data$rank
  y <- data$mean_log2FC
  ylim <- range(c(y, -1, 1), na.rm = TRUE)
  pad <- max(0.1, diff(ylim) * 0.06)
  ylim <- ylim + c(-pad, pad)
  plot(
    x, y,
    pch = 16,
    cex = qc_plot_point_cex(0.65),
    col = adjustcolor("#3f566c", alpha.f = 0.48),
    xlab = "Genes ranked by log2FC",
    ylab = "log2FC (treated / control)",
    main = "Ranked treated-versus-control log2FC",
    ylim = ylim,
    cex.main = scales$title,
    cex.lab = scales$axis,
    cex.axis = scales$sample
  )
  if (qc_plot_grid_enabled()) grid(col = "#e3e9ee", lty = 3, lwd = 1)
  abline(h = c(-1, 1), lty = 2, lwd = qc_plot_lwd, col = adjustcolor("black", alpha.f = 0.55))
  text(max(x, na.rm = TRUE), 1, "+1", pos = 3, xpd = NA, cex = scales$key, font = 2, col = "#666666")
  text(max(x, na.rm = TRUE), -1, "-1", pos = 1, xpd = NA, cex = scales$key, font = 2, col = "#666666")
  ytab_base_axes()
  invisible(data)
}

plot_fitness_mean_log2fc_distribution <- function(result) {
  data <- fitness_mean_log2fc_data(result)
  if (!nrow(data)) return(qc_plot_empty("log2FC values are unavailable for this fitness result."))
  scales <- qc_plot_text_sizes()
  old <- qc_plot_par(mar = c(5.6, 6.0, 4.3, 2.2), mgp = c(3.2, 0.9, 0))
  on.exit(par(old), add = TRUE)
  h <- hist(
    data$mean_log2FC,
    breaks = 80,
    plot = FALSE
  )
  border_par <- par(lwd = ytab_base_bar_lwd())  # rect() borders take par("lwd"), not the lwd argument
  plot(
    h,
    col = qc_plot_fill,
    border = "black",
    lwd = ytab_base_bar_lwd(),
    main = "Distribution of mean treated-versus-control log2FC",
    xlab = "log2FC (treated / control)",
    ylab = "Number of genes/features",
    cex.main = scales$title,
    cex.lab = scales$axis,
    cex.axis = scales$sample
  )
  par(border_par)
  if (qc_plot_grid_enabled()) abline(h = pretty(par("usr")[3:4]), col = "#e3e9ee", lty = 3, lwd = 1)
  abline(v = c(-1, 1), lty = 2, lwd = qc_plot_lwd, col = adjustcolor("black", alpha.f = 0.55))
  text(1, par("usr")[[4L]], "+1", pos = 1, xpd = NA, cex = scales$key, font = 2, col = "#666666")
  text(-1, par("usr")[[4L]], "-1", pos = 1, xpd = NA, cex = scales$key, font = 2, col = "#666666")
  ytab_base_axes()
  invisible(data)
}

fitness_library_size_data <- function(result) {
  path <- file.path(result$output_dir %||% "", "tables", "library_sizes.feature_reads.csv")
  if (!file.exists(path)) return(data.frame())
  data <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  required <- c("sample", "background", "total_feature_reads")
  if (!nrow(data) || !all(required %in% names(data))) return(data.frame())
  data$sample <- as.character(data$sample)
  data$background <- as.character(data$background)
  data$total_feature_reads <- qc_plot_numeric(data$total_feature_reads)
  if (!"condition" %in% names(data)) data$condition <- NA_character_
  if (!"pool" %in% names(data)) data$pool <- NA
  condition <- tolower(as.character(data$condition))
  sample <- tolower(data$sample)
  # Role from the sample sheet's condition, else generic words in the sample name (no project-specific names).
  role <- ifelse(condition %in% c("treated", "treatment") | (is.na(condition) | !nzchar(condition)) & grepl("treated", sample), "Treated",
                 ifelse(condition %in% c("parent", "control", "mock") | grepl("mock|control|parent", sample), "Parent/mock", "Other"))
  data$display_condition <- role
  data <- data[is.finite(data$total_feature_reads) & nzchar(data$sample), , drop = FALSE]
  pool_num <- suppressWarnings(as.numeric(as.character(data$pool)))
  data$.pool_order <- ifelse(is.na(pool_num), seq_len(nrow(data)), pool_num)
  data[order(data$background, data$.pool_order, data$display_condition, data$sample), , drop = FALSE]
}

fitness_library_size_scope_choices <- function(result) {
  data <- fitness_library_size_data(result)
  if (!nrow(data)) return(list())
  backgrounds <- sort(unique(as.character(data$background[nzchar(as.character(data$background))])))
  as.list(c("Combined backgrounds" = "combined", stats::setNames(backgrounds, backgrounds)))
}

plot_fitness_library_sizes_feature_reads <- function(result, scope = "combined") {
  data <- fitness_library_size_data(result)
  if (!nrow(data)) return(qc_plot_empty("Feature-level library-size data are unavailable for this fitness result."))
  scope <- as.character(scope %||% "combined")
  if (!identical(scope, "combined")) data <- data[as.character(data$background) == scope, , drop = FALSE]
  if (!nrow(data)) return(qc_plot_empty("No feature-level library-size values match the selected background."))
  data <- data[order(data$background, data$.pool_order, data$display_condition, data$sample), , drop = FALSE]
  labels <- qc_plot_display_labels(data$sample, dirname(dirname(result$output_dir)))
  values <- data$total_feature_reads
  groups <- unique(data$display_condition)
  palette <- c("Parent/mock" = "#4b728f", "Treated" = "#b45f45", "Other" = "grey55")
  missing_groups <- setdiff(groups, names(palette))
  if (length(missing_groups)) palette <- c(palette, stats::setNames(qc_plot_palette(length(missing_groups)), missing_groups))
  fills <- unname(palette[data$display_condition])
  horizontal <- qc_plot_bar_horizontal(labels)
  old <- par(no.readonly = TRUE)
  on.exit({ qc_plot_reset_layout(); par(old) }, add = TRUE)
  qc_plot_begin_key_layout(0.18)
  plot_margins <- qc_plot_bar_margins(labels, horizontal, top = 4.3, right = if (horizontal) 5.2 else 3.2)
  if (horizontal) plot_margins[[1L]] <- max(plot_margins[[1L]], 7.0)
  else {
    plot_margins[[1L]] <- max(plot_margins[[1L]], 13.5)
    plot_margins[[2L]] <- max(plot_margins[[2L]], 9.2)
  }
  qc_plot_par(mar = plot_margins)
  limit <- max(values, na.rm = TRUE) * 1.12
  border_par <- par(lwd = ytab_base_bar_lwd())  # rect() borders take par("lwd"), not the lwd argument
  at <- barplot(
    values,
    names.arg = rep("", length(values)),
    axes = FALSE,
    horiz = horizontal,
    xlim = if (horizontal) c(0, limit) else NULL,
    ylim = if (horizontal) NULL else c(0, limit),
    col = fills,
    border = "black",
    lwd = ytab_base_bar_lwd(),
    main = if (identical(scope, "combined")) "Feature-level library sizes — combined backgrounds" else paste("Feature-level library sizes —", scope),
    xlab = if (horizontal) "Total feature reads" else "",
    ylab = ""
  )
  par(border_par)
  if (!horizontal) title(ylab = "Total feature reads", line = 7.0)
  qc_plot_add_grid(horizontal)
  if (horizontal) {
    axis(1, labels = format(pretty(c(0, limit)), big.mark = ",", scientific = FALSE), at = pretty(c(0, limit)), lwd = 0, lwd.ticks = ytab_mm_lwd(get0("YTAB_TICK_MM", ifnotfound = .6)))
    qc_plot_draw_horizontal_labels(at, labels)
  } else {
    axis(2, labels = format(pretty(c(0, limit)), big.mark = ",", scientific = FALSE), at = pretty(c(0, limit)), las = 1, lwd = 0, lwd.ticks = ytab_mm_lwd(get0("YTAB_TICK_MM", ifnotfound = .6)))
    qc_plot_draw_vertical_labels(at, labels)
  }
  if (qc_plot_show_value_labels()) {
    value_labels <- format(round(values), big.mark = ",", scientific = FALSE)
    if (horizontal) text(values, at, labels = value_labels, pos = 4, xpd = NA, font = 2, cex = qc_plot_value_cex(0.72))
    else text(at, values + max(values, na.rm = TRUE) * 0.018, labels = value_labels, pos = 3, xpd = NA, font = 2, cex = qc_plot_value_cex(0.72), srt = 30)
  }
  ytab_base_axes(ticks = integer())
  par(mar = c(0.2, 0.5, 0.2, 0.5))
  plot.new()
  text(0.5, 0.82, paste0("Reads per library before normalization", if (nzchar(result$library_sizes_source %||% "")) paste0(" (", result$library_sizes_source, ")") else "", "."),
       cex = qc_plot_key_cex(0.82), col = "#666666", font = 2)
  legend(
    "bottom",
    legend = groups,
    fill = unname(palette[groups]),
    border = qc_plot_border,
    title = "Sample group",
    ncol = length(groups),
    bty = "n",
    text.font = 2,
    cex = qc_plot_key_cex(0.92),
    x.intersp = 0.8,
    y.intersp = 0.95
  )
  invisible(data)
}

fitness_control_z_saved_data <- function(result) {
  path <- file.path(result$output_dir %||% "", "tables", "control_control_z_scores.csv")
  if (!file.exists(path)) return(data.frame())
  data <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (!nrow(data) || !all(c("background", "z") %in% names(data))) return(data.frame())
  data$background <- as.character(data$background)
  data$z <- qc_plot_numeric(data$z)
  data[is.finite(data$z) & nzchar(data$background), , drop = FALSE]
}

# The pipeline's own control-control z (the noise model behind the call threshold); never refitted in the app.
fitness_control_z_data <- function(result) {
  if (is.null(result) || !nzchar(result$output_dir %||% "")) return(data.frame())
  fitness_control_z_saved_data(result)
}

fitness_control_z_scope_choices <- function(result) {
  data <- fitness_control_z_data(result)
  if (!nrow(data)) return(list())
  backgrounds <- sort(unique(as.character(data$background[nzchar(as.character(data$background))])))
  choices <- c("Combined backgrounds" = "combined", stats::setNames(backgrounds, backgrounds))
  as.list(choices)
}

plot_fitness_control_control_z_histogram <- function(result, scope = "combined") {
  data <- fitness_control_z_data(result)
  if (!nrow(data)) return(qc_plot_empty("Control-control z-score values are unavailable for this fitness result."))
  scope <- as.character(scope %||% "combined")
  if (!identical(scope, "combined")) data <- data[as.character(data$background) == scope, , drop = FALSE]
  if (!nrow(data)) return(qc_plot_empty("No control-control z-score values match the selected background."))
  z <- qc_plot_numeric(data$z)
  z <- z[is.finite(z)]
  if (!length(z)) return(qc_plot_empty("No finite control-control z-score values are available."))
  scales <- qc_plot_text_sizes()
  old <- qc_plot_par(mar = c(5.5, 5.8, 4.4, 1.8))
  on.exit(par(old), add = TRUE)
  border_par <- par(lwd = ytab_base_bar_lwd())  # rect() borders take par("lwd"), not the lwd argument
  hist(
    z,
    breaks = 100,
    col = qc_plot_fill,
    border = "black",
    lwd = ytab_base_bar_lwd(),
    main = if (identical(scope, "combined")) "Control-control z histogram — combined backgrounds" else paste("Control-control z histogram —", scope),
    xlab = "Control-control z-score",
    ylab = "Number of genes/features",
    cex.main = scales$title,
    cex.lab = scales$axis,
    cex.axis = scales$sample
  )
  par(border_par)
  if (qc_plot_grid_enabled()) abline(h = pretty(par("usr")[3:4]), col = "#e3e9ee", lty = 3, lwd = 1)
  ytab_base_axes()
  if (identical(scope, "combined")) {
    mtext("Combined view pools control-control z-scores across backgrounds; background-specific views remain available for model diagnostics.", side = 1, line = 4.2, cex = scales$key, col = "#666666")
  } else {
    mtext("Background-specific view reflects the parent-parent local noise model used for this background.", side = 1, line = 4.2, cex = scales$key, col = "#666666")
  }
  invisible(data)
}

fitness_plot_cache_token <- function(path) {
  info <- file.info(path)
  if (!nrow(info) || is.na(info$mtime[[1]])) return("")
  paste0("?v=", as.integer(as.POSIXct(info$mtime[[1]])), "-", as.integer(info$size[[1]] %||% 0L))
}

fitness_combined_ma_plot_file <- function(result) {
  if (is.null(result) || !nzchar(result$output_dir %||% "")) return("")
  candidates <- file.path(result$output_dir, "plots",
                          c("MA_treated_vs_parent_combined.png", "MA_treated_vs_parent_all_pools.png"))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) candidates[[1]] else ""
}

fitness_generated_plot_inventory_for_result <- function(project, result) {
  plots <- fitness_generated_plot_inventory(project)
  if (!nrow(plots) || is.null(result) || !nzchar(result$output_dir %||% "")) return(plots)
  ma_candidates <- file.path(result$output_dir, "plots",
                             c("MA_treated_vs_parent_combined.png", "MA_treated_vs_parent_all_pools.png"))
  dedicated_candidates <- c(ma_candidates, fitness_condition_control_plot_file(result))
  dedicated_candidates <- normalizePath(dedicated_candidates[file.exists(dedicated_candidates)], winslash = "/", mustWork = FALSE)
  if (!length(dedicated_candidates)) return(plots)
  plots[!normalizePath(plots$file, winslash = "/", mustWork = FALSE) %in% dedicated_candidates, , drop = FALSE]
}

fitness_project_served_plot <- function(path, project) {
  if (!nzchar(path %||% "") || !file.exists(path)) return("")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  project_root <- normalizePath(project$project_root, winslash = "/", mustWork = FALSE)
  if (!startsWith(path, paste0(project_root, "/"))) return("")
  rel <- substring(path, nchar(project_root) + 2L)
  paste0("ytab-project-output/", diagnostic_encode_relative_path(rel), fitness_plot_cache_token(path))
}

fitness_combined_ma_plot_card <- function(result, project) {
  path <- fitness_combined_ma_plot_file(result)
  served_url <- fitness_project_served_plot(path, project)
  if (!nzchar(served_url)) {
    return(tags$p(class = "text-muted", "The combined treated-versus-parent MA plot is not available for this selected fitness result."))
  }
  tags$article(
    class = "ytab-static-image-card ytab-release-card",
    tags$div(class = "ytab-plot-card-header",
             tags$div(tags$h4("Generated: Combined treated-versus-parent MA plot"),
                      tags$span(class = "ytab-status-badge", "Static generated image")),
             tags$a(class = "btn btn-secondary btn-sm", href = served_url, target = "_blank", "Open original")),
    tags$div(
      class = "ytab-diagnostic-preview",
      tags$img(src = served_url, alt = basename(path), loading = "lazy",
               style = "width:100%;max-height:620px;object-fit:contain")
    ),
    tags$details(tags$summary("Show filename"), tags$code(relative_project_path(path, project$project_root)))
  )
}

fitness_condition_control_plot_card <- function(result, project) {
  path <- fitness_condition_control_plot_file(result)
  if (!nzchar(path) || !file.exists(path)) {
    path <- tryCatch(save_fitness_condition_control_scatter(result), error = function(e) "")
  }
  served_url <- fitness_project_served_plot(path, project)
  if (!nzchar(served_url)) {
    return(tags$p(class = "text-muted", "The condition-versus-control log-log scatter plot is not available for this selected fitness result."))
  }
  tags$article(
    class = "ytab-static-image-card ytab-release-card",
    tags$div(class = "ytab-plot-card-header",
             tags$div(tags$h4("Generated: Condition-versus-control log-log scatter"),
                      tags$span(class = "ytab-status-badge", "Static generated image")),
             tags$a(class = "btn btn-secondary btn-sm", href = served_url, target = "_blank", "Open original")),
    tags$div(
      class = "ytab-diagnostic-preview",
      tags$img(src = served_url, alt = basename(path), loading = "lazy",
               style = "width:100%;max-height:620px;object-fit:contain")
    ),
    tags$details(tags$summary("Show filename"), tags$code(relative_project_path(path, project$project_root)))
  )
}

fitness_generated_plot_cards <- function(project) {
  plots <- fitness_generated_plot_inventory(project)
  if (!nrow(plots)) return(tags$p(class = "text-muted", "No treated-vs-parent plot files are available."))
  plots <- plots[order(plots$filename), , drop = FALSE]
  tags$div(class = "ytab-plot-grid", lapply(seq_len(nrow(plots)), function(i) {
    tags$article(
      class = "ytab-static-image-card ytab-release-card",
      title = plots$filename[[i]],
      tags$div(class = "ytab-plot-card-header",
               tags$div(tags$h4(paste("Generated:", plots$title[[i]])),
                        tags$span(class = "ytab-status-badge", "Static generated image")),
               tags$a(class = "btn btn-secondary btn-sm", href = plots$served_url[[i]], target = "_blank", "Open original")),
      tags$div(class = "ytab-diagnostic-preview",
        tags$img(src = plots$served_url[[i]], alt = plots$filename[[i]],
                 loading = "lazy", style = "width:100%;max-height:520px;object-fit:contain")
      ),
      tags$details(tags$summary("Show filename"), tags$code(plots$filename[[i]]))
    )
  }))
}
