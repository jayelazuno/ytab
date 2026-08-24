qc_library_midlc_plot_data <- function(project) {
  roots <- c(file.path(project$project_root, "library_diagnostics", "runs", "all_samples"),
             file.path(project$project_root, "library_diagnostics"))
  root_files <- lapply(roots, function(root) {
    if (dir.exists(root)) list.files(root, pattern = "\\.midlc\\.csv$", recursive = TRUE, full.names = TRUE) else character()
  })
  available <- which(vapply(root_files, length, integer(1)) > 0L)
  files <- if (length(available)) root_files[[available[[1]]]] else character()
  rows <- lapply(files, function(file) {
    data <- tryCatch(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(data) || !"Reads Sampled" %in% names(data)) return(NULL)
    sample <- sub("\\.midlc\\.csv$", "", basename(file))
    values <- data[, setdiff(names(data), "Reads Sampled"), drop = FALSE]
    data.frame(
      sample = sample,
      reads_sampled = qc_plot_numeric(data[["Reads Sampled"]]),
      unique_sites = rowMeans(as.data.frame(lapply(values, qc_plot_numeric)), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data) || !nrow(data)) return(data.frame())
  data$sample <- factor(data$sample, levels = qc_plot_sample_order(as.character(data$sample)))
  data[order(data$sample, data$reads_sampled), , drop = FALSE]
}

qc_library_project_samples <- function(project) {
  if (is.data.frame(project$samples) && nrow(project$samples)) return(project$samples)
  path <- as.character(project$sample_sheet %||% "")
  if (nzchar(path)) {
    candidates <- c(path, file.path(project$repo_root %||% "", path),
                    file.path(project$project_root %||% "", "config", basename(path)))
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates)) {
      return(tryCatch(read.csv(candidates[[1]], stringsAsFactors = FALSE, check.names = FALSE),
                      error = function(e) data.frame()))
    }
  }
  data.frame()
}

qc_library_role_from_values <- function(values) {
  values <- tolower(paste(stats::na.omit(as.character(values)), collapse = " "))
  if (!nzchar(values)) return(NA_character_)
  if (grepl("\\btreated\\b|zn-treated|h2o2-treated|1_5mm|h2o2", values)) return("treated")
  if (grepl("\\bcontrol\\b|\\bmock\\b|\\bparent\\b|classifier_control", values)) return("control")
  NA_character_
}

qc_library_sample_roles <- function(project, samples = NULL) {
  if (is.null(samples)) {
    samples <- unique(c(as.character(qc_library_summary_plot_data(project)$sample),
                        as.character(qc_library_midlc_plot_data(project)$sample)))
  }
  samples <- qc_plot_sample_order(unique(as.character(samples)))
  metadata <- qc_library_project_samples(project)
  role_cols <- c("library_role", "control_or_treated", "fitness_role", "classifier_role",
                 "treatment", "condition", "condition_label", "treatment_label",
                 "guessed_condition")
  rows <- lapply(samples, function(sample) {
    row <- if (is.data.frame(metadata) && "sample" %in% names(metadata)) metadata[as.character(metadata$sample) == sample, , drop = FALSE] else data.frame()
    role <- NA_character_
    condition_label <- ""
    pool <- qc_plot_pool(sample)
    if (nrow(row)) {
      for (col in intersect(role_cols, names(row))) {
        role <- qc_library_role_from_values(row[[col]])
        if (!is.na(role)) break
      }
      label_cols <- intersect(c("condition_label", "treatment_label", "condition", "treatment"), names(row))
      labels <- unique(as.character(unlist(row[, label_cols, drop = FALSE])))
      labels <- labels[nzchar(labels) & !is.na(labels)]
      if (length(labels)) condition_label <- labels[[1]]
      if ("pool" %in% names(row) && nzchar(as.character(row$pool[[1]]))) pool <- as.character(row$pool[[1]])
      if ("pool_id" %in% names(row) && nzchar(as.character(row$pool_id[[1]]))) pool <- as.character(row$pool_id[[1]])
    }
    if (is.na(role)) {
      sample_lower <- tolower(sample)
      if (grepl("treated|h2o2|zn-treated|1_5mm", sample_lower)) role <- "treated"
      else if (grepl("mock|control|parent", sample_lower)) role <- "control"
      else role <- "other"
    }
    display_group <- if (identical(role, "control")) {
      if (grepl("mock", tolower(paste(condition_label, sample)))) "Mock control" else "Parent/control"
    } else if (identical(role, "treated")) {
      if (grepl("1\\.5|1_5|zn", tolower(paste(condition_label, sample)))) "1.5 mM Zn-treated"
      else if (grepl("h2o2", tolower(paste(condition_label, sample)))) "H2O2-treated"
      else "Treated"
    } else "Other"
    data.frame(sample = sample, role = role, display_group = display_group,
               pool = pool, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

qc_library_group_choices <- function(project) {
  roles <- qc_library_sample_roles(project)
  has_control <- any(roles$role == "control")
  has_treated <- any(roles$role == "treated")
  control_label <- unique(roles$display_group[roles$role == "control"])[[1]] %||% "Controls"
  treated_label <- unique(roles$display_group[roles$role == "treated"])[[1]] %||% "Treated"
  choices <- list("All eligible samples" = "all")
  if (has_control) choices[[paste0(control_label, " only")]] <- "control"
  if (has_treated) choices[[paste0(treated_label, " only")]] <- "treated"
  if (has_control && has_treated) {
    choices[[paste(control_label, treated_label, sep = " + ")]] <- "both"
    if (length(intersect(roles$pool[roles$role == "control"], roles$pool[roles$role == "treated"]))) {
      choices[[paste("Matched", control_label, "/", treated_label, "pairs by pool")]] <- "matched_pairs"
    }
  }
  choices
}

qc_library_color_choices <- function() {
  c("Sample group" = "group", "Pool" = "pool", "Sample" = "sample", "None / grayscale" = "none")
}

qc_library_group_color_map <- function() {
  c(control = "#2c7fb8", treated = "#d95f02", other = "grey45")
}

qc_library_group_color <- function(role) {
  colors <- qc_library_group_color_map()
  out <- unname(colors[as.character(role)])
  out[is.na(out)] <- unname(colors[["other"]])
  out
}

qc_library_plot_aesthetics <- function(visualization) {
  defaults <- list(
    use_lines = TRUE,
    use_points = FALSE,
    use_bars = FALSE,
    line_by_sample = TRUE,
    sample_linetype = FALSE,
    sample_shapes = FALSE,
    pair_connectors = FALSE,
    individual_alpha = 0.45,
    summary_line = TRUE,
    legend_mode = "group"
  )
  switch(
    as.character(visualization),
    midlc = modifyList(defaults, list(use_points = TRUE, individual_alpha = 0.7, summary_line = FALSE)),
    centromere = modifyList(defaults, list(use_points = FALSE, individual_alpha = 0.28, summary_line = TRUE)),
    metaplot = modifyList(defaults, list(use_points = FALSE, individual_alpha = 0.25, summary_line = TRUE)),
    sequence_bias = modifyList(defaults, list(use_lines = FALSE, use_points = FALSE, use_bars = TRUE,
                                              line_by_sample = FALSE, summary_line = FALSE,
                                              legend_mode = "group")),
    jackpot = modifyList(defaults, list(use_lines = FALSE, use_points = TRUE, use_bars = TRUE,
                                        line_by_sample = FALSE, summary_line = FALSE)),
    defaults
  )
}

qc_library_filter_data <- function(project, data, group = "all") {
  if (!is.data.frame(data) || !nrow(data) || !"sample" %in% names(data)) return(data)
  roles <- qc_library_sample_roles(project, unique(as.character(data$sample)))
  if (identical(group, "control")) roles <- roles[roles$role == "control", , drop = FALSE]
  else if (identical(group, "treated")) roles <- roles[roles$role == "treated", , drop = FALSE]
  else if (identical(group, "both")) roles <- roles[roles$role %in% c("control", "treated"), , drop = FALSE]
  else if (identical(group, "matched_pairs")) {
    pools <- intersect(roles$pool[roles$role == "control"], roles$pool[roles$role == "treated"])
    roles <- roles[roles$pool %in% pools & roles$role %in% c("control", "treated"), , drop = FALSE]
  }
  roles$order <- seq_len(nrow(roles))
  out <- merge(data, roles, by = "sample", all.x = FALSE, all.y = FALSE, sort = FALSE)
  out <- unique(out)
  out <- out[order(out$order, out$sample), , drop = FALSE]
  out$sample <- factor(out$sample, levels = unique(roles$sample))
  out$role <- factor(out$role, levels = intersect(c("control", "treated", "other"), unique(as.character(out$role))))
  out$display_group <- factor(out$display_group, levels = unique(as.character(out$display_group)))
  out <- droplevels(out)
  out
}

qc_library_sample_style <- function(data, color_by = "group") {
  samples <- unique(as.character(data$sample))
  roles <- unique(data[, intersect(c("sample", "role", "display_group", "pool"), names(data)), drop = FALSE])
  roles <- roles[match(samples, roles$sample), , drop = FALSE]
  if (identical(color_by, "none")) {
    colors <- rep("grey35", length(samples))
    key <- data.frame(label = "All samples", color = "grey35", pch = 16, lty = 1, stringsAsFactors = FALSE)
  } else if (identical(color_by, "sample")) {
    colors <- qc_plot_palette(length(samples))
    key <- data.frame(label = qc_plot_display_labels(samples), color = colors, pch = 16, lty = 1, stringsAsFactors = FALSE)
  } else if (identical(color_by, "pool")) {
    pools <- unique(as.character(roles$pool))
    pool_colors <- stats::setNames(qc_plot_palette(length(pools)), pools)
    colors <- unname(pool_colors[as.character(roles$pool)])
    key <- data.frame(label = paste("Pool", pools), color = unname(pool_colors), pch = 16, lty = 1, stringsAsFactors = FALSE)
  } else {
    group_colors <- qc_library_group_color_map()
    colors <- qc_library_group_color(roles$role)
    groups <- unique(roles[, c("role", "display_group"), drop = FALSE])
    key <- data.frame(label = as.character(groups$display_group),
                      color = qc_library_group_color(groups$role),
                      pch = 16, lty = 1, stringsAsFactors = FALSE)
  }
  pch <- rep(c(16, 17, 15, 18), length.out = length(samples))
  lty <- rep(c(1, 2, 3, 4), length.out = length(samples))
  list(samples = samples, colors = colors, pch = pch, lty = lty, key = key)
}

qc_library_curve_color <- function(row, color_by = "group", style = NULL) {
  if (identical(color_by, "group")) return(qc_library_group_color(row$role))
  if (identical(color_by, "none")) return(rep("grey35", nrow(row)))
  if (!is.null(style)) {
    return(style$colors[match(as.character(row$sample), style$samples)])
  }
  qc_library_group_color(row$role)
}

qc_library_plot_cache_key <- function(project, visualization, sample_group = "all",
                                      color_by = "group", metaplot_panel = "",
                                      selected_samples = character(),
                                      diagnostics_result_id = "") {
  paste(
    as.character(project$project_id %||% ""),
    as.character(diagnostics_result_id %||% ""),
    as.character(visualization %||% ""),
    as.character(sample_group %||% "all"),
    as.character(color_by %||% "group"),
    paste(sort(as.character(selected_samples %||% character())), collapse = ","),
    as.character(metaplot_panel %||% ""),
    as.character(getOption("ytab.plot.style", "clean")),
    as.character(getOption("ytab.plot.text_size", "medium")),
    as.character(getOption("ytab.plot.label_mode", "full")),
    as.character(getOption("ytab.plot.grid", "show")),
    as.character(getOption("ytab.plot.bar_orientation", "vertical")),
    as.character(getOption("ytab.plot.show_value_labels", TRUE)),
    sep = "|"
  )
}

qc_library_sample_key_ui <- function(project, data, group = "all") {
  if (!is.data.frame(data) || !nrow(data) || !"sample" %in% names(data)) return(NULL)
  roles <- unique(qc_library_filter_data(project, data.frame(sample = unique(as.character(data$sample))), group))
  if (!nrow(roles)) return(NULL)
  tags$details(
    class = "ytab-more-options",
    tags$summary("Sample key"),
    tags$table(
      class = "table table-sm",
      tags$thead(tags$tr(tags$th("Sample"), tags$th("Group"), tags$th("Pool"))),
      tags$tbody(lapply(seq_len(nrow(roles)), function(i) {
        tags$tr(tags$td(roles$sample[[i]]), tags$td(roles$display_group[[i]]), tags$td(roles$pool[[i]]))
      }))
    )
  )
}

qc_library_plot_provenance_ui <- function(source_label, source_files = character()) {
  tags$details(
    class = "ytab-more-options",
    tags$summary("Plot provenance"),
    tags$p(class = "text-muted", source_label),
    if (length(source_files)) tags$ul(lapply(head(source_files, 8), function(path) tags$li(tags$code(basename(path)))))
  )
}

plot_qc_library_midlc <- function(project, group = "all", color_by = "group") {
  data <- qc_library_midlc_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("MidLC curves are not available."))
  data <- qc_library_filter_data(project, data, group)
  if (!nrow(data)) return(qc_plot_empty("No MidLC curves are available for the selected sample group."))
  old <- par(no.readonly = TRUE)
  on.exit({ qc_plot_reset_layout(); par(old) }, add = TRUE)
  qc_plot_begin_key_layout(0.16)
  qc_plot_par(mar = c(5, 5, 3.5, 1))
  plot(NA, xlim = range(data$reads_sampled, na.rm = TRUE),
       ylim = range(data$unique_sites, na.rm = TRUE), log = "x",
       xlab = "Reads sampled", ylab = "Unique insertion sites",
       main = "MidLC saturation")
  style <- qc_library_sample_style(data, color_by)
  aesthetic <- qc_library_plot_aesthetics("midlc")
  if (qc_plot_grid_enabled()) {
    abline(h = pretty(range(data$unique_sites, na.rm = TRUE)), col = "#e3e9ee", lty = 3, lwd = 1)
    abline(v = pretty(range(data$reads_sampled, na.rm = TRUE)), col = "#e3e9ee", lty = 3, lwd = 1)
  }
  for (i in seq_along(style$samples)) {
    row <- data[as.character(data$sample) == style$samples[[i]], , drop = FALSE]
    row <- row[order(row$reads_sampled), , drop = FALSE]
    if (nrow(row)) lines(row$reads_sampled, row$unique_sites,
                         type = if (isTRUE(aesthetic$use_points)) "b" else "l",
                         pch = 16, col = style$colors[[i]],
                         lty = 1, lwd = qc_plot_lwd)
  }
  qc_plot_metric_key_row(style$key$label, col = style$key$color, pch = style$key$pch,
                         lty = style$key$lty, lwd = qc_plot_lwd,
                         ncol = min(3L, nrow(style$key)), title = "Key")
}

qc_library_summary_plot_data <- function(project) {
  candidates <- c(
    file.path(project$project_root, "library_diagnostics", "runs", "all_samples", "library_diagnostics.summary.csv"),
    file.path(project$project_root, "library_diagnostics", "library_diagnostics.summary.csv")
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) return(data.frame())
  data <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (nrow(data) && "sample" %in% names(data)) data$sample <- factor(data$sample, levels = qc_plot_sample_order(as.character(data$sample)))
  data
}

plot_qc_library_jackpot_depth <- function(project, group = "all", color_by = "group") {
  data <- qc_library_summary_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("Library diagnostic summary is not available."))
  data <- qc_library_filter_data(project, data, group)
  if (!nrow(data)) return(qc_plot_empty("No library diagnostic summary is available for the selected sample group."))
  raw_labels <- as.character(data$sample)
  labels <- qc_plot_display_labels(raw_labels)
  jackpot <- qc_plot_numeric(qc_plot_column(data, c("jackpot_frac_reads", "jackpot_top_frac")))
  depth <- qc_plot_numeric(qc_plot_column(data, "depth_ratio_R_over_midlc"))
  horizontal <- qc_plot_bar_horizontal(labels)
  style <- qc_library_sample_style(data, color_by)
  bar_cols <- style$colors[match(as.character(data$sample), style$samples)]
  old <- par(no.readonly = TRUE)
  on.exit({ qc_plot_reset_layout(); par(old) }, add = TRUE)
  metric_labels <- c(style$key$label, "Triangle = depth / MidLC ratio, scaled to plot range")
  qc_plot_begin_key_layout(0.24)
  qc_plot_par(mar = qc_plot_bar_margins(labels, horizontal, top = 3.8, right = if (horizontal) 7 else 5.5))
  values <- jackpot * 100
  limit <- max(values, na.rm = TRUE)
  at <- barplot(values, names.arg = rep("", length(labels)), axes = FALSE,
                horiz = horizontal,
                xlim = if (horizontal) c(0, limit * 1.18) else NULL,
                ylim = if (horizontal) NULL else c(0, limit * 1.18),
                col = bar_cols, border = qc_plot_border, lwd = qc_plot_lwd,
                xlab = if (horizontal) "Jackpot reads (%)" else "",
                ylab = if (horizontal) "" else "Jackpot reads (%)",
                main = "Jackpots and library depth")
  qc_plot_add_grid(horizontal)
  if (horizontal) axis(1) else axis(2)
  if (horizontal) qc_plot_draw_horizontal_labels(at, labels) else qc_plot_draw_vertical_labels(at, labels)
  if (qc_plot_show_value_labels()) {
    text(if (horizontal) values else at, if (horizontal) at else values,
         labels = sprintf("%.1f%%", values),
         pos = if (horizontal) 4 else 3, xpd = NA, font = 2,
         cex = qc_plot_label_cex(0.82))
  }
  if (any(!is.na(depth))) {
    depth_scaled <- if (max(depth, na.rm = TRUE) > 0) depth / max(depth, na.rm = TRUE) * limit else depth
    points(if (horizontal) depth_scaled else at, if (horizontal) at else depth_scaled,
           pch = 17, col = "black", cex = qc_plot_label_cex(1.05), lwd = qc_plot_lwd)
  }
  qc_plot_metric_key_row(metric_labels,
                         fill = c(style$key$color, NA), col = c(style$key$color, "black"),
                         border = c(rep(qc_plot_border, nrow(style$key)), NA),
                         pch = c(rep(NA, nrow(style$key)), 17), ncol = 1L,
                         title = "Metrics")
}

qc_library_sequence_bias_plot_data <- function(project) {
  roots <- c(file.path(project$project_root, "library_diagnostics", "runs", "all_samples"),
             file.path(project$project_root, "library_diagnostics"))
  root_files <- lapply(roots, function(root) {
    if (dir.exists(root)) list.files(root, pattern = "\\.seqbias_2_7\\.tsv$", recursive = TRUE, full.names = TRUE) else character()
  })
  available <- which(vapply(root_files, length, integer(1)) > 0L)
  files <- if (length(available)) root_files[[available[[1]]]] else character()
  rows <- lapply(files, function(file) {
    data <- tryCatch(read.delim(file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(data) || !"pair" %in% names(data) || !"enrichment" %in% names(data)) return(NULL)
    sample <- sub("\\.seqbias_2_7\\.tsv$", "", basename(file))
    data.frame(sample = sample, pair = data$pair,
               enrichment = qc_plot_numeric(data$enrichment), stringsAsFactors = FALSE)
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data)) data.frame() else data
}

plot_qc_library_sequence_bias <- function(project, group = "all") {
  data <- qc_library_sequence_bias_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("Sequence-bias tables are not available."))
  data <- qc_library_filter_data(project, data, group)
  if (!nrow(data)) return(qc_plot_empty("No sequence-bias tables are available for the selected sample group."))
  summarized <- aggregate(enrichment ~ role + display_group + pair, data, mean, na.rm = TRUE)
  pairs <- sort(unique(summarized$pair))
  group_table <- unique(summarized[, c("role", "display_group"), drop = FALSE])
  group_table <- group_table[order(match(as.character(group_table$role), c("control", "treated", "other")),
                                   as.character(group_table$display_group)), , drop = FALSE]
  conditions <- as.character(group_table$display_group)
  mat <- sapply(seq_len(nrow(group_table)), function(i) {
    row <- summarized[as.character(summarized$role) == as.character(group_table$role[[i]]) &
                        as.character(summarized$display_group) == as.character(group_table$display_group[[i]]), ]
    stats::setNames(row$enrichment, row$pair)[pairs]
  })
  colnames(mat) <- conditions
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = length(pairs), dimnames = list(pairs, conditions))
  display_pairs <- qc_plot_display_labels(pairs)
  horizontal <- qc_plot_bar_horizontal(display_pairs)
  old <- par(no.readonly = TRUE)
  on.exit({ qc_plot_reset_layout(); par(old) }, add = TRUE)
  qc_plot_begin_key_layout(0.15)
  qc_plot_par(mar = qc_plot_bar_margins(display_pairs, horizontal, top = 3.8, right = 2.6))
  group_cols <- qc_library_group_color(group_table$role)
  at <- barplot(t(mat), beside = TRUE, names.arg = rep("", length(pairs)), axes = FALSE,
          horiz = horizontal,
          col = group_cols,
          border = qc_plot_border, lwd = qc_plot_lwd,
          xlab = if (horizontal) "Mean enrichment" else "Dinucleotide (+2/+7)",
          ylab = if (horizontal) "Dinucleotide (+2/+7)" else "Mean enrichment",
          main = "Insertion-site sequence bias")
  qc_plot_add_grid(horizontal)
  if (horizontal) axis(1) else axis(2)
  label_at <- if (is.matrix(at)) colMeans(at) else at
  if (horizontal) qc_plot_draw_horizontal_labels(label_at, display_pairs) else qc_plot_draw_vertical_labels(label_at, display_pairs)
  if (horizontal) abline(v = 1, lty = 2, lwd = qc_plot_lwd, col = "black") else abline(h = 1, lty = 2, lwd = qc_plot_lwd, col = "black")
  if (qc_plot_show_value_labels()) {
    value_mat <- t(mat)
    text(if (horizontal) as.vector(value_mat) else as.vector(at),
         if (horizontal) as.vector(at) else as.vector(value_mat),
         labels = sprintf("%.2f", as.vector(value_mat)),
         pos = if (horizontal) 4 else 3, xpd = NA, font = 2,
         cex = qc_plot_label_cex(0.7))
  }
  qc_plot_metric_key_row(conditions,
                         fill = group_cols,
                         border = qc_plot_border, ncol = min(3L, length(conditions)))
}

qc_library_centromere_plot_data <- function(project) {
  roots <- c(file.path(project$project_root, "library_diagnostics", "runs", "all_samples"),
             file.path(project$project_root, "library_diagnostics"))
  root_files <- lapply(roots, function(root) {
    if (dir.exists(root)) list.files(root, pattern = "\\.centromere_bins\\.tsv$", recursive = TRUE, full.names = TRUE) else character()
  })
  available <- which(vapply(root_files, length, integer(1)) > 0L)
  files <- if (length(available)) root_files[[available[[1]]]] else character()
  rows <- lapply(files, function(file) {
    data <- tryCatch(read.delim(file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(data)) return(NULL)
    names(data)[names(data) == "bin_start"] <- "bin_start_bp"
    names(data)[names(data) == "mean_reads"] <- "mean_reads_across_arms"
    if (!all(c("bin_start_bp", "mean_reads_across_arms") %in% names(data))) return(NULL)
    sample <- sub("\\.centromere_bins\\.tsv$", "", basename(file))
    data.frame(sample = sample,
               distance_kb = qc_plot_numeric(data$bin_start_bp) / 1000,
               mean_reads = qc_plot_numeric(data$mean_reads_across_arms),
               stringsAsFactors = FALSE)
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data)) data.frame() else data
}

plot_qc_library_centromere_bias <- function(project, group = "all", color_by = "group") {
  data <- qc_library_centromere_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("Centromere-bias tables are not available."))
  data <- qc_library_filter_data(project, data, group)
  if (!nrow(data)) return(qc_plot_empty("No centromere-bias data are available for the selected sample group."))
  effective_color_by <- if (color_by %in% c("sample", "none")) color_by else "group"
  style <- qc_library_sample_style(data, effective_color_by)
  aesthetic <- qc_library_plot_aesthetics("centromere")
  old <- par(no.readonly = TRUE)
  on.exit({ qc_plot_reset_layout(); par(old) }, add = TRUE)
  qc_plot_begin_key_layout(0.16)
  qc_plot_par(mar = c(5, 5.5, 3.5, 1))
  plot(NA, xlim = range(data$distance_kb, na.rm = TRUE), ylim = range(data$mean_reads, na.rm = TRUE),
       xlab = "Distance from centromere (kb)", ylab = "Mean reads per bin",
       main = "Centromere bias")
  if (qc_plot_grid_enabled()) {
    abline(h = pretty(range(data$mean_reads, na.rm = TRUE)), col = "#e3e9ee", lty = 3)
    abline(v = pretty(range(data$distance_kb, na.rm = TRUE)), col = "#e3e9ee", lty = 3)
  }
  if (identical(effective_color_by, "sample")) {
    for (i in seq_along(style$samples)) {
      row <- data[as.character(data$sample) == style$samples[[i]], , drop = FALSE]
      row <- row[order(row$distance_kb), , drop = FALSE]
      if (nrow(row)) lines(row$distance_kb, row$mean_reads,
                           col = style$colors[[i]], lty = 1, lwd = qc_plot_lwd)
    }
  } else if (isTRUE(aesthetic$summary_line)) {
    summary <- aggregate(mean_reads ~ display_group + role + distance_kb, data, mean, na.rm = TRUE)
    for (group_label in unique(as.character(summary$display_group))) {
      row <- summary[as.character(summary$display_group) == group_label, , drop = FALSE]
      row <- row[order(row$distance_kb), , drop = FALSE]
      if (nrow(row)) lines(row$distance_kb, row$mean_reads,
                           col = if (identical(effective_color_by, "none")) "grey35" else qc_library_group_color(row$role[[1]]),
                           lty = 1, lwd = max(2.2, qc_plot_lwd * 1.5))
    }
  }
  qc_plot_metric_key_row(style$key$label, col = style$key$color, pch = style$key$pch,
                         lty = style$key$lty, lwd = qc_plot_lwd,
                         ncol = min(3L, nrow(style$key)), title = "Key")
}

qc_library_metaplot_plot_data <- function(project) {
  roots <- c(file.path(project$project_root, "library_diagnostics", "runs", "all_samples"),
             file.path(project$project_root, "library_diagnostics"))
  patterns <- c(tss = "\\.tss_metaplot\\.tsv$", tts = "\\.tts_metaplot\\.tsv$", trna = "\\.trna_metaplot\\.tsv$")
  rows <- list()
  for (feature in names(patterns)) {
    root_files <- lapply(roots, function(root) {
      if (dir.exists(root)) list.files(root, pattern = patterns[[feature]], recursive = TRUE, full.names = TRUE) else character()
    })
    available <- which(vapply(root_files, length, integer(1)) > 0L)
    files <- if (length(available)) root_files[[available[[1]]]] else character()
    rows <- c(rows, lapply(files, function(file) {
      data <- tryCatch(read.delim(file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
      if (is.null(data) || !all(c("rel_bp", "reads_sum", "sites_n") %in% names(data))) return(NULL)
      sample <- sub(paste0("\\.", feature, "_metaplot\\.tsv$"), "", basename(file))
      data.frame(sample = sample, feature = feature,
                 rel_bp = qc_plot_numeric(data$rel_bp),
                 reads_sum = qc_plot_numeric(data$reads_sum),
                 sites_n = qc_plot_numeric(data$sites_n),
                 stringsAsFactors = FALSE)
    }))
  }
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data) || !nrow(data)) return(data.frame())
  data$mean_reads_per_site <- ifelse(data$sites_n > 0, data$reads_sum / data$sites_n, NA_real_)
  data
}

plot_qc_library_metaplot <- function(project, panel = "tss", group = "all", color_by = "group") {
  data <- qc_library_metaplot_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("Feature-metaplot tables are not available."))
  data <- data[data$feature == panel, , drop = FALSE]
  data <- qc_library_filter_data(project, data, group)
  if (!nrow(data)) return(qc_plot_empty("No feature-metaplot data are available for the selected sample group."))
  feature_label <- c(tss = "TSS metaplot", tts = "TTS metaplot", trna = "tRNA metaplot")[[panel]] %||% "Feature metaplot"
  effective_color_by <- if (color_by %in% c("sample", "none")) color_by else "group"
  style <- qc_library_sample_style(data, effective_color_by)
  aesthetic <- qc_library_plot_aesthetics("metaplot")
  old <- par(no.readonly = TRUE)
  on.exit({ qc_plot_reset_layout(); par(old) }, add = TRUE)
  qc_plot_begin_key_layout(0.16)
  qc_plot_par(mar = c(5.5, 5.5, 3.5, 1))
  plot(NA, xlim = range(data$rel_bp, na.rm = TRUE), ylim = range(data$mean_reads_per_site, na.rm = TRUE),
       xlab = "Position relative to feature (bp)", ylab = "Mean reads per site",
       main = feature_label)
  if (qc_plot_grid_enabled()) {
    abline(h = pretty(range(data$mean_reads_per_site, na.rm = TRUE)), col = "#e3e9ee", lty = 3)
    abline(v = 0, col = "grey45", lty = 2, lwd = qc_plot_lwd)
  }
  if (identical(effective_color_by, "sample")) {
    for (i in seq_along(style$samples)) {
      row <- data[as.character(data$sample) == style$samples[[i]], , drop = FALSE]
      row <- row[order(row$rel_bp), , drop = FALSE]
      if (nrow(row)) lines(row$rel_bp, row$mean_reads_per_site,
                           col = style$colors[[i]], lty = 1, lwd = qc_plot_lwd)
    }
  } else if (isTRUE(aesthetic$summary_line)) {
    summary <- aggregate(mean_reads_per_site ~ display_group + role + rel_bp, data, mean, na.rm = TRUE)
    for (group_label in unique(as.character(summary$display_group))) {
      row <- summary[as.character(summary$display_group) == group_label, , drop = FALSE]
      row <- row[order(row$rel_bp), , drop = FALSE]
      if (nrow(row)) lines(row$rel_bp, row$mean_reads_per_site,
                           col = if (identical(effective_color_by, "none")) "grey35" else qc_library_group_color(row$role[[1]]),
                           lty = 1, lwd = max(2.2, qc_plot_lwd * 1.5))
    }
  }
  qc_plot_metric_key_row(style$key$label, col = style$key$color, pch = style$key$pch,
                         lty = style$key$lty, lwd = qc_plot_lwd,
                         ncol = min(3L, nrow(style$key)), title = "Key")
}

qc_library_plot_inventory <- function(project, plot_types = c("Centromere bias", "Feature metaplots")) {
  inventory <- build_diagnostic_file_inventory(project)
  if (!nrow(inventory)) return(data.frame())
  plots <- inventory[inventory$is_plot & inventory$plot_type %in% plot_types & inventory$preview_available, , drop = FALSE]
  if (any(plots$run_id != "current_legacy")) plots <- plots[plots$run_id != "current_legacy", , drop = FALSE]
  plots[order(plots$plot_type, plots$sample, decreasing = FALSE), , drop = FALSE]
}

qc_library_plot_gallery <- function(project, plot_types = c("Centromere bias", "Feature metaplots")) {
  plots <- qc_library_plot_inventory(project, plot_types)
  if (!nrow(plots)) return(tags$p(class = "text-muted", "No Library Diagnostics plot files are available for this view."))
  tags$div(class = "ytab-plot-grid ytab-diagnostic-gallery-grid", lapply(seq_len(min(nrow(plots), 12L)), function(i) {
    tags$article(
      class = "ytab-static-image-card ytab-release-card",
      title = plots$filename[[i]],
      tags$div(class = "ytab-plot-card-header",
               tags$div(tags$h4(plots$plot_type[[i]])),
               tags$a(class = "btn btn-secondary btn-sm", href = plots$served_url[[i]], target = "_blank", "Open original")),
      tags$p(class = "ytab-plot-sample", plots$sample[[i]]),
      tags$div(class = "ytab-diagnostic-preview",
        tags$img(src = plots$served_url[[i]], alt = plots$filename[[i]],
                 loading = "lazy", style = "width:100%;max-height:320px;object-fit:contain")
      ),
      tags$details(tags$summary("Plot provenance"), tags$p("Source: generated PNG"), tags$code(plots$filename[[i]]))
    )
  }))
}

qc_library_single_plot_card <- function(plots, index) {
  if (!nrow(plots) || is.na(index) || index < 1L || index > nrow(plots))
    return(tags$p(class = "text-muted", "Selected diagnostic plot is unavailable."))
  tags$article(
    class = "ytab-static-image-card ytab-release-card",
    title = plots$filename[[index]],
    tags$div(class = "ytab-plot-card-header",
             tags$div(tags$h4(plots$plot_type[[index]])),
             tags$a(class = "btn btn-secondary btn-sm", href = plots$served_url[[index]], target = "_blank", "Open original")),
    tags$p(class = "ytab-plot-sample", plots$sample[[index]]),
    tags$div(class = "ytab-diagnostic-preview",
      tags$img(src = plots$served_url[[index]], alt = plots$filename[[index]],
               loading = "lazy", style = "width:100%;max-height:560px;object-fit:contain")
    ),
    tags$details(tags$summary("Plot provenance"), tags$p("Source: generated PNG"), tags$code(plots$filename[[index]]))
  )
}
