qc_summary_stats_plot_data <- function(project) {
  root <- file.path(project$project_root, "summary")
  files <- if (dir.exists(root)) list.files(root, pattern = "^stats\\.csv$", recursive = TRUE, full.names = TRUE) else character()
  rows <- lapply(files, function(file) {
    row <- tryCatch(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(row) || !nrow(row)) return(NULL)
    sample <- basename(dirname(file))
    data.frame(
      sample = sample,
      condition = qc_plot_condition(sample),
      total_hits = qc_plot_numeric(qc_plot_column(row, c("Total Hits", "total_hits"), nrow(row))),
      percent_hits_in_features = qc_plot_numeric(qc_plot_column(row, c("% of hits in features", "percent_hits_in_features"), nrow(row))),
      percent_intergenic_hits = qc_plot_numeric(qc_plot_column(row, c("% of intergenic hits", "percent_intergenic_hits"), nrow(row))),
      percent_features_hit = qc_plot_numeric(qc_plot_column(row, c("% of features hit", "percent_features_hit"), nrow(row))),
      mean_reads_per_hit = qc_plot_numeric(qc_plot_column(row, c("Mean Reads Per Hit", "mean_reads_per_hit"), nrow(row))),
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data) || !nrow(data)) return(data.frame())
  data$sample <- factor(data$sample, levels = qc_plot_sample_order(as.character(data$sample)))
  data[order(data$sample), , drop = FALSE]
}

plot_qc_summary_metric <- function(project, metric = "complexity") {
  data <- qc_summary_stats_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("SummaryTable statistics are not available."))
  raw_labels <- as.character(data$sample)
  display_labels <- qc_plot_display_labels(raw_labels)
  horizontal <- qc_plot_bar_horizontal(display_labels)
  has_key <- identical(metric, "feature_intergenic")
  old <- par(no.readonly = TRUE)
  on.exit({ qc_plot_reset_layout(); par(old) }, add = TRUE)
  if (has_key) qc_plot_begin_key_layout(0.14)
  qc_plot_par(mar = qc_plot_bar_margins(display_labels, horizontal, top = 3.8, right = if (horizontal) 5.2 else 2.8))

  draw_bar_labels <- function(at) {
    if (horizontal) qc_plot_draw_horizontal_labels(at, display_labels)
    else qc_plot_draw_vertical_labels(at, display_labels)
  }
  draw_grid <- function() qc_plot_add_grid(horizontal)
  draw_values <- function(at, values, suffix = "") {
    if (!qc_plot_show_value_labels()) return(invisible())
    labels <- paste0(format(round(values, 1), trim = TRUE, big.mark = ","), suffix)
    if (horizontal) {
      text(values, at, labels = labels, pos = 4, xpd = NA, font = 2, cex = qc_plot_label_cex(0.85))
    } else {
      text(at, values, labels = labels, pos = 3, xpd = NA, font = 2, cex = qc_plot_label_cex(0.85))
    }
  }

  if (identical(metric, "features")) {
    values <- data$percent_features_hit
    limit <- max(100, values, na.rm = TRUE)
    at <- barplot(values, names.arg = rep("", length(values)), axes = FALSE, horiz = horizontal,
            xlim = if (horizontal) c(0, limit * 1.08) else NULL,
            ylim = if (horizontal) NULL else c(0, limit * 1.08),
            col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
            xlab = if (horizontal) "% features hit" else "",
            ylab = if (horizontal) "" else "% features hit", main = "Genomic features hit")
    draw_grid()
    if (horizontal) axis(1) else axis(2)
    draw_bar_labels(at)
    draw_values(at, values, "%")
  } else if (identical(metric, "reads_per_hit")) {
    values <- log10(pmax(data$mean_reads_per_hit, 1))
    limit <- max(values, na.rm = TRUE)
    at <- barplot(values, names.arg = rep("", length(values)), axes = FALSE, horiz = horizontal,
            xlim = if (horizontal) c(0, limit * 1.12) else NULL,
            ylim = if (horizontal) NULL else c(0, limit * 1.12),
            col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
            xlab = if (horizontal) "Mean reads per hit (log10)" else "",
            ylab = if (horizontal) "" else "Mean reads per hit (log10)", main = "Reads per insertion site")
    draw_grid()
    if (horizontal) axis(1) else axis(2)
    draw_bar_labels(at)
    draw_values(at, values)
  } else if (identical(metric, "feature_intergenic")) {
    values <- rbind(data$percent_hits_in_features, data$percent_intergenic_hits)
    totals <- colSums(values, na.rm = TRUE)
    limit <- max(100, totals, na.rm = TRUE)
    at <- barplot(values, names.arg = rep("", length(raw_labels)), axes = FALSE, horiz = horizontal, col = c(qc_plot_fill, qc_plot_fill_light),
            border = qc_plot_border, lwd = qc_plot_lwd,
            xlim = if (horizontal) c(0, limit * 1.08) else NULL,
            ylim = if (horizontal) NULL else c(0, limit * 1.08),
            xlab = if (horizontal) "Percent of hits" else "",
            ylab = if (horizontal) "" else "Percent of hits",
            main = "Feature vs intergenic insertions")
    draw_grid()
    if (horizontal) axis(1) else axis(2)
    draw_bar_labels(at)
    draw_values(at, data$percent_hits_in_features, "% features")
    qc_plot_metric_key_row(c("Features", "Intergenic"),
                           fill = c(qc_plot_fill, qc_plot_fill_light),
                           border = qc_plot_border, ncol = 2)
  } else {
    values <- data$total_hits
    limit <- max(values, na.rm = TRUE)
    at <- barplot(values, names.arg = rep("", length(values)), axes = FALSE, horiz = horizontal,
            xlim = if (horizontal) c(0, limit * 1.12) else NULL,
            ylim = if (horizontal) NULL else c(0, limit * 1.12),
            col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
            xlab = if (horizontal) "Total unique insertion sites" else "",
            ylab = if (horizontal) "" else "Total unique insertion sites", main = "Library complexity")
    draw_grid()
    if (horizontal) axis(1) else axis(2)
    draw_bar_labels(at)
    draw_values(at, values)
  }
}

qc_summary_combined_features_hit_data <- function(project) {
  root <- file.path(project$project_root, "summary")
  files <- if (dir.exists(root)) {
    list.files(root, pattern = "\\.feature_table\\.RDF_1\\.csv$", recursive = TRUE, full.names = TRUE)
  } else character()
  roles <- qc_summary_combined_feature_sample_roles(project)
  rows <- lapply(files, function(file) {
    sample <- basename(dirname(file))
    data <- tryCatch(read.csv(file, skip = 1, stringsAsFactors = FALSE, check.names = FALSE),
                     error = function(e) NULL)
    if (is.null(data) || !"standard_name" %in% names(data) || !"hits" %in% names(data)) return(NULL)
    role <- roles$role[match(sample, roles$sample)]
    data.frame(
      sample = sample,
      group = role,
      feature_id = as.character(data$standard_name),
      hits = qc_plot_numeric(data$hits),
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data) || !nrow(data)) return(data.frame())
  data <- data[data$group %in% c("control", "treated") &
                 !is.na(data$feature_id) & nzchar(data$feature_id), , drop = FALSE]
  if (!nrow(data)) return(data.frame())
  summed <- aggregate(hits ~ group + feature_id, data, sum, na.rm = TRUE)
  labels <- qc_summary_combined_feature_group_labels(project)
  out <- do.call(rbind, lapply(split(summed, summed$group), function(group_data) {
    total <- length(unique(group_data$feature_id))
    hit <- sum(qc_plot_numeric(group_data$hits) > 0, na.rm = TRUE)
    group_id <- unique(group_data$group)[[1]]
    data.frame(
      group = group_id,
      display_group = labels[[group_id]] %||% group_id,
      features_hit = hit,
      total_features = total,
      percent_features_hit = if (total > 0) 100 * hit / total else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(out)) return(data.frame())
  out$group <- factor(out$group, levels = c("control", "treated"))
  out[order(out$group), , drop = FALSE]
}

qc_summary_project_samples <- function(project) {
  if (is.data.frame(project$samples) && nrow(project$samples)) return(project$samples)
  path <- as.character(project$sample_sheet %||% "")
  if (nzchar(path)) {
    candidates <- c(path, file.path(project$repo_root %||% "", path), file.path(project$project_root %||% "", "config", basename(path)))
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates)) {
      return(tryCatch(read.csv(candidates[[1]], stringsAsFactors = FALSE, check.names = FALSE),
                      error = function(e) data.frame()))
    }
  }
  data.frame()
}

qc_summary_role_from_values <- function(values) {
  values <- tolower(paste(stats::na.omit(as.character(values)), collapse = " "))
  if (!nzchar(values)) return(NA_character_)
  if (grepl("\\btreated\\b|zn-treated|h2o2-treated|1_5mm|h2o2", values)) return("treated")
  if (grepl("\\bcontrol\\b|\\bmock\\b|\\bparent\\b|classifier_control", values)) return("control")
  NA_character_
}

qc_summary_combined_feature_sample_roles <- function(project) {
  root <- file.path(project$project_root, "summary")
  samples <- if (dir.exists(root)) basename(dirname(list.files(root, pattern = "\\.feature_table\\.RDF_1\\.csv$", recursive = TRUE, full.names = TRUE))) else character()
  samples <- qc_plot_sample_order(unique(samples))
  metadata <- qc_summary_project_samples(project)
  role_cols <- c("library_role", "control_or_treated", "fitness_role", "classifier_role",
                 "condition", "condition_label", "treatment", "treatment_label",
                 "guessed_condition")
  roles <- vapply(samples, function(sample) {
    row <- if (is.data.frame(metadata) && "sample" %in% names(metadata)) metadata[as.character(metadata$sample) == sample, , drop = FALSE] else data.frame()
    if (nrow(row)) {
      for (col in intersect(role_cols, names(row))) {
        role <- qc_summary_role_from_values(row[[col]])
        if (!is.na(role)) return(role)
      }
    }
    sample_lower <- tolower(sample)
    if (grepl("treated|h2o2|zn-treated|1_5mm", sample_lower)) return("treated")
    if (grepl("mock|control|parent", sample_lower)) return("control")
    NA_character_
  }, character(1))
  data.frame(sample = samples, role = roles, stringsAsFactors = FALSE)
}

qc_summary_combined_feature_group_labels <- function(project) {
  metadata <- qc_summary_project_samples(project)
  controls <- if (is.data.frame(metadata) && nrow(metadata)) {
    roles <- qc_summary_combined_feature_sample_roles(project)
    metadata[match(roles$sample[roles$role == "control"], as.character(metadata$sample)), , drop = FALSE]
  } else data.frame()
  treated <- if (is.data.frame(metadata) && nrow(metadata)) {
    roles <- qc_summary_combined_feature_sample_roles(project)
    metadata[match(roles$sample[roles$role == "treated"], as.character(metadata$sample)), , drop = FALSE]
  } else data.frame()

  control_label <- "Controls combined"
  if (nrow(controls)) {
    control_text <- tolower(paste(unique(unlist(controls[, intersect(c("condition", "condition_label", "treatment", "treatment_label"), names(controls)), drop = FALSE])), collapse = " "))
    control_label <- if (grepl("mock", control_text)) "Mock controls combined" else if (grepl("parent", control_text)) "Parents combined" else "Controls combined"
  }

  treated_label <- "Treated combined"
  if (nrow(treated)) {
    label_cols <- intersect(c("condition_label", "treatment_label", "condition", "treatment"), names(treated))
    labels <- unique(as.character(unlist(treated[, label_cols, drop = FALSE])))
    labels <- labels[nzchar(labels) & !is.na(labels)]
    if (length(labels)) {
      if (any(grepl("H2O2", labels, ignore.case = TRUE))) {
        treated_label <- "H2O2-treated combined"
      } else {
        preferred <- labels[grepl("treated", labels, ignore.case = TRUE) &
                              !tolower(labels) %in% c("treated", "zn-treated")]
        if (length(preferred)) treated_label <- paste0(preferred[[1]], " combined")
        else {
          specific <- labels[!tolower(labels) %in% c("treated", "control", "mock", "parent")]
          treated_label <- paste0(specific[[1]] %||% "Treated", "-treated combined")
        }
      }
    }
  }
  list(control = control_label, treated = treated_label)
}

qc_summary_combined_feature_group_choices <- function(project) {
  data <- qc_summary_combined_features_hit_data(project)
  labels <- qc_summary_combined_feature_group_labels(project)
  base_label <- function(label) {
    label <- sub(" combined$", "", label)
    if (grepl("^[A-Z][a-z]", label)) {
      paste0(tolower(substr(label, 1, 1)), substr(label, 2, nchar(label)))
    } else label
  }
  choices <- list()
  if (any(data$group == "control")) {
    choices[[paste("Combined", base_label(labels$control %||% "controls"))]] <- "control"
  }
  if (any(data$group == "treated")) {
    choices[[paste("Combined", base_label(labels$treated %||% "treated"))]] <- "treated"
  }
  if (all(c("control", "treated") %in% as.character(data$group))) {
    both_label <- paste("Combined", paste(base_label(labels$control),
                                          base_label(labels$treated), sep = " + "))
    choices[[both_label]] <- "both"
  }
  choices
}

qc_summary_combined_feature_default_group <- function(project) {
  choices <- unname(unlist(qc_summary_combined_feature_group_choices(project), use.names = FALSE))
  if ("both" %in% choices) "both" else choices[[1]] %||% "control"
}

plot_qc_summary_combined_features_hit <- function(project, group = "both") {
  data <- qc_summary_combined_features_hit_data(project)
  if (!nrow(data)) return(qc_plot_empty("Combined feature-hit data are not available."))
  available <- as.character(data$group)
  if (identical(group, "both") && all(c("control", "treated") %in% available)) {
    data <- data[data$group %in% c("control", "treated"), , drop = FALSE]
    title <- "Combined features hit: controls vs treated"
  } else {
    group <- if (group %in% available) group else qc_summary_combined_feature_default_group(project)
    if (identical(group, "both")) group <- available[[1]]
    data <- data[as.character(data$group) == group, , drop = FALSE]
    title <- if (identical(group, "treated")) "Combined features hit: treated" else "Combined features hit: controls"
  }
  if (!nrow(data)) return(qc_plot_empty("Selected combined feature-hit group is not available."))
  names <- data$display_group
  values <- data$percent_features_hit
  display_names <- qc_plot_display_labels(names)
  horizontal <- qc_plot_bar_horizontal(display_names)
  old <- qc_plot_par(mar = qc_plot_bar_margins(display_names, horizontal, top = 3.8, right = 5))
  on.exit(par(old), add = TRUE)
  limit <- max(100, values, na.rm = TRUE)
  n <- length(values)
  bar_half <- if (n == 1L) 0.22 else 0.30
  pos <- seq_len(n)
  if (horizontal) {
    plot(NA, xlim = c(0, limit * 1.18), ylim = c(0.45, n + 0.55),
         axes = FALSE, xlab = "% features hit", ylab = "", main = title)
    qc_plot_add_grid(TRUE)
    rect(0, pos - bar_half, values, pos + bar_half,
         col = qc_plot_fill, border = qc_plot_border, lwd = max(1.6, qc_plot_lwd))
    axis(1)
    qc_plot_draw_horizontal_labels(pos, display_names, cex = qc_plot_label_cex(1.05))
  } else {
    plot(NA, xlim = if (n == 1L) c(0.35, 1.65) else c(0.45, n + 0.55),
         ylim = c(0, limit * 1.18), axes = FALSE,
         xlab = "", ylab = "% features hit", main = title)
    qc_plot_add_grid(FALSE)
    rect(pos - bar_half, 0, pos + bar_half, values,
         col = qc_plot_fill, border = qc_plot_border, lwd = max(1.6, qc_plot_lwd))
    axis(2)
    qc_plot_draw_vertical_labels(pos, display_names, cex = qc_plot_label_cex(1.05))
  }
  if (qc_plot_show_value_labels()) {
    text(if (horizontal) values else pos, if (horizontal) pos else values,
         labels = sprintf("%s / %s\n%.1f%%",
                          format(data$features_hit, big.mark = ","),
                          format(data$total_features, big.mark = ","),
                          data$percent_features_hit),
         pos = if (horizontal) 4 else 3, font = 2, xpd = NA,
         cex = qc_plot_label_cex(0.95))
  }
}

qc_summary_binned_plot_data <- function(project) {
  root <- file.path(project$project_root, "summary")
  files <- if (dir.exists(root)) list.files(root, pattern = "^binned_hits\\.RDF_1\\.csv$", recursive = TRUE, full.names = TRUE) else character()
  rows <- lapply(files, function(file) {
    data <- tryCatch(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(data) || ncol(data) < 2L) return(NULL)
    bin_id <- as.character(data[[1L]])
    chrom <- sub("-[0-9]+$", "", bin_id)
    out <- aggregate(qc_plot_numeric(data[[2L]]), list(chrom_raw = chrom), function(x) mean(x, na.rm = TRUE))
    out$chrom <- qc_plot_chromosome_display(out$chrom_raw, project)
    out
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data) || !nrow(data)) return(data.frame())
  names(data)[names(data) == "x"] <- "mean_bin_signal"
  data
}

plot_qc_summary_genome_bins <- function(project) {
  data <- qc_summary_binned_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("Genome-wide binned hit files are not available."))
  data$chrom <- factor(data$chrom, levels = qc_plot_chromosome_order(data$chrom_raw, project))
  labels <- levels(data$chrom)
  display_labels <- qc_plot_display_labels(labels)
  old <- qc_plot_par(mar = c(qc_plot_label_margin_lines(display_labels,
                                                        angle = qc_plot_label_angle(45L),
                                                        orientation = "vertical"),
                            5, 3, 1))
  on.exit(par(old), add = TRUE)
  boxplot(log10(mean_bin_signal + 1) ~ chrom, data = data, xaxt = "n",
          col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
          xlab = "",
          ylab = "log10(mean binned signal + 1)",
          main = "Genome-wide binned insertion signal")
  qc_plot_add_grid(FALSE)
  qc_plot_draw_vertical_labels(seq_along(labels), labels, angle = qc_plot_label_angle(45L))
  mtext("Chromosome", side = 1, line = max(3.6, par("mar")[[1L]] - 1.4),
        font = 2, cex = qc_plot_text_sizes()$lab)
}

qc_summary_pairwise_plot_data <- function(project) {
  root <- file.path(project$project_root, "summary")
  files <- if (dir.exists(root)) list.files(root, pattern = "\\.feature_table\\.RDF_1\\.csv$", recursive = TRUE, full.names = TRUE) else character()
  rows <- lapply(files, function(file) {
    sample <- basename(dirname(file))
    data <- tryCatch(read.csv(file, skip = 1, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(data) || !"standard_name" %in% names(data) || !"reads" %in% names(data)) return(NULL)
    data.frame(sample = sample, standard_name = data$standard_name,
               reads = qc_plot_numeric(data$reads), stringsAsFactors = FALSE)
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data) || !nrow(data)) return(data.frame())
  groups <- split(unique(as.character(data$sample)), qc_plot_condition(unique(as.character(data$sample))))
  out <- list()
  for (group_name in names(groups)) {
    samples <- qc_plot_sample_order(groups[[group_name]])
    if (length(samples) < 2L) next
    out <- c(out, lapply(combn(samples, 2, simplify = FALSE), function(pair) {
      merged <- merge(data[data$sample == pair[[1]], c("standard_name", "reads")],
                      data[data$sample == pair[[2]], c("standard_name", "reads")],
                      by = "standard_name", suffixes = c("_1", "_2"))
      data.frame(condition = group_name,
                 pair = paste(qc_plot_pool(pair[[1]]), qc_plot_pool(pair[[2]]), sep = " vs "),
                 spearman = suppressWarnings(cor(merged$reads_1, merged$reads_2, method = "spearman", use = "pairwise.complete.obs")),
                 stringsAsFactors = FALSE)
    }))
  }
  data <- do.call(rbind, out)
  if (is.null(data)) data.frame() else data
}

plot_qc_summary_pairwise_correlations <- function(project) {
  data <- qc_summary_pairwise_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("Feature-table correlations are not available."))
  labels <- qc_plot_display_labels(paste(data$condition, data$pair, sep = ": "))
  old <- qc_plot_par(mar = c(5, qc_plot_label_margin_lines(labels, orientation = "horizontal"), 3, 1))
  on.exit(par(old), add = TRUE)
  y <- seq_len(nrow(data))
  plot(data$spearman, y, xlim = c(0, 1), yaxt = "n", pch = 16, col = "black",
       xlab = "Spearman correlation", ylab = "", main = "Library concordance")
  qc_plot_draw_horizontal_labels(y, labels)
  if (qc_plot_grid_enabled()) abline(v = seq(0, 1, by = 0.25), col = "#e3e9ee", lty = 3)
  if (qc_plot_show_value_labels()) text(data$spearman, y, labels = sprintf("%.2f", data$spearman),
                                        pos = 4, xpd = NA, font = 2,
                                        cex = qc_plot_label_cex(0.85))
}
