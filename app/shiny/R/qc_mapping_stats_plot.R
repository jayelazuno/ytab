qc_mapping_stats_plot_data <- function(project) {
  data <- read_qc_csvs(file.path(project$project_root, "mapfastq"), "mapping_stats.*\\.csv$")
  if (!nrow(data)) return(data)
  if (!"sample" %in% names(data)) data$sample <- basename(dirname(data$.source))
  data$sample <- factor(as.character(data$sample), levels = qc_plot_sample_order(as.character(data$sample)))
  data[order(data$sample), , drop = FALSE]
}

plot_qc_mapping_stats_legacy_stacked <- function(project) {
  data <- qc_mapping_stats_plot_data(project)
  if (!nrow(data) || !"total_records" %in% names(data))
    return(qc_plot_empty("Mapping statistics are not available."))

  total <- qc_plot_numeric(data$total_records)
  primary_mapped <- qc_plot_numeric(qc_plot_column(data, "primary_mapped"))
  percent_mapped <- qc_plot_numeric(qc_plot_column(data, "percent_mapped"))
  missing_primary <- !is.finite(primary_mapped) & is.finite(total) & is.finite(percent_mapped)
  primary_mapped[missing_primary] <- total[missing_primary] * percent_mapped[missing_primary] / 100
  missing_percent <- !is.finite(percent_mapped) & is.finite(total) & total > 0 & is.finite(primary_mapped)
  percent_mapped[missing_percent] <- 100 * primary_mapped[missing_percent] / total[missing_percent]
  mapq_source <- intersect(
    c("avg_mapq_mapped_primary", "mean_mapq_mapped_primary", "avg_mapq_primary", "mean_mapq"),
    names(data)
  )
  mapped_mapq <- qc_plot_numeric(qc_plot_column(data, mapq_source))
  raw_labels <- as.character(data$sample)
  labels <- qc_plot_display_labels(raw_labels)
  horizontal <- qc_plot_bar_horizontal(raw_labels)
  show_values <- qc_plot_show_value_labels()
  label_angle <- qc_plot_label_angle(0L)
  label_mode <- qc_plot_label_mode()
  mapq_available <- length(mapq_source) > 0L && any(is.finite(mapped_mapq))
  left_margin <- if (horizontal && !identical(label_mode, "hide")) {
    min(qc_plot_label_margin_lines(labels, orientation = "horizontal"), 22)
  } else 6
  bottom_margin <- if (horizontal) {
    3.2
  } else if (identical(label_mode, "hide")) {
    4.2
  } else {
    min(qc_plot_label_margin_lines(labels, angle = label_angle, orientation = "vertical") + 2, 7)
  }

  old <- par(no.readonly = TRUE)
  on.exit({ par(old); layout(1) }, add = TRUE)
  layout(matrix(c(1, 2, 3, 4), ncol = 1), heights = c(1.35, 1, 1, 0.18))

  total_millions <- total / 1e6
  mapped_millions <- primary_mapped / 1e6
  total_millions[!is.finite(total_millions)] <- NA_real_
  mapped_millions[!is.finite(mapped_millions)] <- NA_real_
  percent_mapped[!is.finite(percent_mapped)] <- NA_real_
  mapped_mapq[!is.finite(mapped_mapq)] <- NA_real_

  read_cols <- c("#d8e4ef", "#2f6f9f")
  read_border <- c("#8ca8bf", "#1f4f73")
  read_matrix <- rbind(`Total reads` = total_millions, `Mapped reads` = mapped_millions)
  if (horizontal) {
    qc_plot_par(mar = c(3.2, left_margin, 3.1, 2.4))
    x_limit_counts <- c(0, max(read_matrix, na.rm = TRUE) * 1.12)
    if (!all(is.finite(x_limit_counts))) x_limit_counts <- c(0, 1)
    count_pos <- barplot(
      read_matrix, beside = TRUE, names.arg = labels, horiz = TRUE, las = 1,
      xlim = x_limit_counts, col = read_cols, border = read_border,
      lwd = qc_plot_lwd, xlab = "Reads (millions)", ylab = "",
      main = "Mapping QC: total and mapped reads"
    )
    if (qc_plot_grid_enabled()) abline(v = pretty(x_limit_counts), col = "#e3e9ee", lty = 3)
    if (show_values) {
      text(mapped_millions + max(x_limit_counts) * 0.015, count_pos[2, ],
           labels = sprintf("%.1fM", mapped_millions),
           pos = 4, cex = qc_plot_value_cex(0.76), font = 2, xpd = NA)
    }

    qc_plot_par(mar = c(3.0, left_margin, 2.5, 2.4))
    pct_pos <- barplot(
      percent_mapped, names.arg = labels, horiz = TRUE, las = 1,
      xlim = c(0, 104), col = qc_plot_fill, border = qc_plot_border,
      lwd = qc_plot_lwd, xlab = "Mapped reads (%)", ylab = "",
      main = "Percent mapped reads"
    )
    if (qc_plot_grid_enabled()) abline(v = seq(0, 100, by = 20), col = "#e3e9ee", lty = 3)
    if (show_values) {
      text(pmin(percent_mapped + 1.8, 102), pct_pos,
           labels = sprintf("%.1f%%", percent_mapped),
           pos = 4, cex = qc_plot_value_cex(0.78), font = 2, xpd = NA)
    }

    qc_plot_par(mar = c(3.0, left_margin, 2.5, 2.4))
    mapq_limit <- c(0, max(45, max(mapped_mapq, na.rm = TRUE) * 1.08))
    if (!mapq_available || !all(is.finite(mapq_limit))) mapq_limit <- c(0, 45)
    mapq_pos <- barplot(
      mapped_mapq, names.arg = labels, horiz = TRUE, las = 1,
      xlim = mapq_limit, col = "#c9d7c2", border = "#5b7654",
      lwd = qc_plot_lwd, xlab = "Average MAPQ", ylab = "",
      main = "Mapping quality score"
    )
    if (qc_plot_grid_enabled()) abline(v = pretty(mapq_limit), col = "#e3e9ee", lty = 3)
    if (show_values && mapq_available) {
      text(mapped_mapq + max(mapq_limit) * 0.015, mapq_pos,
           labels = sprintf("%.1f", mapped_mapq),
           pos = 4, cex = qc_plot_value_cex(0.78), font = 2, xpd = NA)
    }
  } else {
    qc_plot_par(mar = c(2.5, 5.3, 3.1, 2.0))
    y_limit_counts <- c(0, max(read_matrix, na.rm = TRUE) * 1.14)
    if (!all(is.finite(y_limit_counts))) y_limit_counts <- c(0, 1)
    count_pos <- barplot(
      read_matrix, beside = TRUE, names.arg = rep("", length(labels)),
      ylim = y_limit_counts, col = read_cols, border = read_border,
      lwd = qc_plot_lwd, ylab = "Reads (millions)", xlab = "",
      main = "Mapping QC: total and mapped reads"
    )
    if (qc_plot_grid_enabled()) abline(h = pretty(y_limit_counts), col = "#e3e9ee", lty = 3)
    if (show_values) {
      text(count_pos[2, ], mapped_millions + max(y_limit_counts) * 0.025,
           labels = sprintf("%.1fM", mapped_millions),
           cex = qc_plot_value_cex(0.72), font = 2, xpd = NA)
    }

    qc_plot_par(mar = c(2.5, 5.3, 2.5, 2.0))
    pct_pos <- barplot(
      percent_mapped, names.arg = rep("", length(labels)), ylim = c(0, 104),
      col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
      ylab = "Mapped reads (%)", xlab = "", main = "Percent mapped reads"
    )
    if (qc_plot_grid_enabled()) abline(h = seq(0, 100, by = 20), col = "#e3e9ee", lty = 3)
    if (show_values) {
      text(pct_pos, pmin(percent_mapped + 2, 102), labels = sprintf("%.1f%%", percent_mapped),
           cex = qc_plot_value_cex(0.78), font = 2, xpd = NA)
    }

    qc_plot_par(mar = c(bottom_margin, 5.3, 2.5, 2.0))
    mapq_limit <- c(0, max(45, max(mapped_mapq, na.rm = TRUE) * 1.08))
    if (!mapq_available || !all(is.finite(mapq_limit))) mapq_limit <- c(0, 45)
    mapq_pos <- barplot(
      mapped_mapq, names.arg = rep("", length(labels)), ylim = mapq_limit,
      col = "#c9d7c2", border = "#5b7654", lwd = qc_plot_lwd,
      ylab = "Average MAPQ", xlab = "", main = "Mapping quality score"
    )
    if (qc_plot_grid_enabled()) abline(h = pretty(mapq_limit), col = "#e3e9ee", lty = 3)
    if (!identical(label_mode, "hide")) {
      axis(1, at = mapq_pos, labels = FALSE)
      axis_labels <- if (label_angle == 0L) qc_plot_wrapped_labels(labels) else labels
      label_y <- par("usr")[3] - max(mapq_limit) * 0.055
      if (label_angle %in% c(30L, 45L)) {
        text(mapq_pos, label_y, labels = axis_labels, srt = label_angle,
             adj = 1, xpd = NA, cex = qc_plot_label_cex(0.82), font = 2)
      } else if (label_angle == 90L) {
        text(mapq_pos, label_y, labels = labels, srt = 90,
             adj = 1, xpd = NA, cex = qc_plot_label_cex(0.82), font = 2)
      } else {
        axis(1, at = mapq_pos, labels = axis_labels, las = 1)
      }
    }
    if (show_values && mapq_available) {
      text(mapq_pos, mapped_mapq + max(mapq_limit) * 0.025,
           labels = sprintf("%.1f", mapped_mapq),
           cex = qc_plot_value_cex(0.78), font = 2, xpd = NA)
    }
  }
  par(mar = c(0, left_margin, 0, 3.2), font.axis = 2, font.lab = 2, lwd = 1.35)
  plot.new()
  legend("center", horiz = TRUE,
         legend = c("Total reads", "Mapped reads", "% reads mapped", "Average MAPQ mapped primary"),
         fill = c(read_cols, qc_plot_fill, "#c9d7c2"),
         border = c(read_border, qc_plot_border, "#5b7654"),
         bty = "n", text.font = 2, cex = qc_plot_key_cex(0.95),
         x.intersp = 0.8, y.intersp = 0.9)
}

plot_qc_mapping_stats <- function(project, metric = "read_counts") {
  data <- qc_mapping_stats_plot_data(project)
  if (!nrow(data) || !"total_records" %in% names(data))
    return(qc_plot_empty("Mapping statistics are not available."))

  metric <- as.character(metric %||% "read_counts")
  if (!metric %in% c("read_counts", "percent_mapped", "mapq")) metric <- "read_counts"

  total <- qc_plot_numeric(data$total_records)
  primary_mapped <- qc_plot_numeric(qc_plot_column(data, "primary_mapped"))
  percent_mapped <- qc_plot_numeric(qc_plot_column(data, "percent_mapped"))
  missing_primary <- !is.finite(primary_mapped) & is.finite(total) & is.finite(percent_mapped)
  primary_mapped[missing_primary] <- total[missing_primary] * percent_mapped[missing_primary] / 100
  missing_percent <- !is.finite(percent_mapped) & is.finite(total) & total > 0 & is.finite(primary_mapped)
  percent_mapped[missing_percent] <- 100 * primary_mapped[missing_percent] / total[missing_percent]

  mapq_source <- intersect(
    c("avg_mapq_mapped_primary", "mean_mapq_mapped_primary", "avg_mapq_primary", "mean_mapq"),
    names(data)
  )
  mapped_mapq <- qc_plot_numeric(qc_plot_column(data, mapq_source))

  raw_labels <- as.character(data$sample)
  labels <- qc_plot_display_labels(raw_labels)
  horizontal <- qc_plot_bar_horizontal(raw_labels)
  show_values <- qc_plot_show_value_labels()
  label_angle <- qc_plot_label_angle(0L)
  label_mode <- qc_plot_label_mode()
  left_margin <- if (horizontal && !identical(label_mode, "hide")) {
    qc_plot_label_margin_lines(labels, orientation = "horizontal")
  } else 6
  bottom_margin <- if (horizontal) {
    6
  } else if (identical(label_mode, "hide")) {
    5.5
  } else {
    qc_plot_label_margin_lines(labels, angle = label_angle, orientation = "vertical") + 2
  }

  draw_vertical_labels <- function(at, plot_labels, value_range) {
    if (identical(label_mode, "hide")) return(invisible())
    axis(1, at = at, labels = FALSE)
    axis_labels <- if (label_angle == 0L) qc_plot_wrapped_labels(plot_labels) else plot_labels
    label_y <- par("usr")[3] - max(value_range, na.rm = TRUE) * 0.035
    if (label_angle %in% c(30L, 45L)) {
      text(at, label_y, labels = axis_labels, srt = label_angle,
           adj = 1, xpd = NA, cex = qc_plot_label_cex(0.82), font = 2)
    } else if (label_angle == 90L) {
      text(at, label_y, labels = plot_labels, srt = 90,
           adj = 1, xpd = NA, cex = qc_plot_label_cex(0.82), font = 2)
    } else {
      axis(1, at = at, labels = axis_labels, las = 1)
    }
    invisible()
  }

  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)

  if (identical(metric, "read_counts")) {
    total_millions <- total / 1e6
    mapped_millions <- primary_mapped / 1e6
    total_millions[!is.finite(total_millions)] <- NA_real_
    mapped_millions[!is.finite(mapped_millions)] <- NA_real_
    read_matrix <- rbind(`Total reads` = total_millions, `Mapped reads` = mapped_millions)
    read_cols <- c("#d8e4ef", "#2f6f9f")
    read_border <- c("#8ca8bf", "#1f4f73")
    if (horizontal) {
      qc_plot_par(mar = c(6.2, left_margin, 4.2, 2.8))
      x_limit <- c(0, max(read_matrix, na.rm = TRUE) * 1.12)
      if (!all(is.finite(x_limit))) x_limit <- c(0, 1)
      pos <- barplot(read_matrix, beside = TRUE, names.arg = labels, horiz = TRUE,
                     las = 1, xlim = x_limit, col = read_cols, border = read_border,
                     lwd = qc_plot_lwd, xlab = "Reads (millions)", ylab = "",
                     main = "Mapping summary — total and mapped reads")
      if (qc_plot_grid_enabled()) abline(v = pretty(x_limit), col = "#e3e9ee", lty = 3)
      if (show_values) {
        text(mapped_millions + max(x_limit) * 0.015, pos[2, ],
             labels = sprintf("%.1fM", mapped_millions), pos = 4,
             cex = qc_plot_value_cex(0.82), font = 2, xpd = NA)
      }
    } else {
      qc_plot_par(mar = c(bottom_margin, 6, 4.2, 2.8))
      y_limit <- c(0, max(read_matrix, na.rm = TRUE) * 1.14)
      if (!all(is.finite(y_limit))) y_limit <- c(0, 1)
      pos <- barplot(read_matrix, beside = TRUE, names.arg = rep("", length(labels)),
                     ylim = y_limit, col = read_cols, border = read_border,
                     lwd = qc_plot_lwd, ylab = "Reads (millions)", xlab = "",
                     main = "Mapping summary — total and mapped reads")
      if (qc_plot_grid_enabled()) abline(h = pretty(y_limit), col = "#e3e9ee", lty = 3)
      draw_vertical_labels(colMeans(pos), labels, y_limit)
      if (show_values) {
        text(pos[2, ], mapped_millions + max(y_limit) * 0.025,
             labels = sprintf("%.1fM", mapped_millions),
             cex = qc_plot_value_cex(0.78), font = 2, xpd = NA)
      }
    }
    return(invisible())
  }

  if (identical(metric, "percent_mapped")) {
    percent_mapped[!is.finite(percent_mapped)] <- NA_real_
    if (horizontal) {
      qc_plot_par(mar = c(6.2, left_margin, 4.2, 2.8))
      pos <- barplot(percent_mapped, names.arg = labels, horiz = TRUE, las = 1,
                     xlim = c(0, 104), col = qc_plot_fill, border = qc_plot_border,
                     lwd = qc_plot_lwd, xlab = "Mapped reads (%)", ylab = "",
                     main = "Mapping summary — percent mapped reads")
      if (qc_plot_grid_enabled()) abline(v = seq(0, 100, by = 20), col = "#e3e9ee", lty = 3)
      if (show_values) {
        text(pmin(percent_mapped + 1.8, 102), pos,
             labels = sprintf("%.1f%%", percent_mapped), pos = 4,
             cex = qc_plot_value_cex(0.82), font = 2, xpd = NA)
      }
    } else {
      qc_plot_par(mar = c(bottom_margin, 6, 4.2, 2.8))
      pos <- barplot(percent_mapped, names.arg = rep("", length(labels)),
                     ylim = c(0, 104), col = qc_plot_fill, border = qc_plot_border,
                     lwd = qc_plot_lwd, ylab = "Mapped reads (%)", xlab = "",
                     main = "Mapping summary — percent mapped reads")
      if (qc_plot_grid_enabled()) abline(h = seq(0, 100, by = 20), col = "#e3e9ee", lty = 3)
      draw_vertical_labels(pos, labels, c(0, 104))
      if (show_values) {
        text(pos, pmin(percent_mapped + 2, 102), labels = sprintf("%.1f%%", percent_mapped),
             cex = qc_plot_value_cex(0.82), font = 2, xpd = NA)
      }
    }
    return(invisible())
  }

  mapped_mapq[!is.finite(mapped_mapq)] <- NA_real_
  mapq_limit <- c(0, max(45, max(mapped_mapq, na.rm = TRUE) * 1.08))
  if (!length(mapq_source) || !all(is.finite(mapq_limit))) mapq_limit <- c(0, 45)
  if (horizontal) {
    qc_plot_par(mar = c(6.2, left_margin, 4.2, 2.8))
    pos <- barplot(mapped_mapq, names.arg = labels, horiz = TRUE, las = 1,
                   xlim = mapq_limit, col = "#c9d7c2", border = "#5b7654",
                   lwd = qc_plot_lwd, xlab = "Average MAPQ", ylab = "",
                   main = "Mapping summary — mapping quality score")
    if (qc_plot_grid_enabled()) abline(v = pretty(mapq_limit), col = "#e3e9ee", lty = 3)
    if (show_values && length(mapq_source)) {
      text(mapped_mapq + max(mapq_limit) * 0.015, pos,
           labels = sprintf("%.1f", mapped_mapq), pos = 4,
           cex = qc_plot_value_cex(0.82), font = 2, xpd = NA)
    }
  } else {
    qc_plot_par(mar = c(bottom_margin, 6, 4.2, 2.8))
    pos <- barplot(mapped_mapq, names.arg = rep("", length(labels)),
                   ylim = mapq_limit, col = "#c9d7c2", border = "#5b7654",
                   lwd = qc_plot_lwd, ylab = "Average MAPQ", xlab = "",
                   main = "Mapping summary — mapping quality score")
    if (qc_plot_grid_enabled()) abline(h = pretty(mapq_limit), col = "#e3e9ee", lty = 3)
    draw_vertical_labels(pos, labels, mapq_limit)
    if (show_values && length(mapq_source)) {
      text(pos, mapped_mapq + max(mapq_limit) * 0.025,
           labels = sprintf("%.1f", mapped_mapq),
           cex = qc_plot_value_cex(0.82), font = 2, xpd = NA)
    }
  }
  invisible()
}
