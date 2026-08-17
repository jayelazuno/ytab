qc_mapping_stats_plot_data <- function(project) {
  data <- read_qc_csvs(file.path(project$project_root, "mapfastq"), "mapping_stats.*\\.csv$")
  if (!nrow(data)) return(data)
  if (!"sample" %in% names(data)) data$sample <- basename(dirname(data$.source))
  data$sample <- factor(as.character(data$sample), levels = qc_plot_sample_order(as.character(data$sample)))
  data[order(data$sample), , drop = FALSE]
}

plot_qc_mapping_stats <- function(project) {
  data <- qc_mapping_stats_plot_data(project)
  if (!nrow(data) || !"total_records" %in% names(data))
    return(qc_plot_empty("Mapping statistics are not available."))

  total <- qc_plot_numeric(data$total_records)
  mapped <- qc_plot_numeric(qc_plot_column(data, "percent_mapped"))
  hq <- qc_plot_numeric(qc_plot_column(data, c("percent_mapq_ge20", "percent_mapq_ge_threshold")))
  raw_labels <- as.character(data$sample)
  labels <- qc_plot_display_labels(raw_labels)
  horizontal <- qc_plot_bar_horizontal(raw_labels)
  show_values <- qc_plot_show_value_labels()
  label_angle <- qc_plot_label_angle(0L)
  label_mode <- qc_plot_label_mode()
  left_margin <- if (horizontal && !identical(label_mode, "hide")) {
    max(18, min(34, max(nchar(labels), na.rm = TRUE) * 0.8))
  } else 5
  bottom_margin <- if (horizontal) 5.5 else if (label_angle %in% c(30L, 45L)) 9 else if (label_angle == 90L) 11 else 7

  old <- qc_plot_par(mar = c(bottom_margin, left_margin, 4.2, 16))
  on.exit(par(old), add = TRUE)

  mapped[!is.finite(mapped)] <- NA_real_
  hq[!is.finite(hq)] <- NA_real_
  x_limit <- c(0, 104)
  if (horizontal) {
    ypos <- barplot(mapped, names.arg = labels, horiz = TRUE, las = 1,
                    xlim = x_limit, col = qc_plot_fill, border = qc_plot_border,
                    lwd = qc_plot_lwd, xlab = "% reads mapped", ylab = "",
                    main = "Mapping summary")
    if (qc_plot_grid_enabled()) abline(v = seq(0, 100, by = 20), col = "#e3e9ee", lty = 3)
    points(hq, ypos, pch = 17, col = "black", cex = 1.15, lwd = qc_plot_lwd)
    if (show_values) {
      text(pmin(mapped + 2, 102), ypos,
           labels = sprintf("%.1f%% mapped · %.1f%% HQ · %s reads",
                            mapped, hq, format(round(total), big.mark = ",")),
           pos = 4, cex = 0.78 * qc_plot_text_scale(), font = 2, xpd = NA)
    }
  } else {
    ypos <- barplot(mapped, names.arg = rep("", length(labels)), ylim = x_limit,
                    col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
                    ylab = "% reads mapped", xlab = "", main = "Mapping summary")
    if (qc_plot_grid_enabled()) abline(h = seq(0, 100, by = 20), col = "#e3e9ee", lty = 3)
    points(ypos, hq, pch = 17, col = "black", cex = 1.15, lwd = qc_plot_lwd)
    if (!identical(label_mode, "hide")) {
      axis(1, at = ypos, labels = FALSE)
      if (label_angle %in% c(30L, 45L)) {
        text(ypos, par("usr")[3] - 4, labels = labels, srt = label_angle,
             adj = 1, xpd = NA, cex = 0.8 * qc_plot_text_scale(), font = 2)
      } else {
        axis(1, at = ypos, labels = labels, las = if (label_angle == 90L) 2 else 1)
      }
    }
    if (show_values) {
      text(ypos, pmin(mapped + 2, 102), labels = sprintf("%.1f%%", mapped),
           cex = 0.78 * qc_plot_text_scale(), font = 2, xpd = NA)
    }
  }
  legend("right", inset = c(-0.48, 0), xpd = NA, horiz = FALSE,
         legend = c("% mapped (bar)", "% HQ aligned (triangle)"),
         fill = c(qc_plot_fill, NA), border = c(qc_plot_border, NA),
         pch = c(NA, 17), col = c(qc_plot_border, "black"),
         bty = "n", text.font = 2, cex = 0.9 * qc_plot_text_scale())
}
