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
  hq_source <- intersect(c("percent_mapq_ge_threshold", "percent_mapq_ge20"), names(data))
  hq <- qc_plot_numeric(qc_plot_column(data, hq_source))
  raw_labels <- as.character(data$sample)
  labels <- qc_plot_display_labels(raw_labels)
  horizontal <- qc_plot_bar_horizontal(raw_labels)
  show_values <- qc_plot_show_value_labels()
  label_angle <- qc_plot_label_angle(0L)
  label_mode <- qc_plot_label_mode()
  hq_available <- length(hq_source) > 0L && any(is.finite(hq))
  metric_labels <- c("% reads mapped", if (hq_available) "% HQ aligned" else character())
  metric_pch <- c(15, if (hq_available) 17 else integer())
  metric_col <- c(qc_plot_border, if (hq_available) "black" else character())
  left_margin <- if (horizontal && !identical(label_mode, "hide")) {
    qc_plot_label_margin_lines(labels, orientation = "horizontal")
  } else 6
  bottom_margin <- if (horizontal) {
    8.5
  } else if (identical(label_mode, "hide")) {
    5.5
  } else {
    qc_plot_label_margin_lines(labels, angle = label_angle, orientation = "vertical") + 2
  }

  old <- par(no.readonly = TRUE)
  on.exit({ par(old); layout(1) }, add = TRUE)
  layout(matrix(c(1, 2), ncol = 1), heights = c(1, 0.12))
  qc_plot_par(mar = c(bottom_margin, left_margin, 4.5, 3.2))

  mapped[!is.finite(mapped)] <- NA_real_
  hq[!is.finite(hq)] <- NA_real_
  x_limit <- c(0, 104)
  if (horizontal) {
    ypos <- barplot(mapped, names.arg = labels, horiz = TRUE, las = 1,
                    xlim = x_limit, col = qc_plot_fill, border = qc_plot_border,
                    lwd = qc_plot_lwd, xlab = "% reads mapped", ylab = "",
                    main = "Mapping summary")
    if (qc_plot_grid_enabled()) abline(v = seq(0, 100, by = 20), col = "#e3e9ee", lty = 3)
    if (hq_available) points(hq, ypos, pch = 17, col = "black", cex = 1.15, lwd = qc_plot_lwd)
    if (show_values) {
      value_text <- if (hq_available) {
        sprintf("%.1f%% mapped · %.1f%% HQ · %s reads", mapped, hq, format(round(total), big.mark = ","))
      } else {
        sprintf("%.1f%% mapped · %s reads", mapped, format(round(total), big.mark = ","))
      }
      text(pmin(mapped + 2, 102), ypos,
           labels = value_text,
           pos = 4, cex = qc_plot_value_cex(0.82), font = 2, xpd = NA)
    }
  } else {
    ypos <- barplot(mapped, names.arg = rep("", length(labels)), ylim = x_limit,
                    col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
                    ylab = "% reads mapped", xlab = "", main = "Mapping summary")
    if (qc_plot_grid_enabled()) abline(h = seq(0, 100, by = 20), col = "#e3e9ee", lty = 3)
    if (hq_available) points(ypos, hq, pch = 17, col = "black", cex = 1.15, lwd = qc_plot_lwd)
    if (!identical(label_mode, "hide")) {
      axis(1, at = ypos, labels = FALSE)
      axis_labels <- if (label_angle == 0L) qc_plot_wrapped_labels(labels) else labels
      if (label_angle %in% c(30L, 45L)) {
        text(ypos, par("usr")[3] - 4, labels = axis_labels, srt = label_angle,
             adj = 1, xpd = NA, cex = qc_plot_label_cex(0.82), font = 2)
      } else if (label_angle == 90L) {
        text(ypos, par("usr")[3] - 4, labels = labels, srt = 90,
             adj = 1, xpd = NA, cex = qc_plot_label_cex(0.82), font = 2)
      } else {
        axis(1, at = ypos, labels = axis_labels, las = 1)
      }
    }
    if (show_values) {
      text(ypos, pmin(mapped + 2, 102), labels = sprintf("%.1f%%", mapped),
           cex = qc_plot_value_cex(0.82), font = 2, xpd = NA)
    }
  }
  par(mar = c(0, left_margin, 0, 3.2), font.axis = 2, font.lab = 2, lwd = 1.35)
  plot.new()
  legend("center", horiz = TRUE,
         legend = metric_labels, pch = metric_pch, col = metric_col,
         pt.bg = c(qc_plot_fill, if (hq_available) "black" else character()),
         pt.cex = c(1.35, if (hq_available) 1.25 else numeric()),
         bty = "n", text.font = 2, cex = qc_plot_key_cex(0.95),
         x.intersp = 0.8, y.intersp = 0.9)
}
