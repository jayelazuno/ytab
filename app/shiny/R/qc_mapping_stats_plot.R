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
  labels <- as.character(data$sample)
  scale_factor <- max(total, na.rm = TRUE) / 100
  if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1

  old <- qc_plot_par(mar = c(5, 12, 3, 4))
  on.exit(par(old), add = TRUE)
  ypos <- barplot(total / scale_factor, names.arg = labels, horiz = TRUE, las = 1,
                  xlim = c(0, 100), col = qc_plot_fill, border = qc_plot_border,
                  lwd = qc_plot_lwd, xlab = "Percent", main = "Mapping summary")
  points(mapped, ypos, pch = 16, col = "black", cex = 1.25, lwd = qc_plot_lwd)
  points(hq, ypos, pch = 17, col = "black", cex = 1.25, lwd = qc_plot_lwd)
  legend("bottomright", legend = c("% mapped", "% HQ"), pch = c(16, 17),
         col = c("black", "black"), bty = "n", text.font = 2, cex = 1)
  ticks <- pretty(total / scale_factor)
  axis(3, at = ticks, labels = format(round(ticks * scale_factor), big.mark = ","))
  mtext("Total reads", side = 3, line = 2.2)
}
