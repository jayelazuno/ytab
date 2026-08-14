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
  labels <- as.character(data$sample)
  old <- qc_plot_par(mar = c(10, 5, 3, 1))
  on.exit(par(old), add = TRUE)

  if (identical(metric, "features")) {
    barplot(data$percent_features_hit, names.arg = labels, las = 2,
            col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
            ylab = "% features hit", main = "Genomic features hit")
  } else if (identical(metric, "reads_per_hit")) {
    barplot(log10(pmax(data$mean_reads_per_hit, 1)), names.arg = labels, las = 2,
            col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
            ylab = "Mean reads per hit (log10)", main = "Reads per insertion site")
  } else if (identical(metric, "feature_intergenic")) {
    values <- rbind(data$percent_hits_in_features, data$percent_intergenic_hits)
    barplot(values, names.arg = labels, las = 2, col = c(qc_plot_fill, qc_plot_fill_light),
            border = qc_plot_border, lwd = qc_plot_lwd, ylab = "Percent of hits",
            main = "Feature vs intergenic insertions")
    legend("topright", legend = c("Features", "Intergenic"),
           fill = c(qc_plot_fill, qc_plot_fill_light), bty = "n", text.font = 2, cex = 1)
  } else {
    barplot(data$total_hits, names.arg = labels, las = 2,
            col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
            ylab = "Total unique insertion sites", main = "Library complexity")
  }
}

qc_summary_combined_features_hit_data <- function(project) {
  root <- file.path(project$project_root, "summary")
  files <- if (dir.exists(root)) {
    list.files(root, pattern = "\\.feature_table\\.RDF_1\\.csv$", recursive = TRUE, full.names = TRUE)
  } else character()
  rows <- lapply(files, function(file) {
    sample <- basename(dirname(file))
    data <- tryCatch(read.csv(file, skip = 1, stringsAsFactors = FALSE, check.names = FALSE),
                     error = function(e) NULL)
    if (is.null(data) || !"standard_name" %in% names(data) || !"hits" %in% names(data)) return(NULL)
    data.frame(
      condition = qc_plot_condition(sample),
      feature_id = as.character(data$standard_name),
      hits = qc_plot_numeric(data$hits),
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data) || !nrow(data)) return(data.frame())
  data <- data[data$condition %in% c("parent", "H2O2-treated-facs") &
                 !is.na(data$feature_id) & nzchar(data$feature_id), , drop = FALSE]
  if (!nrow(data)) return(data.frame())
  summed <- aggregate(hits ~ condition + feature_id, data, sum, na.rm = TRUE)
  out <- do.call(rbind, lapply(split(summed, summed$condition), function(group) {
    total <- length(unique(group$feature_id))
    hit <- sum(qc_plot_numeric(group$hits) > 0, na.rm = TRUE)
    data.frame(
      condition = unique(group$condition)[[1]],
      features_hit = hit,
      total_features = total,
      percent_features_hit = if (total > 0) 100 * hit / total else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(out)) return(data.frame())
  out$condition <- factor(out$condition, levels = c("parent", "H2O2-treated-facs"))
  out[order(out$condition), , drop = FALSE]
}

plot_qc_summary_combined_features_hit <- function(project) {
  data <- qc_summary_combined_features_hit_data(project)
  if (!nrow(data)) return(qc_plot_empty("Combined feature-hit data are not available."))
  labels <- c(parent = "Parents combined", `H2O2-treated-facs` = "H2O2-treated combined")
  names <- unname(labels[as.character(data$condition)])
  values <- data$percent_features_hit
  old <- qc_plot_par(mar = c(7, 5, 3, 1))
  on.exit(par(old), add = TRUE)
  ypos <- barplot(values, names.arg = names, las = 2, ylim = c(0, max(100, values, na.rm = TRUE)),
                  col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
                  ylab = "% features hit", main = "Combined features hit by condition")
  text(ypos, values, labels = sprintf("%s / %s\n%.1f%%",
                                      format(data$features_hit, big.mark = ","),
                                      format(data$total_features, big.mark = ","),
                                      data$percent_features_hit),
       pos = 3, font = 2, cex = 0.95)
}

qc_summary_binned_plot_data <- function(project) {
  root <- file.path(project$project_root, "summary")
  files <- if (dir.exists(root)) list.files(root, pattern = "^binned_hits\\.RDF_1\\.csv$", recursive = TRUE, full.names = TRUE) else character()
  rows <- lapply(files, function(file) {
    data <- tryCatch(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(data) || ncol(data) < 2L) return(NULL)
    bin_id <- as.character(data[[1L]])
    chrom <- sub("-[0-9]+$", "", bin_id)
    aggregate(qc_plot_numeric(data[[2L]]), list(chrom = chrom), function(x) mean(x, na.rm = TRUE))
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data) || !nrow(data)) return(data.frame())
  names(data)[names(data) == "x"] <- "mean_bin_signal"
  data
}

plot_qc_summary_genome_bins <- function(project) {
  data <- qc_summary_binned_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("Genome-wide binned hit files are not available."))
  data$chrom <- factor(data$chrom, levels = unique(data$chrom))
  old <- qc_plot_par(mar = c(8, 5, 3, 1))
  on.exit(par(old), add = TRUE)
  boxplot(log10(data$mean_bin_signal + 1) ~ data$chrom, las = 2,
          col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
          ylab = "log10(mean binned signal + 1)",
          main = "Genome-wide binned insertion signal")
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
  labels <- paste(data$condition, data$pair, sep = ": ")
  old <- qc_plot_par(mar = c(5, 12, 3, 1))
  on.exit(par(old), add = TRUE)
  y <- seq_len(nrow(data))
  plot(data$spearman, y, xlim = c(0, 1), yaxt = "n", pch = 16, col = "black",
       xlab = "Spearman correlation", ylab = "", main = "Library concordance")
  axis(2, at = y, labels = labels, las = 1)
  abline(v = seq(0, 1, by = 0.25), col = "#e3e9ee", lty = 3)
}
