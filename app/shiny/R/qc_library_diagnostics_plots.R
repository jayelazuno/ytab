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

plot_qc_library_midlc <- function(project) {
  data <- qc_library_midlc_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("MidLC curves are not available."))
  old <- qc_plot_par(mar = c(5, 5, 3, 1))
  on.exit(par(old), add = TRUE)
  plot(NA, xlim = range(data$reads_sampled, na.rm = TRUE),
       ylim = range(data$unique_sites, na.rm = TRUE), log = "x",
       xlab = "Reads sampled", ylab = "Unique insertion sites",
       main = "MidLC saturation curves")
  samples <- levels(data$sample)
  colors <- grDevices::rainbow(length(samples))
  for (i in seq_along(samples)) {
    row <- data[as.character(data$sample) == samples[[i]], , drop = FALSE]
    if (nrow(row)) lines(row$reads_sampled, row$unique_sites, type = "b", pch = 16, col = "black", lwd = qc_plot_lwd)
  }
  legend("bottomright", legend = samples, col = "black", pch = 16, lty = 1, lwd = qc_plot_lwd, cex = 0.75, bty = "n", text.font = 2)
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

plot_qc_library_jackpot_depth <- function(project) {
  data <- qc_library_summary_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("Library diagnostic summary is not available."))
  labels <- as.character(data$sample)
  jackpot <- qc_plot_numeric(qc_plot_column(data, c("jackpot_frac_reads", "jackpot_top_frac")))
  depth <- qc_plot_numeric(qc_plot_column(data, "depth_ratio_R_over_midlc"))
  old <- qc_plot_par(mar = c(10, 5, 3, 4))
  on.exit(par(old), add = TRUE)
  barplot(jackpot * 100, names.arg = labels, las = 2, col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
          ylab = "Jackpot reads (%)", main = "Jackpots and library depth")
  par(new = TRUE)
  plot(seq_along(depth), depth, type = "b", pch = 16, col = "black", lwd = qc_plot_lwd,
       axes = FALSE, xlab = "", ylab = "")
  axis(4)
  mtext("Reads / MidLC", side = 4, line = 2.4)
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
    data.frame(sample = sample, condition = qc_plot_condition(sample), pair = data$pair,
               enrichment = qc_plot_numeric(data$enrichment), stringsAsFactors = FALSE)
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data)) data.frame() else data
}

plot_qc_library_sequence_bias <- function(project) {
  data <- qc_library_sequence_bias_plot_data(project)
  if (!nrow(data)) return(qc_plot_empty("Sequence-bias tables are not available."))
  summarized <- aggregate(enrichment ~ condition + pair, data, mean, na.rm = TRUE)
  pairs <- sort(unique(summarized$pair))
  conditions <- c("parent", "H2O2-treated-facs", "other")
  conditions <- intersect(conditions, unique(summarized$condition))
  mat <- sapply(conditions, function(condition) {
    row <- summarized[summarized$condition == condition, ]
    stats::setNames(row$enrichment, row$pair)[pairs]
  })
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = length(pairs), dimnames = list(pairs, conditions))
  old <- qc_plot_par(mar = c(8, 5, 3, 1))
  on.exit(par(old), add = TRUE)
  barplot(t(mat), beside = TRUE, names.arg = pairs, las = 2,
          col = c(qc_plot_fill, qc_plot_fill_light, "grey60"),
          border = qc_plot_border, lwd = qc_plot_lwd,
          xlab = "Dinucleotide (+2/+7)", ylab = "Mean enrichment",
          main = "Insertion-site sequence bias")
  abline(h = 1, lty = 2, lwd = qc_plot_lwd, col = "black")
  legend("topright", legend = conditions,
         fill = c(qc_plot_fill, qc_plot_fill_light, "grey60")[seq_along(conditions)],
         bty = "n", text.font = 2, cex = 1)
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
  tags$div(class = "ytab-plot-grid", lapply(seq_len(min(nrow(plots), 12L)), function(i) {
    tags$article(
      class = "ytab-plot-card",
      title = plots$filename[[i]],
      tags$h4(plots$plot_type[[i]]),
      tags$p(class = "ytab-plot-sample", plots$sample[[i]]),
      tags$div(class = "ytab-diagnostic-preview",
        tags$img(src = plots$served_url[[i]], alt = plots$filename[[i]],
                 loading = "lazy", style = "width:100%;height:170px;object-fit:contain")
      ),
      tags$details(tags$summary("Show filename"), tags$code(plots$filename[[i]]))
    )
  }))
}

qc_library_single_plot_card <- function(plots, index) {
  if (!nrow(plots) || is.na(index) || index < 1L || index > nrow(plots))
    return(tags$p(class = "text-muted", "Selected diagnostic plot is unavailable."))
  tags$article(
    class = "ytab-plot-card",
    title = plots$filename[[index]],
    tags$h4(plots$plot_type[[index]]),
    tags$p(class = "ytab-plot-sample", plots$sample[[index]]),
    tags$div(class = "ytab-diagnostic-preview",
      tags$img(src = plots$served_url[[index]], alt = plots$filename[[index]],
               loading = "lazy", style = "width:100%;max-height:560px;object-fit:contain")
    ),
    tags$details(tags$summary("Show filename"), tags$code(plots$filename[[index]]))
  )
}
