fitness_condition_control_plot_data <- function(result) {
  if (is.null(result) || !nzchar(result$output_dir %||% "")) return(data.frame())
  reads_path <- file.path(result$output_dir, "tables", "feature_table.long.raw_reads.csv")
  if (!file.exists(reads_path))
    reads_path <- file.path(result$output_dir, "tables", "feature_reads_cpm.long.csv")
  if (!file.exists(reads_path) || !file.exists(result$table)) return(data.frame())

  reads <- tryCatch(read.csv(reads_path, stringsAsFactors = FALSE, check.names = FALSE),
                    error = function(e) data.frame())
  calls <- tryCatch(read.csv(result$table, stringsAsFactors = FALSE, check.names = FALSE),
                    error = function(e) data.frame())
  if (!nrow(reads) || !nrow(calls) || !"condition" %in% names(reads)) return(data.frame())

  id_col <- intersect(c("feature_id", "standard_name", "gene_id"), names(reads))
  read_col <- intersect(c("reads", "raw_reads", "count", "cpm"), names(reads))
  if (!length(id_col) || !length(read_col)) return(data.frame())
  id_col <- id_col[[1]]
  read_col <- read_col[[1]]

  reads$condition_group <- ifelse(grepl("parent", reads$condition, ignore.case = TRUE), "parent",
                                  ifelse(grepl("treated|H2O2|facs", reads$condition, ignore.case = TRUE),
                                         "treated", "other"))
  reads$feature_id_tmp <- as.character(reads[[id_col]])
  reads$value <- qc_plot_numeric(reads[[read_col]])
  summarized <- aggregate(value ~ feature_id_tmp + condition_group,
                          reads,
                          mean, na.rm = TRUE)
  names(summarized)[names(summarized) == "feature_id_tmp"] <- "feature_id"
  wide <- reshape(summarized, idvar = "feature_id", timevar = "condition_group", direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  if (!all(c("parent", "treated") %in% names(wide))) return(data.frame())

  call_id <- intersect(c("feature_id", "standard_name", "gene_id"), names(calls))
  if (!length(call_id)) return(data.frame())
  call_id <- call_id[[1]]
  keep <- intersect(c(call_id, "final_call", "mean_log2FC", "standard_name", "common_name"), names(calls))
  calls <- calls[, keep, drop = FALSE]
  names(calls)[names(calls) == call_id] <- "feature_id"
  calls <- calls[!duplicated(calls$feature_id), , drop = FALSE]

  data <- merge(wide, calls, by = "feature_id", all.x = TRUE)
  data$final_call <- as.character(data$final_call %||% "")
  data$mean_log2FC <- qc_plot_numeric(data$mean_log2FC)
  data$label <- if ("common_name" %in% names(data)) as.character(data$common_name) else ""
  data$label[is.na(data$label) | !nzchar(trimws(data$label))] <-
    if ("standard_name" %in% names(data)) as.character(data$standard_name[is.na(data$label) | !nzchar(trimws(data$label))]) else data$feature_id[is.na(data$label) | !nzchar(trimws(data$label))]
  data$label[is.na(data$label) | !nzchar(trimws(data$label))] <- data$feature_id[is.na(data$label) | !nzchar(trimws(data$label))]
  data
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

plot_fitness_condition_control_scatter <- function(result) {
  data <- fitness_condition_control_plot_data(result)
  if (!nrow(data)) return(qc_plot_empty("Condition-versus-control read data are not available."))

  parent <- qc_plot_numeric(data$parent) + 1
  treated <- qc_plot_numeric(data$treated) + 1
  keep <- is.finite(parent) & is.finite(treated) & parent > 0 & treated > 0
  if (!any(keep)) return(qc_plot_empty("Condition-versus-control read data are not available."))

  calls <- as.character(data$final_call)
  point_col <- ifelse(grepl("^consistently_", calls), "black",
                      ifelse(grepl("depleted", calls), "#2f6fb5",
                             ifelse(grepl("enriched", calls), "#c83f3f", "grey70")))
  point_col[is.na(point_col)] <- "grey70"
  highlight <- fitness_condition_control_highlights(data, 10L)
  highlight <- intersect(highlight, which(keep))

  old <- qc_plot_par(mar = c(5.5, 5.5, 3.5, 1))
  on.exit(par(old), add = TRUE)
  lim <- range(c(parent[keep], treated[keep]), na.rm = TRUE)
  plot(parent[keep], treated[keep], log = "xy", pch = 16, cex = 0.65,
       col = adjustcolor(point_col[keep], alpha.f = 0.65),
       xlim = lim, ylim = lim,
       xlab = "Parent reads per gene",
       ylab = "H2O2-treated reads per gene",
       main = "H2O2-treated versus parent insertion reads")
  abline(0, 1, lwd = qc_plot_lwd, lty = 2, col = "black")
  if (length(highlight)) {
    points(parent[highlight], treated[highlight], pch = 21, bg = "black",
           col = "black", cex = 1.05, lwd = qc_plot_lwd)
    text(parent[highlight], treated[highlight], labels = data$label[highlight],
         pos = rep(c(2, 4, 3), length.out = length(highlight)),
         cex = 0.85, font = 2, col = "black", xpd = NA)
  }
  legend(
    "topleft",
    legend = c("not significant / unchanged",
               "depleted in H2O2",
               "enriched in H2O2",
               "consistently changed or highlighted"),
    col = c("grey70", "#2f6fb5", "#c83f3f", "black"),
    pch = 16, bty = "n", text.font = 2, cex = 0.9
  )
  mtext("Below diagonal: mutants depleted in H2O2; above diagonal: mutants enriched in H2O2.",
        side = 3, line = 0.15, cex = 0.9, font = 2)
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
