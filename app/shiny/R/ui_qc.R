relative_project_path <- function(path, project_root) {
  path <- normalizePath(path, winslash="/", mustWork=FALSE); root <- normalizePath(project_root,winslash="/",mustWork=FALSE)
  if (startsWith(path,paste0(root,"/"))) substring(path,nchar(root)+2L) else basename(path)
}

read_qc_csvs <- function(root, pattern) {
  files <- if(dir.exists(root)) list.files(root,pattern=pattern,recursive=TRUE,full.names=TRUE) else character()
  rows <- lapply(files,function(file)tryCatch(transform(read.csv(file,stringsAsFactors=FALSE,check.names=FALSE),.source=file),error=function(e)NULL))
  rows <- Filter(Negate(is.null),rows); if(length(rows))do.call(rbind,rows)else data.frame()
}

format_qc_metrics <- function(data, count_columns=character(), percentage_columns=character(), mean_columns=character()) {
  for(name in intersect(count_columns,names(data))) data[[name]] <- ifelse(is.na(data[[name]]),"",format(round(as.numeric(data[[name]])),big.mark=",",scientific=FALSE,trim=TRUE))
  for(name in intersect(percentage_columns,names(data))) { value<-suppressWarnings(as.numeric(sub("%$","",as.character(data[[name]])))); data[[name]]<-ifelse(is.na(value),"",sprintf("%.1f%%",value)) }
  for(name in intersect(mean_columns,names(data))) data[[name]] <- ifelse(is.na(data[[name]]),"",sprintf("%.2f",as.numeric(data[[name]])))
  data
}

mapping_qc_data <- function(project) {
  raw<-read_qc_csvs(file.path(project$project_root,"mapfastq"),"mapping_stats.*\\.csv$"); if(!nrow(raw))return(raw)
  sample<-raw$sample %||% basename(dirname(raw$.source)); bam<-file.path(project$project_root,"mapfastq",sample,paste0(sample,".sorted.bam"))
  result<-data.frame(Sample=ytab_sample_labels(sample,project),check.names=FALSE,stringsAsFactors=FALSE)
  available<-list(`Reads processed`="total_records",`Aligned reads`="primary_mapped",`Alignment percentage`="percent_mapped")
  for(label in names(available))if(available[[label]]%in%names(raw))result[[label]]<-raw[[available[[label]]]]
  result$BAM<-ifelse(file.exists(bam),"Available","Missing");result$`Mapping statistics`<-"Available"
  attr(result,"details")<-data.frame(Sample=sample,BAM=vapply(bam,relative_project_path,"",project_root=project$project_root),Statistics=vapply(raw$.source,relative_project_path,"",project_root=project$project_root),stringsAsFactors=FALSE)
  format_qc_metrics(result,c("Reads processed","Aligned reads"),"Alignment percentage")
}

summary_qc_data <- function(project) {
  raw<-read_qc_csvs(file.path(project$project_root,"summary"),"summary_stats.*\\.csv$");if(!nrow(raw))return(raw)
  wanted<-c("sample","total_reads","total_hits","percent_hits_in_features","percent_intergenic_hits","percent_features_hit","mean_hits_per_feature"); visible<-raw[intersect(wanted,names(raw))]
  names(visible)<-c(sample="Sample",total_reads="Total reads",total_hits="Total hits",percent_hits_in_features="Hits in features",percent_intergenic_hits="Intergenic hits",percent_features_hit="Features hit",mean_hits_per_feature="Mean hits / feature")[names(visible)]
  detail_names<-setdiff(names(raw),c(wanted,".source"));details<-raw[,c("sample",detail_names),drop=FALSE];details$source<-vapply(raw$.source,relative_project_path,"",project_root=project$project_root);attr(visible,"details")<-details
  if ("Sample" %in% names(visible)) visible$Sample <- ytab_sample_labels(visible$Sample, project)
  format_qc_metrics(visible,c("Total reads","Total hits"),c("Hits in features","Intergenic hits","Features hit"),"Mean hits / feature")
}

summary_qc_cards <- function(project) {
  raw<-read_qc_csvs(file.path(project$project_root,"summary"),"summary_stats.*\\.csv$"); num<-function(name){if(!name%in%names(raw))numeric()else suppressWarnings(as.numeric(sub("%$","",raw[[name]])))}; safe<-function(x,fun,digits=0){x<-x[is.finite(x)];if(!length(x))"Unavailable"else format(round(fun(x),digits),big.mark=",",trim=TRUE)}
  tags$div(
    class = "ytab-qc-summary-list",
    tags$p(class = "ytab-qc-summary-caption", "Per-sample SummaryTable overview"),
    tags$dl(
      tags$div(tags$dt("Included samples"), tags$dd(format(nrow(raw), big.mark = ","))),
      tags$div(tags$dt("Median unique insertion sites / sample"), tags$dd(safe(num("total_hits"), median))),
      tags$div(tags$dt("Median genes hit / sample"), tags$dd(paste0(safe(num("percent_features_hit"), median, 1), "%"))),
      tags$div(tags$dt("Lowest genes-hit sample"), tags$dd(paste0(safe(num("percent_features_hit"), min, 1), "%"))),
      tags$div(tags$dt("Highest genes-hit sample"), tags$dd(paste0(safe(num("percent_features_hit"), max, 1), "%")))
    )
  )
}

summary_qc_role_value <- function(values) {
  values <- tolower(paste(stats::na.omit(as.character(values)), collapse = " "))
  if (!nzchar(values)) return(NA_character_)
  if (grepl("\\btreated\\b", values)) return("treated")
  if (grepl("\\bcontrol\\b|\\bmock\\b|\\bparent\\b|classifier_control", values)) return("control")
  NA_character_
}

summary_qc_project_sample_roles <- function(project) {
  metadata <- if (is.data.frame(project$samples) && nrow(project$samples)) project$samples else data.frame()
  summary_root <- file.path(project$project_root, "summary")
  summary_samples <- if (dir.exists(summary_root)) {
    unique(basename(dirname(list.files(summary_root, pattern = "\\.all_hits\\.csv$|\\.feature_table\\.RDF_1\\.csv$", recursive = TRUE, full.names = TRUE))))
  } else character()
  if (is.data.frame(metadata) && "sample" %in% names(metadata)) {
    samples <- unique(as.character(metadata$sample))
  } else {
    samples <- summary_samples
  }
  if (length(summary_samples)) samples <- intersect(samples, summary_samples)
  samples <- samples[nzchar(samples) & !is.na(samples)]
  if (!length(samples)) return(data.frame(sample = character(), role = character(), stringsAsFactors = FALSE))
  role_cols <- c("library_role", "control_or_treated", "fitness_role", "classifier_role",
                 "condition", "condition_label", "treatment", "treatment_label",
                 "guessed_condition")
  roles <- vapply(samples, function(sample) {
    row <- if (is.data.frame(metadata) && "sample" %in% names(metadata)) metadata[as.character(metadata$sample) == sample, , drop = FALSE] else data.frame()
    if (nrow(row)) {
      for (col in intersect(role_cols, names(row))) {
        role <- summary_qc_role_value(row[[col]])
        if (!is.na(role)) return(role)
      }
    }
    summary_qc_role_value(sample)
  }, character(1))
  included <- rep(TRUE, length(samples))
  if (is.data.frame(metadata) && "sample" %in% names(metadata)) {
    idx <- match(samples, as.character(metadata$sample))
    include_cols <- intersect(c("include", "included"), names(metadata))
    include_col <- if (length(include_cols)) include_cols[[1]] else ""
    if (nzchar(include_col)) {
      include_values <- tolower(as.character(metadata[[include_col]][idx]))
      included <- is.na(include_values) | include_values %in% c("true", "1", "yes", "")
    }
  }
  data.frame(sample = samples, role = roles, included = included, stringsAsFactors = FALSE)
}

summary_qc_read_feature_hits <- function(file, sample, role) {
  data <- tryCatch(read.csv(file, skip = 1, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(data) || !"standard_name" %in% names(data) || !"hits" %in% names(data)) return(NULL)
  feature <- as.character(data$standard_name)
  keep <- !is.na(feature) & nzchar(feature)
  if (!any(keep)) return(NULL)
  data.frame(
    sample = sample,
    role = role,
    feature_id = feature[keep],
    hits = suppressWarnings(as.numeric(data$hits[keep])),
    stringsAsFactors = FALSE
  )
}

summary_qc_read_all_hit_sites <- function(file, sample, role) {
  data <- tryCatch(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(data) || ncol(data) < 2L || !nrow(data)) return(NULL)
  chrom_cols <- intersect(c("Chromosome", "chromosome", "chrom", "Chrom", "seqnames"), names(data))
  pos_cols <- intersect(c("Position", "position", "pos", "Hit position", "hit_position", "start"), names(data))
  chrom_col <- if (length(chrom_cols)) chrom_cols[[1]] else names(data)[[1]]
  pos_col <- if (length(pos_cols)) pos_cols[[1]] else names(data)[[2]]
  chrom <- as.character(data[[chrom_col]])
  pos <- as.character(data[[pos_col]])
  keep <- !is.na(chrom) & nzchar(chrom) & !is.na(pos) & nzchar(pos)
  if (!any(keep)) return(NULL)
  data.frame(
    sample = sample,
    role = role,
    insertion_site = paste(chrom[keep], pos[keep], sep = ":"),
    stringsAsFactors = FALSE
  )
}

summary_qc_read_hit_count_total <- function(project, sample) {
  stats_file <- file.path(project$project_root, "summary", sample, "stats.csv")
  if (file.exists(stats_file)) {
    stats <- tryCatch(read.csv(stats_file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.data.frame(stats) && nrow(stats) && "Total Hits" %in% names(stats)) {
      total_hits <- suppressWarnings(as.numeric(stats[["Total Hits"]][[1]]))
      if (is.finite(total_hits)) return(total_hits)
    }
  }
  candidates <- c(
    file.path(project$project_root, "create_hit_file", sample, paste0(sample, "_hits.txt")),
    file.path(project$project_root, "create_hit_file", sample, paste0(sample, ".hits.txt"))
  )
  metadata <- if (is.data.frame(project$samples) && "sample" %in% names(project$samples)) project$samples else data.frame()
  if (nrow(metadata)) {
    row <- metadata[as.character(metadata$sample) == sample, , drop = FALSE]
    if (nrow(row) && "hit_file" %in% names(row)) {
      hit_file <- as.character(row$hit_file[[1]] %||% "")
      if (nzchar(hit_file)) {
        candidates <- c(candidates, hit_file, file.path(project$repo_root %||% "", hit_file))
      }
    }
  }
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  hit_file <- if (length(existing)) existing[[1]] else ""
  if (!nzchar(hit_file)) return(NA_real_)
  data <- tryCatch(read.delim(hit_file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(data) || !nrow(data)) return(NA_real_)
  nrow(data)
}

summary_qc_insertion_role_stats_raw <- function(project) {
  roles <- summary_qc_project_sample_roles(project)
  roles <- roles[roles$included & roles$role %in% c("control", "treated"), , drop = FALSE]
  if (!nrow(roles)) return(data.frame())
  summary_root <- file.path(project$project_root, "summary")
  feature_rows <- lapply(seq_len(nrow(roles)), function(i) {
    sample <- roles$sample[[i]]
    file <- file.path(summary_root, sample, paste0(sample, ".feature_table.RDF_1.csv"))
    if (!file.exists(file)) return(NULL)
    summary_qc_read_feature_hits(file, sample, roles$role[[i]])
  })
  site_rows <- lapply(seq_len(nrow(roles)), function(i) {
    sample <- roles$sample[[i]]
    file <- file.path(summary_root, sample, paste0(sample, ".all_hits.csv"))
    if (!file.exists(file)) return(NULL)
    summary_qc_read_all_hit_sites(file, sample, roles$role[[i]])
  })
  feature_data <- do.call(rbind, Filter(Negate(is.null), feature_rows))
  site_data <- do.call(rbind, Filter(Negate(is.null), site_rows))
  role_labels <- c(control = "Parents / controls", treated = "Treated")
  out <- lapply(c("control", "treated"), function(role) {
    role_samples <- roles$sample[roles$role == role]
    features <- if (is.data.frame(feature_data) && nrow(feature_data)) feature_data[feature_data$role == role, , drop = FALSE] else data.frame()
    sites <- if (is.data.frame(site_data) && nrow(site_data)) site_data[site_data$role == role, , drop = FALSE] else data.frame()
    total_genes <- if (nrow(features)) length(unique(features$feature_id)) else NA_integer_
    feature_totals <- if (nrow(features)) aggregate(hits ~ feature_id, features, sum, na.rm = TRUE) else data.frame()
    hit_genes <- if (nrow(feature_totals)) sum(feature_totals$hits > 0, na.rm = TRUE) else NA_integer_
    feature_hit_total <- if (nrow(feature_totals)) sum(feature_totals$hits, na.rm = TRUE) else NA_real_
    unique_insertions <- if (nrow(sites)) length(unique(sites$insertion_site)) else NA_integer_
    data.frame(
      role = role,
      `Sample group` = role_labels[[role]],
      Samples = length(role_samples),
      `Total unique insertions` = unique_insertions,
      `Total genes evaluated` = total_genes,
      `Genes with insertions` = hit_genes,
      `Average insertions / gene` = if (is.finite(total_genes) && total_genes > 0) feature_hit_total / total_genes else NA_real_,
      `Average insertions / hit gene` = if (is.finite(hit_genes) && hit_genes > 0) feature_hit_total / hit_genes else NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, out)
  result[result$Samples > 0, , drop = FALSE]
}

summary_qc_insertion_role_stats <- function(project) {
  data <- summary_qc_insertion_role_stats_raw(project)
  if (!nrow(data)) return(data)
  data <- data[, setdiff(names(data), "role"), drop = FALSE]
  data <- format_qc_metrics(data, count_columns = c("Samples", "Total unique insertions", "Total genes evaluated", "Genes with insertions"),
                            mean_columns = c("Average insertions / gene", "Average insertions / hit gene"))
  data
}

summary_qc_sample_label <- function(sample, role = "", pool = "", project = NULL) {
  alias <- ytab_sample_labels(sample, project)
  if (!identical(alias, as.character(sample))) return(alias)
  role_label <- if (identical(role, "control")) "Parent/control" else if (identical(role, "treated")) "Treated" else "Sample"
  pool <- as.character(pool %||% "")
  if (nzchar(pool) && !identical(pool, "NA")) return(paste(role_label, paste0("pool ", pool)))
  paste(role_label, sample)
}

summary_qc_pool_value <- function(metadata, sample, cols) {
  if (!is.data.frame(metadata) || !"sample" %in% names(metadata)) return("")
  row <- metadata[as.character(metadata$sample) == sample, , drop = FALSE]
  if (!nrow(row)) return("")
  for (col in intersect(cols, names(row))) {
    value <- as.character(row[[col]][[1]] %||% "")
    if (nzchar(value) && !is.na(value)) return(sub("^pool", "", value, ignore.case = TRUE))
  }
  ""
}

summary_qc_pool_insertion_stats_raw <- function(project) {
  roles <- summary_qc_project_sample_roles(project)
  roles <- roles[roles$included & roles$role %in% c("control", "treated"), , drop = FALSE]
  if (!nrow(roles)) return(data.frame())
  metadata <- if (is.data.frame(project$samples) && nrow(project$samples)) project$samples else data.frame()
  summary_root <- file.path(project$project_root, "summary")
  rows <- lapply(seq_len(nrow(roles)), function(i) {
    sample <- roles$sample[[i]]
    role <- roles$role[[i]]
    file <- file.path(summary_root, sample, paste0(sample, ".all_hits.csv"))
    if (!file.exists(file)) return(NULL)
    sites <- summary_qc_read_all_hit_sites(file, sample, role)
    if (is.null(sites) || !nrow(sites)) return(NULL)
    pool <- summary_qc_pool_value(metadata, sample, c("pool", "pool_id", "guessed_pool", "replicate_id"))
    role_label <- if (identical(role, "control")) "Parents / controls" else "Treated"
    data.frame(
      role = role,
      `Sample group` = role_label,
      sample = sample,
      Pool = pool,
      `Display sample` = summary_qc_sample_label(sample, role, pool, project),
      `Total hit count` = summary_qc_read_hit_count_total(project, sample),
      `Unique insertion sites` = length(unique(sites$insertion_site)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  data <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(data) || !nrow(data)) return(data.frame())
  role_order <- c(control = 1L, treated = 2L)
  pool_number <- suppressWarnings(as.numeric(data$Pool))
  data[order(role_order[data$role], pool_number, data$sample, na.last = TRUE), , drop = FALSE]
}

summary_qc_pool_insertion_stats <- function(project) {
  data <- summary_qc_pool_insertion_stats_raw(project)
  if (!nrow(data)) return(data)
  data <- data[, setdiff(names(data), "role"), drop = FALSE]
  format_qc_metrics(data, count_columns = c("Total hit count", "Unique insertion sites"))
}

plot_qc_summary_insertion_pool_stats <- function(project, roles = c("control", "treated"),
                                                 samples = NULL, metric = "unique") {
  data <- summary_qc_pool_insertion_stats_raw(project)
  if (!nrow(data)) return(qc_plot_empty("Per-pool insertion data are not available."))
  roles <- roles %||% c("control", "treated")
  roles <- intersect(as.character(roles), c("control", "treated"))
  if (!length(roles)) return(qc_plot_empty("Select at least one sample group."))
  data <- data[data$role %in% roles, , drop = FALSE]
  if (!is.null(samples) && length(samples)) {
    samples <- as.character(samples)
    data <- data[data$sample %in% samples, , drop = FALSE]
  }
  if (!nrow(data)) return(qc_plot_empty("No selected pools are available."))
  metric <- if (identical(metric, "total")) "total" else "unique"
  value_col <- switch(metric,
                      total = "Total hit count",
                      "Unique insertion sites")
  values <- as.numeric(data[[value_col]])
  labels <- as.character(data$`Display sample`)
  fills <- c(control = "#d8e4ef", treated = "#2f6f9f")[as.character(data$role)]
  fills[is.na(fills)] <- qc_plot_fill
  horizontal <- qc_plot_bar_horizontal(labels) || length(labels) > 4L
  old <- par(no.readonly = TRUE)
  on.exit({ qc_plot_reset_layout(); par(old) }, add = TRUE)
  qc_plot_reset_layout()
  qc_plot_par(mar = qc_plot_bar_margins(labels, horizontal, top = 4.0, right = if (horizontal) 5.6 else 2.8))
  limit <- max(values, na.rm = TRUE)
  if (!is.finite(limit) || limit <= 0) limit <- 1
  at <- barplot(values, names.arg = rep("", length(values)), axes = FALSE, horiz = horizontal,
                xlim = if (horizontal) c(0, limit * 1.14) else NULL,
                ylim = if (horizontal) NULL else c(0, limit * 1.16),
                col = fills, border = qc_plot_border, lwd = max(1.7, qc_plot_lwd),
                xlab = if (horizontal) value_col else "",
                ylab = if (horizontal) "" else value_col,
                main = paste(value_col, "per pool"))
  qc_plot_add_grid(horizontal)
  axis_values <- pretty(c(0, limit))
  axis_labels <- format(axis_values, big.mark = ",", scientific = FALSE, trim = TRUE)
  if (horizontal) axis(1, at = axis_values, labels = axis_labels) else axis(2, at = axis_values, labels = axis_labels)
  if (horizontal) qc_plot_draw_horizontal_labels(at, labels)
  else qc_plot_draw_vertical_labels(at, labels, angle = qc_plot_label_angle(30L))
  if (qc_plot_show_value_labels()) {
    value_labels <- format(round(values), big.mark = ",", trim = TRUE)
    text(if (horizontal) values else at,
         if (horizontal) at else values,
         labels = value_labels, pos = if (horizontal) 4 else 3, xpd = NA,
         font = 2, cex = qc_plot_value_cex(0.88))
  }
}

summary_qc_gene_hit_bins_raw <- function(project, roles = c("control", "treated"),
                                         samples = NULL, mode = "per_pool") {
  role_data <- summary_qc_project_sample_roles(project)
  role_data <- role_data[role_data$included & role_data$role %in% c("control", "treated"), , drop = FALSE]
  roles <- intersect(as.character(roles %||% c("control", "treated")), c("control", "treated"))
  if (length(roles)) role_data <- role_data[role_data$role %in% roles, , drop = FALSE]
  if (!is.null(samples) && length(samples)) {
    role_data <- role_data[role_data$sample %in% as.character(samples), , drop = FALSE]
  }
  if (!nrow(role_data)) return(data.frame())
  metadata <- if (is.data.frame(project$samples) && nrow(project$samples)) project$samples else data.frame()
  summary_root <- file.path(project$project_root, "summary")
  feature_rows <- lapply(seq_len(nrow(role_data)), function(i) {
    sample <- role_data$sample[[i]]
    role <- role_data$role[[i]]
    file <- file.path(summary_root, sample, paste0(sample, ".feature_table.RDF_1.csv"))
    if (!file.exists(file)) return(NULL)
    data <- summary_qc_read_feature_hits(file, sample, role)
    if (is.null(data) || !nrow(data)) return(NULL)
    pool <- summary_qc_pool_value(metadata, sample, c("pool", "pool_id", "guessed_pool", "replicate_id"))
    data$Pool <- pool
    data$`Sample group` <- if (identical(role, "control")) "Parents / controls" else "Treated"
    data$`Display sample` <- summary_qc_sample_label(sample, role, pool, project)
    data
  })
  features <- do.call(rbind, Filter(Negate(is.null), feature_rows))
  if (is.null(features) || !nrow(features)) return(data.frame())
  features$hits[is.na(features$hits)] <- 0
  mode <- if (identical(mode, "combined")) "combined" else "per_pool"
  # Insertion sites per gene: 0, 1-3, 4-7, 8-9, 10-19, >=20. Wider upper bins keep dense libraries from piling up in the top bin.
  bin_levels <- c("0", "1\u20133", "4\u20137", "8\u20139", "10\u201319", "\u226520")
  make_bins <- function(data, group_label, role = "", sample = "", sample_group = "") {
    bins <- as.character(cut(floor(data$hits), c(-Inf, 0, 3, 7, 9, 19, Inf), labels = bin_levels))
    counts <- as.data.frame(table(factor(bins, levels = bin_levels)), stringsAsFactors = FALSE)
    names(counts) <- c("Hit bin", "Genes")
    counts$Genes <- as.integer(counts$Genes)
    counts$group_label <- group_label
    counts$role <- role
    counts$sample <- sample
    counts$`Sample group` <- sample_group
    counts
  }
  if (identical(mode, "combined")) {
    rows <- lapply(split(features, features$role), function(group_data) {
      role <- unique(group_data$role)[[1]]
      label <- if (identical(role, "control")) "Parents / controls combined" else "Treated combined"
      summed <- aggregate(hits ~ feature_id, group_data, sum, na.rm = TRUE)
      make_bins(summed, label, role = role, sample = "", sample_group = sub(" combined$", "", label))
    })
  } else {
    rows <- lapply(split(features, features$sample), function(sample_data) {
      make_bins(sample_data, unique(sample_data$`Display sample`)[[1]],
                role = unique(sample_data$role)[[1]],
                sample = unique(sample_data$sample)[[1]],
                sample_group = unique(sample_data$`Sample group`)[[1]])
    })
  }
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out)) return(data.frame())
  out$`Hit bin` <- factor(out$`Hit bin`, levels = bin_levels)
  role_order <- c(control = 1L, treated = 2L)
  out[order(role_order[out$role], out$group_label, out$`Hit bin`, na.last = TRUE), , drop = FALSE]
}

plot_qc_summary_gene_hit_bins <- function(project, roles = c("control", "treated"),
                                          samples = NULL, mode = "per_pool") {
  data <- summary_qc_gene_hit_bins_raw(project, roles = roles, samples = samples, mode = mode)
  if (!nrow(data)) return(qc_plot_empty("Gene insertion-count distribution is not available."))
  groups <- unique(as.character(data$group_label))
  bins <- levels(data$`Hit bin`)
  matrix_data <- matrix(0, nrow = length(bins), ncol = length(groups),
                        dimnames = list(bins, groups))
  for (i in seq_len(nrow(data))) {
    matrix_data[as.character(data$`Hit bin`[[i]]), as.character(data$group_label[[i]])] <- data$Genes[[i]]
  }
  display_groups <- qc_plot_display_labels(groups)
  horizontal <- qc_plot_bar_horizontal(display_groups) || length(groups) > 4L
  colors <- c("0" = "#F7FBFF", "1" = "#D0E1F2", "2" = "#9CC8E8",
              "3" = "#5DA5D8", "4" = "#1F78B4", "\u22655" = "#08306B")
  old <- par(no.readonly = TRUE)
  on.exit({ qc_plot_reset_layout(); par(old) }, add = TRUE)
  qc_plot_begin_key_layout(0.14)
  vertical_label_angle <- if (!horizontal && length(groups) <= 3L) 0L else qc_plot_label_angle(30L)
  plot_margins <- if (horizontal) {
    qc_plot_bar_margins(display_groups, TRUE, top = 4.0, right = 5.4)
  } else {
    c(qc_plot_label_margin_lines(display_groups, angle = vertical_label_angle, orientation = "vertical") + 1.6,
      7.4, 4.0, 2.8)
  }
  qc_plot_par(mar = plot_margins)
  totals <- colSums(matrix_data, na.rm = TRUE)
  limit <- max(totals, na.rm = TRUE)
  if (!is.finite(limit) || limit <= 0) limit <- 1
  bar_width <- 0.42
  bar_space <- 0.62
  at_preview <- barplot(matrix_data, plot = FALSE, horiz = horizontal,
                        width = bar_width, space = bar_space)
  minimum_category_span <- if (identical(mode, "combined")) 3.4 else 5.8
  padded_category_limits <- function(at, minimum_span, edge_pad = 0.75) {
    upper <- max(max(at, na.rm = TRUE) + edge_pad, minimum_span)
    c(0, upper)
  }
  if (horizontal) {
    ylim <- padded_category_limits(at_preview, minimum_category_span)
    at <- barplot(matrix_data, names.arg = rep("", length(groups)), axes = FALSE,
                  horiz = TRUE, col = colors[rownames(matrix_data)],
                  border = "#1f2d3a", lwd = max(1.4, qc_plot_lwd * 0.85),
                  width = bar_width, space = bar_space,
                  xlim = c(0, limit * 1.14), ylim = ylim,
                  xlab = "", ylab = "",
                  main = "Genes by insertion count")
  } else {
    xlim <- padded_category_limits(at_preview, minimum_category_span)
    at <- barplot(matrix_data, names.arg = rep("", length(groups)), axes = FALSE,
                  horiz = FALSE, col = colors[rownames(matrix_data)],
                  border = "#1f2d3a", lwd = max(1.4, qc_plot_lwd * 0.85),
                  width = bar_width, space = bar_space,
                  xlim = xlim, ylim = c(0, limit * 1.16),
                  xlab = "", ylab = "",
                  main = "Genes by insertion count")
  }
  qc_plot_add_grid(horizontal)
  axis_values <- pretty(c(0, limit))
  axis_labels <- format(axis_values, big.mark = ",", scientific = FALSE, trim = TRUE)
  if (horizontal) {
    axis(1, at = axis_values, labels = axis_labels, line = 0.45)
    mtext("Number of genes", side = 1, line = 3.4, font = 2, cex = qc_plot_text_sizes()$lab)
  } else {
    axis(2, at = axis_values, labels = axis_labels, las = 1)
    mtext("Number of genes", side = 2, line = 4.2, font = 2, cex = qc_plot_text_sizes()$lab)
  }
  if (horizontal) qc_plot_draw_horizontal_labels(at, groups)
  else qc_plot_draw_vertical_labels(at, groups, angle = vertical_label_angle)
  if (qc_plot_show_value_labels()) {
    text(if (horizontal) totals else at,
         if (horizontal) at else totals,
         labels = format(totals, big.mark = ",", trim = TRUE),
         pos = if (horizontal) 4 else 3, xpd = NA,
         font = 2, cex = qc_plot_value_cex(0.82))
  }
  qc_plot_metric_key_row(rownames(matrix_data), fill = colors[rownames(matrix_data)],
                         border = qc_plot_border, ncol = 6,
                         title = "Hits per gene")
}

plot_qc_summary_insertion_role_stats <- function(project) {
  data <- summary_qc_insertion_role_stats_raw(project)
  if (!nrow(data)) return(qc_plot_empty("Parent/treated insertion summary is not available."))
  data <- data[data$role %in% c("control", "treated"), , drop = FALSE]
  if (!nrow(data)) return(qc_plot_empty("Parent/treated insertion summary is not available."))
  data$role <- factor(data$role, levels = c("control", "treated"))
  data <- data[order(data$role), , drop = FALSE]
  labels <- as.character(data$`Sample group`)
  fills <- c(control = "#d8e4ef", treated = "#2f6f9f")[as.character(data$role)]
  fills[is.na(fills)] <- qc_plot_fill
  show_values <- qc_plot_show_value_labels()
  old <- par(no.readonly = TRUE)
  on.exit({ qc_plot_reset_layout(); par(old) }, add = TRUE)
  qc_plot_reset_layout()
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))

  draw_metric <- function(values, title, ylab, value_digits = 0L, count_axis = FALSE) {
    display_labels <- qc_plot_display_labels(labels)
    old_par <- qc_plot_par(mar = c(qc_plot_label_margin_lines(display_labels, angle = 0L, orientation = "vertical") + 1.0,
                                  8.0, 4.2, 2.2), mgp = c(6.0, 0.9, 0))
    on.exit(par(old_par), add = TRUE)
    values <- as.numeric(values)
    display <- if (count_axis) qc_plot_count_display(values, getOption("ytab.plot.count_scale", "thousands")) else
      list(values = values, divisor = 1, suffix = "")
    plotted <- display$values
    limit <- max(plotted, na.rm = TRUE)
    if (!is.finite(limit) || limit <= 0) limit <- 1
    at <- barplot(plotted, names.arg = rep("", length(values)), axes = FALSE,
                  ylim = c(0, limit * 1.18),
                  col = fills, border = qc_plot_border, lwd = max(1.7, qc_plot_lwd),
                  ylab = paste0(ylab, display$suffix), main = title)
    qc_plot_add_grid(FALSE)
    ticks <- axTicks(2)
    axis(2, at = ticks, labels = qc_plot_axis_labels(ticks / display$divisor), las = 1)
    qc_plot_draw_vertical_labels(at, labels, angle = 0L, cex = qc_plot_label_cex(0.92))
    if (show_values) {
      value_labels <- if (value_digits > 0L) {
        sprintf(paste0("%.", value_digits, "f"), values)
      } else {
        format(round(values), big.mark = ",", trim = TRUE)
      }
      text(at, plotted, labels = value_labels, pos = 3, xpd = NA,
           font = 2, cex = qc_plot_value_cex(0.88))
    }
    box()
  }

  draw_metric(data$`Total unique insertions`,
              "Total unique insertions",
              "Unique insertion sites",
              value_digits = 0L, count_axis = TRUE)
  draw_metric(data$`Average insertions / gene`,
              "Average insertions per gene",
              "Insertions per gene",
              value_digits = 2L)
}

compact_qc_table <- function(data) DT::datatable(data,rownames=FALSE,selection="none",options=list(pageLength=10,lengthMenu=c(10,25,50),autoWidth=FALSE,ordering=TRUE,searching=TRUE,scrollX=TRUE,columnDefs=list(list(className="dt-right",targets=which(vapply(data,is.numeric,FALSE))-1L),list(className="ytab-nowrap",targets=0))))

qc_sample_eligibility <- function(project,status) {
  samples<-project$samples; included<-if("include"%in%names(samples))tolower(as.character(samples$include))%in%c("true","1","yes")else rep(TRUE,nrow(samples));idx<-match(samples$sample,status$sample);hit_status<-status$hit_file[idx];path<-status$hits_path[idx];valid<-included & hit_status=="success" & vapply(path,nonempty_file,FALSE)
  data.frame(Sample=as.character(samples$sample),Condition=as.character(samples$guessed_condition%||%""),Background=as.character(samples$guessed_background%||%""),Pool=as.character(samples$guessed_pool%||%""),`Hit-file status`=ifelse(valid,"Available",ifelse(included,"Unavailable","Excluded")),Eligible=ifelse(valid,"Yes","No"),Reason=ifelse(valid,"",ifelse(!included,"Not included in project","Valid successful hit file unavailable")),hit_path=path,stringsAsFactors=FALSE,check.names=FALSE)
}

diagnostic_cache_manifest <- function(project,run_id) { path<-file.path(project$project_root,"manifests","library_diagnostics","runs",run_id,"manifest.json");if(!file.exists(path))return(NULL);tryCatch(jsonlite::fromJSON(path,simplifyVector=FALSE),error=function(e)NULL) }
diagnostic_action_label <- function(mode="preview",cache_match=FALSE,force=FALSE) if(mode=="preview")"Preview diagnostics"else if(force)"Rerun diagnostics"else if(cache_match)"Use cached diagnostics"else"Run diagnostics"
diagnostic_result_heading <- function(status) switch(status,preview="Preview complete",cached="Cached result reused",success="Library Diagnostics complete",failed="Library Diagnostics failed","Library Diagnostics")

human_file_size <- function(bytes) {bytes<-as.numeric(bytes);vapply(bytes,function(x)if(is.na(x))"Unknown"else if(x<1024)sprintf("%d B",round(x))else if(x<1024^2)sprintf("%.0f KB",x/1024)else if(x<1024^3)sprintf("%.1f MB",x/1024^2)else sprintf("%.1f GB",x/1024^3),"")}
diagnostic_plot_type <- function(filename) {value<-tolower(filename);if(grepl("centromere_bias",value))"Centromere bias"else if(grepl("metaplots",value))"Feature metaplots"else if(grepl("tss_metaplot",value))"TSS metaplot"else if(grepl("trna_metaplot",value))"tRNA metaplot"else if(grepl("seqbias",value))"Sequence bias"else if(grepl("midlc",value))"MidLC"else"Diagnostic plot"}
diagnostic_plot_title <- function(filename) diagnostic_plot_type(filename)

diagnostic_resource_roots <- function(project) {
  roots <- c(
    `ytab-diagnostics-project` = file.path(project$project_root, "library_diagnostics"),
    `ytab-diagnostics-export` = file.path(project$export_root, "qc", "diagnostics")
  )
  roots[dir.exists(roots)]
}

diagnostic_path_within <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  startsWith(path, paste0(root, "/"))
}

diagnostic_encode_relative_path <- function(path) {
  parts <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1L]]
  if (!length(parts) || any(parts %in% c("", ".", ".."))) return("")
  paste(vapply(parts, URLencode, character(1), reserved = FALSE), collapse = "/")
}

diagnostic_served_file <- function(path, project) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  display <- relative_project_path(path, project$project_root)
  result <- list(served_url = "", preview_available = FALSE,
                 relative_path = display,
                 preview_message = "Preview unavailable; file exists but could not be served.")
  if (!file.exists(path) || dir.exists(path)) {
    result$preview_message <- "Preview unavailable; file is missing."
    return(result)
  }
  roots <- diagnostic_resource_roots(project)
  for (route in names(roots)) {
    root <- roots[[route]]
    if (!diagnostic_path_within(path, root)) next
    relative <- substring(path, nchar(normalizePath(root, winslash = "/", mustWork = FALSE)) + 2L)
    encoded <- diagnostic_encode_relative_path(relative)
    if (!nzchar(encoded)) return(result)
    result$served_url <- paste0(route, "/", encoded)
    result$preview_available <- TRUE
    result$relative_path <- if (identical(route, "ytab-diagnostics-project"))
      file.path("library_diagnostics", relative) else
      file.path("exports", "qc", "diagnostics", relative)
    return(result)
  }
  result
}

build_diagnostic_file_inventory <- function(project) {
  roots<-c(file.path(project$project_root,"library_diagnostics"),file.path(project$export_root,"qc","diagnostics"),file.path(project$project_root,"manifests","library_diagnostics"),file.path(project$project_root,"logs","library_diagnostics"));files<-unique(unlist(lapply(roots,function(root)if(dir.exists(root))list.files(root,recursive=TRUE,full.names=TRUE)else character())));files<-files[file.exists(files)&!dir.exists(files)];if(!length(files))return(data.frame())
  served<-lapply(files,diagnostic_served_file,project=project);filename<-basename(files);ext<-tolower(tools::file_ext(files));is_plot<-ext%in%c("png","jpg","jpeg","svg");type<-ifelse(is_plot,"Plot",ifelse(ext%in%c("csv","tsv"),"Table",ifelse(ext=="json","Manifest",ifelse(ext=="html","HTML",ifelse(ext%in%c("txt","log"),"Text","Log")))));relative<-vapply(served,`[[`,"","relative_path");parts<-strsplit(relative,"/");sample<-vapply(parts,function(x)if(length(x)>2&&x[[1]]%in%c("library_diagnostics","exports")&&!x[[length(x)-1]]%in%c("runs","library_diagnostics","diagnostics"))x[[length(x)-1]]else"Project","");run<-vapply(parts,function(x){index<-match("runs",x);if(!is.na(index)&&length(x)>index)x[[index+1L]]else"current_legacy"},"");info<-file.info(files)
  data.frame(file=files,filename=filename,extension=ext,display_type=type,plot_type=vapply(filename,diagnostic_plot_type,""),display_title=vapply(filename,diagnostic_plot_title,""),sample=sample,sample_set=ifelse(run=="current_legacy","Legacy diagnostics result",tools::toTitleCase(gsub("_"," ",run))),run_id=run,size_bytes=as.numeric(info$size),size_display=human_file_size(info$size),modified=format(info$mtime,"%Y-%m-%d %H:%M"),relative_path=relative,served_url=vapply(served,`[[`,"","served_url"),preview_available=vapply(served,`[[`,FALSE,"preview_available"),preview_message=vapply(served,`[[`,"","preview_message"),viewable=ext%in%c("png","jpg","jpeg","svg","csv","tsv","txt","html"),is_plot=is_plot,stringsAsFactors=FALSE)
}

filter_diagnostic_inventory <- function(data,sample="All",plot_type="All",sample_set="All",run_id="All",type="All") {
  if(!nrow(data))return(data);keep<-rep(TRUE,nrow(data));if(sample!="All")keep<-keep&data$sample==sample;if(plot_type!="All")keep<-keep&data$plot_type==plot_type;if(sample_set!="All")keep<-keep&data$sample_set==sample_set;if(run_id!="All")keep<-keep&data$run_id==run_id;if(type!="All")keep<-keep&data$display_type==type;data[keep,,drop=FALSE]
}
paginate_diagnostic_inventory <- function(data,page=1L,page_size=8L) {pages<-max(1L,ceiling(nrow(data)/page_size));page<-max(1L,min(as.integer(page),pages));start<-(page-1L)*page_size+1L;rows<-if(nrow(data)&&start<=nrow(data))seq.int(start,min(start+page_size-1L,nrow(data)))else integer();list(data=data[rows,,drop=FALSE],page=page,pages=pages,total=nrow(data))}

# The Mapping QC, Summary QC and Library Diagnostics panels were replaced by the Quality Control plot menus (qc_module.R);
# the data functions above feed those plots.
# quality_control_ui() is in qc_module.R (Overview table + themed subtabs).
