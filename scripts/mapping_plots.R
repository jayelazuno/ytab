
#!/usr/bin/env Rscript
## usage 
# Rscript scripts/mapping_plots.R --project_root=. 
#selected samples:
# Rscript scripts/plot_mapping_qc.R \
 # --project_root=. \
  #--samples=yH298-parent-pool1,yH299-parent-pool3
# change the export version or output PNG name:
 # Rscript scripts/plot_mapping_qc.R \
  #--project_root=. \
  #--export_version=smoke_test_v1 \
  #--figure_name=mapping_qc_selected.png



required_pkgs <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "tidyr",
  "stringr",
  "purrr",
  "forcats",
  "scales",
  "patchwork"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      "\nRun scripts/install_plot_packages.R first."
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(forcats)
  library(scales)
  library(patchwork)
})

parse_args <- function(args) {
  defaults <- list(
    project_root   = ".",
    export_version = "smoke_test_v1",
    samples        = NULL,
    figure_name    = "mapping_qc_overview.png"
  )

  if (length(args) == 0) return(defaults)

  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    key <- sub("^--", "", strsplit(arg, "=", fixed = TRUE)[[1]][1])
    value <- sub("^[^=]+=", "", arg)
    if (key %in% names(defaults)) defaults[[key]] <- value
  }

  defaults
}

standardize_mapping_cols <- function(df) {
  rename_pairs <- c(
    total_reads = "total_records",
    primary_total = "primary_records",
    mapped_reads = "primary_mapped",
    unmapped_reads = "primary_unmapped",
    hq_reads = "mapq_ge20",
    pct_mapped = "percent_mapped",
    pct_mapq_ge20 = "percent_mapq_ge20",
    pct_duplicates = "percent_duplicates",
    mean_mapq = "avg_mapq_mapped_primary"
  )

  for (old_name in names(rename_pairs)) {
    new_name <- rename_pairs[[old_name]]
    if (old_name %in% names(df) && !new_name %in% names(df)) {
      names(df)[names(df) == old_name] <- new_name
    }
  }

  df
}

find_mapping_files <- function(smoketest_dir, selected_samples = NULL) {
  if (!dir.exists(smoketest_dir)) {
    stop("Smoketest mapfastq directory was not found: ", smoketest_dir, call. = FALSE)
  }

  all_files <- list.files(
    smoketest_dir,
    pattern = "\\.mapping_stats\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(all_files) == 0) {
    stop("No mapping stats files were found under: ", smoketest_dir, call. = FALSE)
  }

  file_tbl <- tibble(
    path = all_files,
    sample = basename(dirname(all_files))
  ) %>%
    distinct(sample, .keep_all = TRUE)

  if (!is.null(selected_samples) && nzchar(selected_samples)) {
    wanted <- str_split(selected_samples, pattern = ",", simplify = FALSE)[[1]] %>%
      str_trim() %>%
      discard(~ .x == "")

    missing <- setdiff(wanted, file_tbl$sample)
    if (length(missing) > 0) {
      stop(
        "Selected sample(s) not found in output/smoketests/mapfastq: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }

    file_tbl <- file_tbl %>%
      filter(sample %in% wanted) %>%
      mutate(sample = factor(sample, levels = wanted)) %>%
      arrange(sample) %>%
      mutate(sample = as.character(sample))
  } else {
    file_tbl <- file_tbl %>% arrange(sample)
  }

  file_tbl
}

read_mapping_stats <- function(file_tbl) {
  map_dfr(seq_len(nrow(file_tbl)), function(i) {
    path_i <- file_tbl$path[i]
    sample_i <- file_tbl$sample[i]

    read_csv(path_i, show_col_types = FALSE) %>%
      standardize_mapping_cols() %>%
      mutate(sample = sample_i, .before = 1)
  }) %>%
    distinct(sample, .keep_all = TRUE)
}

compute_missing_metrics <- function(df) {
  out <- df

  if (!"primary_records" %in% names(out) && "total_records" %in% names(out)) {
    out$primary_records <- out$total_records
  }

  if (!"percent_mapped" %in% names(out) && all(c("primary_mapped", "primary_records") %in% names(out))) {
    out$percent_mapped <- ifelse(out$primary_records > 0, 100 * out$primary_mapped / out$primary_records, NA_real_)
  }

  if (!"percent_duplicates" %in% names(out) && all(c("primary_duplicates", "primary_mapped") %in% names(out))) {
    out$percent_duplicates <- ifelse(out$primary_mapped > 0, 100 * out$primary_duplicates / out$primary_mapped, NA_real_)
  }

  if (!"percent_mapq_ge20" %in% names(out) && all(c("mapq_ge20", "primary_mapped") %in% names(out))) {
    out$percent_mapq_ge20 <- ifelse(out$primary_mapped > 0, 100 * out$mapq_ge20 / out$primary_mapped, NA_real_)
  }

  out
}

make_total_reads_panel <- function(df) {
  if (!"total_records" %in% names(df)) return(NULL)

  df_reads <- df %>%
    mutate(sample = factor(sample, levels = sample))

  ggplot(df_reads, aes(x = sample, y = total_records)) +
    geom_col(width = 0.72, fill = "grey55") +
    coord_flip() +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    labs(
      title = "Total reads",
      x = NULL,
      y = "Reads"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )
}

make_percent_panel <- function(df) {
  metric_map <- c(
    percent_mapped = "% mapped",
    percent_mapq_ge20 = "% HQ (MAPQ≥20)",
    percent_duplicates = "% duplicates"
  )

  keep_metrics <- names(metric_map)[names(metric_map) %in% names(df)]
  if (length(keep_metrics) == 0) return(NULL)

  long_df <- df %>%
    select(sample, all_of(keep_metrics)) %>%
    pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
    mutate(
      metric = recode(metric, !!!metric_map),
      sample = factor(sample, levels = unique(df$sample))
    )

  p <- ggplot(long_df, aes(x = sample, y = value, color = metric, group = metric))

  if (dplyr::n_distinct(long_df$sample) > 1) {
    p <- p + geom_line(linewidth = 0.8, alpha = 0.75)
  }

  p +
    geom_point(size = 3) +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, 100),
      labels = label_number(suffix = "%")
    ) +
    labs(
      title = "Mapping and alignment rates",
      x = NULL,
      y = "Percent",
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 11),
      legend.position = "top",
      legend.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )
}

make_mapq_panel <- function(df) {
  if (!"avg_mapq_mapped_primary" %in% names(df)) return(NULL)

  df_mapq <- df %>%
    mutate(sample = factor(sample, levels = sample))

  ggplot(df_mapq, aes(x = sample, y = avg_mapq_mapped_primary)) +
    geom_col(width = 0.72, fill = "grey35") +
    coord_flip() +
    labs(
      title = "Average MAPQ of mapped primary reads",
      x = NULL,
      y = "Average MAPQ"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  project_root <- normalizePath(args$project_root, winslash = "/", mustWork = FALSE)
  smoketest_dir <- file.path(project_root, "output", "smoketests", "mapfastq")
  export_img_dir <- file.path(project_root, "output", "exports", args$export_version, "qc", "images")
  dir.create(export_img_dir, recursive = TRUE, showWarnings = FALSE)

  file_tbl <- find_mapping_files(smoketest_dir, args$samples)
  mapping_stats <- read_mapping_stats(file_tbl) %>%
    compute_missing_metrics()

  mapping_stats <- mapping_stats %>%
    mutate(sample = factor(sample, levels = file_tbl$sample)) %>%
    arrange(sample) %>%
    mutate(sample = as.character(sample))

  panels <- compact(list(
    make_total_reads_panel(mapping_stats),
    make_percent_panel(mapping_stats),
    make_mapq_panel(mapping_stats)
  ))

  if (length(panels) == 0) {
    stop("No plottable mapping QC columns were found in the mapping stats files.", call. = FALSE)
  }

  combined_plot <- wrap_plots(panels, ncol = 1, guides = "collect") +
    plot_annotation(
      title = "Mapping QC",
      subtitle = paste0("Samples: ", paste(file_tbl$sample, collapse = ", "))
    ) &
    theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 11)
    )

  out_file <- file.path(export_img_dir, args$figure_name)
  ggsave(
    filename = out_file,
    plot = combined_plot,
    width = 11,
    height = max(6.5, 2.8 * length(panels) + 1.2),
    dpi = 300,
    bg = "white"
  )

  message("Saved Mapping QC figure: ", out_file)
}

main()


