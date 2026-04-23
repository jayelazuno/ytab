#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# plot_feature_summary_qc.R
#
# Purpose:
#   Read per-sample stats.csv files from:
#     output/smoketests/summary/<sample>/stats.csv
#
#   combine selected samples,
#   build one clean feature-summary / classifier-style QC figure,
#   and save a PNG into:
#     output/exports/<export_version>/qc/images/
#
# Notes:
#   - Library complexity was handled separately, so this script does NOT plot
#     Total Hits.
#   - This script focuses on feature-summary QC from stats.csv:
#       1) % of genomic features hit
#       2) Mean reads per hit
#       3) Feature vs intergenic hit distribution
#
# Usage examples:
#   Rscript scripts/install_plot_packages.R
#
#   Rscript scripts/plot_feature_summary_qc.R \
#     --project_root=.
#
#   Rscript scripts/plot_feature_summary_qc.R \
#     --project_root=. \
#     --samples=yH298-parent-pool1,yH298-parent-pool2
#
#   Rscript scripts/plot_feature_summary_qc.R \
#     --project_root=. \
#     --export_version=smoke_test_v1 \
#     --figure_name=feature_summary_qc.png
# -----------------------------------------------------------------------------

required_pkgs <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "tidyr",
  "stringr",
  "purrr",
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
  library(scales)
  library(patchwork)
})

parse_args <- function(args) {
  defaults <- list(
    project_root   = ".",
    export_version = "smoke_test_v1",
    samples        = NULL,
    figure_name    = "feature_summary_qc.png"
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

canonical_sample_order <- function(samples) {
  base_order <- c(
    "yH298-parent-pool1",
    "yH298-parent-pool2",
    "yH299-parent-pool3",
    "yH299-parent-pool4",
    "yH298-H2O2-treated-facs-pool1",
    "yH298-H2O2-treated-facs-pool2",
    "yH299-H2O2-treated-facs-pool3",
    "yH299-H2O2-treated-facs-pool4"
  )

  extras <- setdiff(samples, base_order)
  c(base_order, sort(extras)) |> unique() |> intersect(samples)
}

find_stats_files <- function(summary_dir, selected_samples = NULL) {
  if (!dir.exists(summary_dir)) {
    stop("Summary directory was not found: ", summary_dir, call. = FALSE)
  }

  stats_files <- list.files(
    summary_dir,
    pattern = "^stats\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  stats_files <- stats_files[dirname(stats_files) != summary_dir]

  if (length(stats_files) == 0) {
    stop("No stats.csv files were found under: ", summary_dir, call. = FALSE)
  }

  file_tbl <- tibble(
    path = stats_files,
    sample = basename(dirname(stats_files))
  ) %>%
    distinct(sample, .keep_all = TRUE)

  if (!is.null(selected_samples) && nzchar(selected_samples)) {
    wanted <- str_split(selected_samples, pattern = ",", simplify = FALSE)[[1]] %>%
      str_trim() %>%
      keep(~ .x != "")

    missing <- setdiff(wanted, file_tbl$sample)
    if (length(missing) > 0) {
      stop(
        "Selected sample(s) not found in output/smoketests/summary: ",
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
    ord <- canonical_sample_order(file_tbl$sample)
    file_tbl <- file_tbl %>%
      mutate(sample = factor(sample, levels = ord)) %>%
      arrange(sample) %>%
      mutate(sample = as.character(sample))
  }

  file_tbl
}

read_stats_data <- function(file_tbl) {
  stats_tbl <- map_dfr(seq_len(nrow(file_tbl)), function(i) {
    sample_i <- file_tbl$sample[i]
    path_i <- file_tbl$path[i]

    dat <- read_csv(path_i, show_col_types = FALSE) %>%
      mutate(sample = sample_i)

    dat
  })

  required_cols <- c(
    "% of hits in features",
    "% of intergenic hits",
    "% of features hit",
    "Mean Reads Per Hit"
  )

  missing_cols <- setdiff(required_cols, names(stats_tbl))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "stats.csv is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  percent_cols <- c(
    "% of hits in features",
    "% of intergenic hits",
    "% of features hit"
  )

  stats_tbl %>%
    mutate(
      condition = if_else(
        str_detect(sample, "parent"),
        "parent",
        "H2O2-treated-facs"
      ),
      condition = factor(condition, levels = c("parent", "H2O2-treated-facs")),
      sample = factor(sample, levels = file_tbl$sample)
    ) %>%
    mutate(
      across(
        all_of(percent_cols),
        ~ round(as.numeric(str_remove(as.character(.x), "%")), 1)
      ),
      `Mean Reads Per Hit` = round(as.numeric(`Mean Reads Per Hit`), 2)
    )
}

make_qc_theme <- function() {
  theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust = 1,
        size = 10, color = "black", face = "bold"
      ),
      axis.text.y = element_text(size = 10, color = "black", face = "bold"),
      axis.title  = element_text(face = "bold", size = 13),
      strip.text  = element_text(face = "bold", size = 11),
      legend.title = element_blank(),
      legend.text  = element_text(size = 11, color = "black", face = "bold"),
      axis.ticks.length = unit(2.5, "mm"),
      panel.grid.minor = element_blank()
    )
}

make_feature_summary_plot <- function(stats_tbl) {
  qc_theme <- make_qc_theme()

  p_features_hit <- ggplot(
    stats_tbl,
    aes(x = sample, y = `% of features hit`)
  ) +
    geom_col(
      width = 0.72,
      fill = "grey96",
      color = "black",
      linewidth = 0.9
    ) +
    facet_wrap(~ condition, scales = "free_x") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = "% genomic features hit",
      x = NULL,
      y = "% features hit"
    ) +
    qc_theme +
    theme(legend.position = "none")

  p_reads_per_hit <- ggplot(
    stats_tbl,
    aes(x = sample, y = `Mean Reads Per Hit`)
  ) +
    geom_col(
      width = 0.72,
      fill = "grey96",
      color = "black",
      linewidth = 0.9
    ) +
    facet_wrap(~ condition, scales = "free_x") +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    labs(
      title = "Mean reads per insertion site",
      x = NULL,
      y = "Mean reads per hit"
    ) +
    qc_theme +
    theme(legend.position = "none")

  stats_long <- stats_tbl %>%
    select(
      sample, condition,
      `% of hits in features`,
      `% of intergenic hits`
    ) %>%
    pivot_longer(
      cols = starts_with("%"),
      names_to = "category",
      values_to = "percent"
    ) %>%
    mutate(
      category = recode(
        category,
        `% of hits in features` = "Features",
        `% of intergenic hits` = "Intergenic"
      ),
      category = factor(category, levels = c("Features", "Intergenic"))
    )

  p_feature_intergenic <- ggplot(
    stats_long,
    aes(x = sample, y = percent, fill = category)
  ) +
    geom_col(
      width = 0.72,
      color = "black",
      linewidth = 0.9
    ) +
    facet_wrap(~ condition, scales = "free_x") +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = "Feature vs intergenic hit distribution",
      x = NULL,
      y = "% of hits"
    ) +
    qc_theme

  p_features_hit / p_reads_per_hit / p_feature_intergenic +
    plot_annotation(
      title = "Feature-summary / classifier-style QC",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
      )
    )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  project_root <- normalizePath(args$project_root, winslash = "/", mustWork = FALSE)
  summary_dir <- file.path(project_root, "output", "smoketests", "summary")
  export_img_dir <- file.path(project_root, "output", "exports", args$export_version, "qc", "images")
  dir.create(export_img_dir, recursive = TRUE, showWarnings = FALSE)

  file_tbl <- find_stats_files(summary_dir, args$samples)
  stats_tbl <- read_stats_data(file_tbl)

  if (nrow(stats_tbl) == 0) {
    stop("No feature-summary QC data could be read after filtering.", call. = FALSE)
  }

  p <- make_feature_summary_plot(stats_tbl)
  out_file <- file.path(export_img_dir, args$figure_name)

  n_samples <- dplyr::n_distinct(stats_tbl$sample)

  ggsave(
    filename = out_file,
    plot = p,
    width = max(12, 0.95 * n_samples + 6),
    height = 14,
    dpi = 300,
    bg = "white"
  )

  message("Saved feature-summary QC figure: ", out_file)
}

main()
