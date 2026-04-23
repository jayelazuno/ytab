#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# plot_genome_wide_hit_distribution.R
#
# Purpose:
#   Read per-sample *.all_hits.csv files from:
#     output/smoketests/summary/<sample>/<sample>.all_hits.csv
#
#   bin insertion data by user-selected bin size,
#   summarize either reads or hits,
#   and save genome-wide chromosome-resolved plots into:
#     output/exports/<export_version>/qc/images/
#
# Expected input columns in each .all_hits.csv file:
#   Chromosome, Position, Reads
#
# Flexible options:
#   - users can choose metric: reads or hits
#   - users can choose bin size: e.g. 1000 or 10000
#   - users can choose samples
#   - optionally generate all default combinations at once
#
# Usage examples:
#   Rscript scripts/install_plot_packages.R
#
#   # Make one plot: reads in 1 kb bins for all samples
#   Rscript scripts/plot_genome_wide_hit_distribution.R \
#     --project_root=. \
#     --metric=reads \
#     --bin_size=1000
#
#   # Make one plot: hits in 10 kb bins for selected samples
#   Rscript scripts/plot_genome_wide_hit_distribution.R \
#     --project_root=. \
#     --metric=hits \
#     --bin_size=10000 \
#     --samples=yH298-parent-pool1,yH298-parent-pool2
#
#   # Make all four default plots: reads/hits x 1 kb/10 kb
#   Rscript scripts/plot_genome_wide_hit_distribution.R \
#     --project_root=. \
#     --make_all=true
# -----------------------------------------------------------------------------

required_pkgs <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "stringr",
  "purrr",
  "scales"
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
  library(stringr)
  library(purrr)
  library(scales)
})

parse_args <- function(args) {
  defaults <- list(
    project_root   = ".",
    export_version = "smoke_test_v1",
    samples        = NULL,
    metric         = "reads",
    bin_size       = "1000",
    figure_name    = NULL,
    make_all       = "false"
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

parse_bool <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
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

find_all_hits_files <- function(summary_dir, selected_samples = NULL) {
  if (!dir.exists(summary_dir)) {
    stop("Summary directory was not found: ", summary_dir, call. = FALSE)
  }

  all_hits_files <- list.files(
    summary_dir,
    pattern = "\\.all_hits\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  all_hits_files <- all_hits_files[dirname(all_hits_files) != summary_dir]

  if (length(all_hits_files) == 0) {
    stop(
      "No .all_hits.csv files were found under: ", summary_dir,
      call. = FALSE
    )
  }

  file_tbl <- tibble(
    path = all_hits_files,
    sample = basename(dirname(all_hits_files))
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

read_all_hits_data <- function(file_tbl) {
  hits_tbl <- map_dfr(seq_len(nrow(file_tbl)), function(i) {
    sample_i <- file_tbl$sample[i]
    path_i <- file_tbl$path[i]

    dat <- read_csv(path_i, show_col_types = FALSE) %>%
      mutate(sample = sample_i)

    needed <- c("Chromosome", "Position", "Reads")
    missing <- setdiff(needed, names(dat))
    if (length(missing) > 0) {
      stop(
        paste0(
          "The file is missing required columns: ", paste(missing, collapse = ", "),
          "\nFile: ", path_i,
          "\nExpected columns include: Chromosome, Position, Reads"
        ),
        call. = FALSE
      )
    }

    dat
  })

  chrom_map <- tibble(
    Chromosome = c(
      "CP048230.1","CP048231.1","CP048232.1","CP048233.1",
      "CP048234.1","CP048235.1","CP048236.1","CP048237.1",
      "CP048238.1","CP048239.1","CP048240.1","CP048241.1",
      "CP048242.1"
    ),
    ChromLabel = paste("Chr", LETTERS[1:13])
  )

  hits_tbl <- hits_tbl %>%
    left_join(chrom_map, by = "Chromosome") %>%
    mutate(
      sample = factor(sample, levels = file_tbl$sample),
      Chromosome = if_else(is.na(ChromLabel), as.character(Chromosome), ChromLabel),
      Chromosome = factor(Chromosome, levels = unique(c(chrom_map$ChromLabel, sort(setdiff(unique(Chromosome), chrom_map$ChromLabel))))),
      Position = suppressWarnings(as.numeric(Position)),
      Reads = suppressWarnings(as.numeric(Reads))
    ) %>%
    filter(!is.na(Position), !is.na(Reads)) %>%
    select(sample, Chromosome, Position, Reads)

  hits_tbl
}

bin_hits_data <- function(hits_tbl, bin_size) {
  hits_tbl %>%
    mutate(
      bin_size = as.numeric(bin_size),
      bin = floor(Position / as.numeric(bin_size))
    ) %>%
    group_by(sample, Chromosome, bin_size, bin) %>%
    summarise(
      hits = dplyr::n(),
      reads = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_chrom_bins <- function(df, value_col, bin_size) {
  pretty_metric <- if (value_col == "reads") "Reads" else "Hits"
  pretty_bin <- if (bin_size >= 1000) paste0(format(bin_size / 1000, trim = TRUE), " kb") else paste0(bin_size, " bp")

  ggplot(df, aes(x = bin, y = .data[[value_col]])) +
    geom_col(width = 1, fill = "steelblue") +
    facet_grid(
      Chromosome ~ sample,
      scales = "free_x",
      space = "free_x"
    ) +
    scale_y_log10(labels = label_number()) +
    labs(
      title = paste0(pretty_metric, " per ", pretty_bin, " bin"),
      subtitle = "Genome-wide chromosome-resolved insertion distribution",
      x = "Genomic bin",
      y = paste0("log10(", value_col, " per bin)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(face = "bold"),
      strip.text.x = element_text(angle = 90, face = "bold", size = 10),
      strip.text.y = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.35, "lines")
    )
}

save_one_plot <- function(binned_tbl, metric, bin_size, export_img_dir, figure_name = NULL) {
  metric <- tolower(metric)
  if (!metric %in% c("reads", "hits")) {
    stop("metric must be either 'reads' or 'hits'", call. = FALSE)
  }

  df <- binned_tbl %>% filter(bin_size == as.numeric(bin_size))
  if (nrow(df) == 0) {
    stop("No binned data found for bin_size = ", bin_size, call. = FALSE)
  }

  p <- plot_chrom_bins(df, metric, bin_size)

  if (is.null(figure_name) || !nzchar(figure_name)) {
    figure_name <- paste0("genome_wide_", metric, "_", bin_size, "bp.png")
  }

  out_file <- file.path(export_img_dir, figure_name)
  n_samples <- dplyr::n_distinct(df$sample)
  n_chr <- dplyr::n_distinct(df$Chromosome)

  ggsave(
    filename = out_file,
    plot = p,
    width = max(12, 1.0 * n_samples + 0.55 * n_chr + 10),
    height = max(8, 0.55 * n_chr + 4),
    dpi = 300,
    bg = "white"
  )

  message("Saved genome-wide distribution figure: ", out_file)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  project_root <- normalizePath(args$project_root, winslash = "/", mustWork = FALSE)
  summary_dir <- file.path(project_root, "output", "smoketests", "summary")
  export_img_dir <- file.path(project_root, "output", "exports", args$export_version, "qc", "images")
  dir.create(export_img_dir, recursive = TRUE, showWarnings = FALSE)

  file_tbl <- find_all_hits_files(summary_dir, args$samples)
  hits_tbl <- read_all_hits_data(file_tbl)

  if (nrow(hits_tbl) == 0) {
    stop("No all-hit data could be read after filtering.", call. = FALSE)
  }

  default_bin_sizes <- c(1000, 10000)
  default_metrics <- c("reads", "hits")
  binned_tbl <- map_dfr(default_bin_sizes, ~ bin_hits_data(hits_tbl, .x))

  if (parse_bool(args$make_all)) {
    for (bin_size_i in default_bin_sizes) {
      for (metric_i in default_metrics) {
        save_one_plot(
          binned_tbl = binned_tbl,
          metric = metric_i,
          bin_size = bin_size_i,
          export_img_dir = export_img_dir,
          figure_name = paste0("genome_wide_", metric_i, "_", bin_size_i, "bp.png")
        )
      }
    }
  } else {
    metric_i <- tolower(args$metric)
    bin_size_i <- as.numeric(args$bin_size)

    if (!metric_i %in% c("reads", "hits")) {
      stop("--metric must be 'reads' or 'hits'", call. = FALSE)
    }
    if (!is.finite(bin_size_i) || bin_size_i <= 0) {
      stop("--bin_size must be a positive number", call. = FALSE)
    }

    if (!(bin_size_i %in% default_bin_sizes)) {
      extra_tbl <- bin_hits_data(hits_tbl, bin_size_i)
      binned_tbl <- bind_rows(binned_tbl, extra_tbl)
    }

    save_one_plot(
      binned_tbl = binned_tbl,
      metric = metric_i,
      bin_size = bin_size_i,
      export_img_dir = export_img_dir,
      figure_name = args$figure_name
    )
  }
}

main()
