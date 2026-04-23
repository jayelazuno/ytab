#!/usr/bin/env Rscript

#Rscript scripts/install_plot_packages.R

#Rscript scripts/plot_centromere_bias.R \
  #--project_root=. \
  #--export_version=smoke_test_v1

#Rscript scripts/plot_centromere_bias.R \
  #--project_root=. \
  #--export_version=smoke_test_v1 \
  #--samples=yH298-parent-pool1,yH298-parent-pool2

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
    figure_name    = "centromere_bias_overview.png"
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

standardize_centromere_cols <- function(df) {
  rename_pairs <- c(
    bin_start = "bin_start_bp",
    start_bp = "bin_start_bp",
    distance_bp = "bin_start_bp",
    mean_reads = "mean_reads_across_arms",
    mean_signal = "mean_reads_across_arms",
    mean_hits = "mean_reads_across_arms",
    arms = "n_arms"
  )

  for (old_name in names(rename_pairs)) {
    new_name <- rename_pairs[[old_name]]
    if (old_name %in% names(df) && !new_name %in% names(df)) {
      names(df)[names(df) == old_name] <- new_name
    }
  }

  df
}

find_centromere_files <- function(diagnostics_dir, selected_samples = NULL) {
  if (!dir.exists(diagnostics_dir)) {
    stop("Library diagnostics directory was not found: ", diagnostics_dir, call. = FALSE)
  }

  all_files <- list.files(
    diagnostics_dir,
    pattern = "\\.centromere_bins\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(all_files) == 0) {
    stop("No .centromere_bins.tsv files were found under: ", diagnostics_dir, call. = FALSE)
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
        "Selected sample(s) not found in output/smoketests/library_diagnostics: ",
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

read_centromere_data <- function(file_tbl) {
  out <- map_dfr(seq_len(nrow(file_tbl)), function(i) {
    path_i <- file_tbl$path[i]
    sample_i <- file_tbl$sample[i]

    read_tsv(path_i, show_col_types = FALSE) %>%
      standardize_centromere_cols() %>%
      mutate(sample = sample_i, .before = 1)
  })

  required_cols <- c("sample", "bin_start_bp", "mean_reads_across_arms")
  missing_cols <- setdiff(required_cols, names(out))
  if (length(missing_cols) > 0) {
    stop(
      "The centromere files are missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  out %>%
    mutate(
      sample = factor(sample, levels = file_tbl$sample),
      bin_start_bp = as.numeric(bin_start_bp),
      mean_reads_across_arms = as.numeric(mean_reads_across_arms),
      distance_kb = bin_start_bp / 1000
    ) %>%
    arrange(sample, bin_start_bp)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  project_root <- normalizePath(args$project_root, winslash = "/", mustWork = FALSE)
  diagnostics_dir <- file.path(project_root, "output", "smoketests", "library_diagnostics")
  export_img_dir <- file.path(project_root, "output", "exports", args$export_version, "qc", "images")
  dir.create(export_img_dir, recursive = TRUE, showWarnings = FALSE)

  file_tbl <- find_centromere_files(diagnostics_dir, args$samples)
  cent_tbl <- read_centromere_data(file_tbl)

  max_y <- suppressWarnings(max(cent_tbl$mean_reads_across_arms, na.rm = TRUE))
  if (!is.finite(max_y)) max_y <- 1

  p <- ggplot(cent_tbl, aes(x = distance_kb, y = mean_reads_across_arms, group = sample)) +
    geom_line(linewidth = 0.9, color = "black") +
    geom_point(size = 1.4, color = "black") +
    facet_wrap(~ sample, scales = "free_y") +
    scale_x_continuous(labels = label_number(suffix = " kb")) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.08))) +
    labs(
      title = "Centromere-proximal insertion bias",
      subtitle = "Mean reads across chromosome arms as a function of distance from centromeres",
      x = "Distance from centromere",
      y = "Mean reads across arms"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 10),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )

  out_file <- file.path(export_img_dir, args$figure_name)
  ggsave(
    filename = out_file,
    plot = p,
    width = max(10, 3.4 * min(nrow(file_tbl), 4)),
    height = max(5.5, 2.8 * ceiling(nrow(file_tbl) / 2)),
    dpi = 300,
    bg = "white"
  )

  message("Saved centromere bias figure: ", out_file)
}

main()
