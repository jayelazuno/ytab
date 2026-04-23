#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# plot_gene_feature_metaplots.R
#
# Purpose:
#   Read per-sample gene-feature metaplot files from:
#     output/smoketests/library_diagnostics/<sample>/<sample>.tss_metaplot.tsv
#     output/smoketests/library_diagnostics/<sample>/<sample>.tts_metaplot.tsv
#     output/smoketests/library_diagnostics/<sample>/<sample>.trna_metaplot.tsv
#
#   combine selected samples,
#   make one clean metaplot figure for TSS, TTS, and tRNA profiles,
#   and save a PNG into:
#     output/exports/<export_version>/qc/images/
#
# Expected input columns in each .tsv file:
#   rel_bp, reads_sum, sites_n
#
# Plotting note:
#   The plotted signal is mean reads per site in each relative-position bin:
#     mean_reads_per_site = reads_sum / sites_n
#   This is more comparable across positions than raw reads_sum alone.
#
# Usage examples:
#   Rscript scripts/install_plot_packages.R
#
#   Rscript scripts/plot_gene_feature_metaplots.R \
#     --project_root=.
#
#   Rscript scripts/plot_gene_feature_metaplots.R \
#     --project_root=. \
#     --export_version=smoke_test_v1
#
#   Rscript scripts/plot_gene_feature_metaplots.R \
#     --project_root=. \
#     --export_version=smoke_test_v1 \
#     --samples=yH298-parent-pool1,yH298-parent-pool2
# -----------------------------------------------------------------------------

required_pkgs <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "tidyr",
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
  library(tidyr)
  library(stringr)
  library(purrr)
  library(scales)
})

parse_args <- function(args) {
  defaults <- list(
    project_root   = ".",
    export_version = "smoke_test_v1",
    samples        = NULL,
    figure_name    = "gene_feature_metaplots.png"
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

find_feature_files <- function(diagnostics_dir, selected_samples = NULL) {
  if (!dir.exists(diagnostics_dir)) {
    stop("Library diagnostics directory was not found: ", diagnostics_dir, call. = FALSE)
  }

  feature_patterns <- c(
    tss = "\\.tss_metaplot\\.tsv$",
    tts = "\\.tts_metaplot\\.tsv$",
    trna = "\\.trna_metaplot\\.tsv$"
  )

  tbls <- lapply(names(feature_patterns), function(feature_name) {
    files <- list.files(
      diagnostics_dir,
      pattern = feature_patterns[[feature_name]],
      recursive = TRUE,
      full.names = TRUE
    )

    if (length(files) == 0) {
      return(tibble(sample = character(), feature = character(), path = character()))
    }

    tibble(
      path = files,
      sample = basename(dirname(files)),
      feature = feature_name
    )
  })

  file_tbl <- bind_rows(tbls) %>%
    distinct(sample, feature, .keep_all = TRUE)

  if (nrow(file_tbl) == 0) {
    stop(
      "No gene-feature metaplot files were found under: ", diagnostics_dir,
      "\nExpected files end with .tss_metaplot.tsv, .tts_metaplot.tsv, or .trna_metaplot.tsv",
      call. = FALSE
    )
  }

  if (!is.null(selected_samples) && nzchar(selected_samples)) {
    wanted <- str_split(selected_samples, pattern = ",", simplify = FALSE)[[1]] %>%
      str_trim() %>%
      discard(~ .x == "")

    missing <- setdiff(wanted, unique(file_tbl$sample))
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
      arrange(sample, feature) %>%
      mutate(sample = as.character(sample))
  } else {
    ord <- canonical_sample_order(unique(file_tbl$sample))
    file_tbl <- file_tbl %>%
      mutate(sample = factor(sample, levels = ord)) %>%
      arrange(sample, feature) %>%
      mutate(sample = as.character(sample))
  }

  missing_combos <- expand_grid(
    sample = unique(file_tbl$sample),
    feature = c("tss", "tts", "trna")
  ) %>%
    anti_join(file_tbl, by = c("sample", "feature"))

  if (nrow(missing_combos) > 0) {
    warning(
      "Some sample-feature metaplot files are missing and will be omitted:\n",
      paste0("  - ", missing_combos$sample, " (", missing_combos$feature, ")", collapse = "\n"),
      call. = FALSE
    )
  }

  file_tbl
}

standardize_metaplot_cols <- function(df) {
  rename_pairs <- c(
    relative_bp = "rel_bp",
    relative_pos = "rel_bp",
    distance_bp = "rel_bp",
    position = "rel_bp",
    reads = "reads_sum",
    reads_total = "reads_sum",
    total_reads = "reads_sum",
    insertions_reads = "reads_sum",
    n_sites = "sites_n",
    sites = "sites_n",
    number_sites = "sites_n",
    windows_n = "sites_n"
  )

  for (old_name in names(rename_pairs)) {
    new_name <- rename_pairs[[old_name]]
    if (old_name %in% names(df) && !new_name %in% names(df)) {
      names(df)[names(df) == old_name] <- new_name
    }
  }

  df
}

pretty_feature_name <- function(x) {
  recode(x,
    tss  = "TSS",
    tts  = "TTS",
    trna = "tRNA"
  )
}

read_one_metaplot <- function(path_i, sample_i, feature_i) {
  raw <- read_tsv(path_i, show_col_types = FALSE) %>%
    standardize_metaplot_cols()

  needed <- c("rel_bp", "reads_sum", "sites_n")
  missing <- setdiff(needed, names(raw))
  if (length(missing) > 0) {
    stop(
      paste0(
        "The metaplot file is missing required columns: ",
        paste(missing, collapse = ", "),
        "\nFile: ", path_i,
        "\nExpected columns include: rel_bp, reads_sum, sites_n"
      ),
      call. = FALSE
    )
  }

  raw %>%
    transmute(
      sample = sample_i,
      feature = feature_i,
      rel_bp = as.numeric(rel_bp),
      reads_sum = as.numeric(reads_sum),
      sites_n = as.numeric(sites_n)
    ) %>%
    filter(!is.na(rel_bp), !is.na(reads_sum), !is.na(sites_n)) %>%
    mutate(
      mean_reads_per_site = if_else(sites_n > 0, reads_sum / sites_n, NA_real_),
      feature_label = pretty_feature_name(feature)
    )
}

read_metaplot_data <- function(file_tbl) {
  map_dfr(seq_len(nrow(file_tbl)), function(i) {
    read_one_metaplot(
      path_i = file_tbl$path[i],
      sample_i = file_tbl$sample[i],
      feature_i = file_tbl$feature[i]
    )
  }) %>%
    mutate(
      sample = factor(sample, levels = unique(file_tbl$sample)),
      feature = factor(feature, levels = c("tss", "tts", "trna")),
      feature_label = factor(feature_label, levels = c("TSS", "TTS", "tRNA")),
      condition = if_else(str_detect(as.character(sample), "parent"), "parent", "treated")
    ) %>%
    arrange(sample, feature, rel_bp)
}

make_metaplot_figure <- function(meta_tbl) {
  max_abs_bp <- suppressWarnings(max(abs(meta_tbl$rel_bp), na.rm = TRUE))
  if (!is.finite(max_abs_bp)) max_abs_bp <- 1000

  x_breaks <- pretty(c(-max_abs_bp, max_abs_bp), n = 5)
  x_breaks <- unique(c(x_breaks, 0))
  x_breaks <- x_breaks[x_breaks >= -max_abs_bp & x_breaks <= max_abs_bp]

  ggplot(meta_tbl, aes(x = rel_bp, y = mean_reads_per_site, group = sample)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 0.9, color = "black") +
    facet_grid(feature_label ~ sample, scales = "free_y") +
    scale_x_continuous(
      breaks = x_breaks,
      labels = label_number(big.mark = ",")
    ) +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    labs(
      title = "Gene-feature metaplots",
      subtitle = "Mean insertion-associated reads per site around TSS, TTS, and tRNA features",
      x = "Position relative to feature (bp)",
      y = "Mean reads per site",
      caption = "Dashed vertical line marks the annotated feature position (0 bp)."
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 9),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.7, "lines")
    )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  project_root <- normalizePath(args$project_root, winslash = "/", mustWork = FALSE)
  diagnostics_dir <- file.path(project_root, "output", "smoketests", "library_diagnostics")
  export_img_dir <- file.path(project_root, "output", "exports", args$export_version, "qc", "images")
  dir.create(export_img_dir, recursive = TRUE, showWarnings = FALSE)

  file_tbl <- find_feature_files(diagnostics_dir, args$samples)
  meta_tbl <- read_metaplot_data(file_tbl)

  if (nrow(meta_tbl) == 0) {
    stop("No metaplot data could be read after filtering.", call. = FALSE)
  }

  p <- make_metaplot_figure(meta_tbl)

  n_samples <- dplyr::n_distinct(meta_tbl$sample)
  out_file <- file.path(export_img_dir, args$figure_name)

  ggsave(
    filename = out_file,
    plot = p,
    width = max(10, 2.9 * n_samples),
    height = max(7.5, 2.5 * 3),
    dpi = 300,
    bg = "white"
  )

  message("Saved gene-feature metaplot figure: ", out_file)
}

main()
