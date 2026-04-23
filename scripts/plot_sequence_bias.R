#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# plot_sequence_bias.R
#
# Purpose:
#   Read per-sample Hermes sequence-bias files from:
#     output/smoketests/library_diagnostics/<sample>/<sample>.seqbias_2_7.tsv
#   combine selected samples,
#   make one clean figure focused on the +2 / +7 Hermes insertion signature,
#   and save a PNG into:
#     output/exports/<export_version>/qc/images/
#
# Usage examples:
#   Rscript scripts/install_plot_packages.R
#
#   Rscript scripts/plot_sequence_bias.R \
#     --project_root=. \
#     --export_version=smoke_test_v1
#
#   Rscript scripts/plot_sequence_bias.R \
#     --project_root=. \
#     --export_version=smoke_test_v1 \
#     --samples=yH298-parent-pool1,yH298-parent-pool2
# -----------------------------------------------------------------------------

required_pkgs <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "stringr",
  "purrr",
  "forcats",
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
  library(forcats)
  library(scales)
})

parse_args <- function(args) {
  defaults <- list(
    project_root   = ".",
    export_version = "smoke_test_v1",
    samples        = NULL,
    figure_name    = "sequence_bias_signature.png"
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

find_seqbias_files <- function(diagnostics_dir, selected_samples = NULL) {
  if (!dir.exists(diagnostics_dir)) {
    stop("Library diagnostics directory was not found: ", diagnostics_dir, call. = FALSE)
  }

  all_files <- list.files(
    diagnostics_dir,
    pattern = "\\.seqbias_2_7\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(all_files) == 0) {
    stop("No .seqbias_2_7.tsv files were found under: ", diagnostics_dir, call. = FALSE)
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

standardize_seqbias_cols <- function(df) {
  rename_pairs <- c(
    motif = "pair",
    dinucleotide = "pair",
    kmer = "pair",
    dimer = "pair",
    observed_reads = "obs_reads",
    reads = "obs_reads",
    observed_fraction = "obs_frac",
    observed_freq = "obs_frac",
    observed_prop = "obs_frac",
    background_sites = "bg_sites",
    genome_sites = "bg_sites",
    background_fraction = "bg_frac",
    background_freq = "bg_frac",
    background_prop = "bg_frac",
    fold_enrichment = "enrichment",
    enrichment_ratio = "enrichment",
    obs_over_bg = "enrichment",
    ratio = "enrichment"
  )

  for (old_name in names(rename_pairs)) {
    new_name <- rename_pairs[[old_name]]
    if (old_name %in% names(df) && !new_name %in% names(df)) {
      names(df)[names(df) == old_name] <- new_name
    }
  }

  df
}

read_one_seqbias <- function(path_i, sample_i) {
  raw <- read_tsv(path_i, show_col_types = FALSE) %>%
    standardize_seqbias_cols()

  needed <- c("pair", "enrichment")
  missing <- setdiff(needed, names(raw))
  if (length(missing) > 0) {
    stop(
      paste0(
        "The sequence-bias file is missing required columns: ",
        paste(missing, collapse = ", "),
        "\nFile: ", path_i,
        "\nExpected columns include: pair, obs_reads, obs_frac, bg_sites, bg_frac, enrichment"
      ),
      call. = FALSE
    )
  }

  raw %>%
    transmute(
      sample = sample_i,
      pair = str_to_upper(as.character(pair)),
      obs_reads = if ("obs_reads" %in% names(raw)) as.numeric(obs_reads) else NA_real_,
      obs_frac  = if ("obs_frac"  %in% names(raw)) as.numeric(obs_frac)  else NA_real_,
      bg_sites  = if ("bg_sites"  %in% names(raw)) as.numeric(bg_sites)  else NA_real_,
      bg_frac   = if ("bg_frac"   %in% names(raw)) as.numeric(bg_frac)   else NA_real_,
      enrichment = as.numeric(enrichment)
    ) %>%
    filter(!is.na(pair), !is.na(enrichment), nchar(pair) == 2) %>%
    mutate(
      base_p2 = substr(pair, 1, 1),
      base_p7 = substr(pair, 2, 2),
      condition = if_else(str_detect(sample, "parent"), "parent", "treated")
    )
}

read_seqbias_data <- function(file_tbl) {
  pair_levels <- as.vector(outer(c("A", "C", "G", "T"), c("A", "C", "G", "T"), paste0))

  map_dfr(seq_len(nrow(file_tbl)), function(i) {
    read_one_seqbias(file_tbl$path[i], file_tbl$sample[i])
  }) %>%
    mutate(
      sample = factor(sample, levels = file_tbl$sample),
      pair = factor(pair, levels = pair_levels),
      condition = factor(condition, levels = c("parent", "treated"))
    ) %>%
    arrange(sample, pair)
}

make_sequence_bias_plot <- function(seq_tbl) {
  ggplot(seq_tbl, aes(x = pair, y = enrichment, fill = condition)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.35) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
    facet_wrap(~ sample, ncol = 2, scales = "fixed") +
    labs(
      title = "Hermes sequence bias signature (+2 / +7)",
      subtitle = "Each 2-base pair encodes the nucleotide at +2 followed by the nucleotide at +7 relative to the insertion site",
      x = "Dinucleotide pair (+2, +7)",
      y = "Observed / background enrichment",
      caption = "Dashed line marks enrichment = 1 (no bias)."
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 9),
      axis.text.y = element_text(face = "bold", size = 10),
      strip.text = element_text(face = "bold", size = 11),
      legend.title = element_blank(),
      legend.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  project_root <- normalizePath(args$project_root, winslash = "/", mustWork = FALSE)
  diagnostics_dir <- file.path(project_root, "output", "smoketests", "library_diagnostics")
  export_img_dir <- file.path(project_root, "output", "exports", args$export_version, "qc", "images")
  dir.create(export_img_dir, recursive = TRUE, showWarnings = FALSE)

  file_tbl <- find_seqbias_files(diagnostics_dir, args$samples)
  seq_tbl <- read_seqbias_data(file_tbl)
  p <- make_sequence_bias_plot(seq_tbl)

  out_file <- file.path(export_img_dir, args$figure_name)
  ggsave(
    filename = out_file,
    plot = p,
    width = 12,
    height = max(6, ceiling(nlevels(seq_tbl$sample) / 2) * 3.6),
    dpi = 300,
    bg = "white"
  )

  message("Saved sequence bias figure: ", out_file)
}

main()
