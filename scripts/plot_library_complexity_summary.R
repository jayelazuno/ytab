#!/usr/bin/env Rscript

# ============================================================
# YTAB - Library complexity summary QC plot
# ------------------------------------------------------------
# Reads:
#   output/smoketests/library_diagnostics/library_diagnostics.summary.csv
#
# Optional sample filtering via:
#   --samples=sample1,sample2
#
# Saves:
#   output/exports/<export_version>/qc/images/library_complexity_summary.png
# ============================================================

# -----------------------------
# 0) Parse command line args
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- grep(paste0('^', flag, '='), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0('^', flag, '='), '', hit[1])
}

project_root   <- arg_value('--project_root', '.')
export_version <- arg_value('--export_version', 'smoke_test_v1')
figure_name    <- arg_value('--figure_name', 'library_complexity_summary.png')
samples_arg    <- arg_value('--samples', NULL)

selected_samples <- NULL
if (!is.null(samples_arg) && nzchar(samples_arg)) {
  selected_samples <- trimws(strsplit(samples_arg, ',', fixed = TRUE)[[1]])
  selected_samples <- selected_samples[nzchar(selected_samples)]
}

# -----------------------------
# 1) Install + load packages
# -----------------------------
cran_pkgs <- c(
  'ggplot2', 'dplyr', 'tidyr', 'readr', 'stringr',
  'forcats', 'scales', 'fs', 'glue', 'patchwork'
)

to_install <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, repos = 'https://cloud.r-project.org')
}

invisible(lapply(cran_pkgs, library, character.only = TRUE))

# -----------------------------
# 2) Define paths
# -----------------------------
summary_file <- file.path(
  project_root,
  'output', 'smoketests', 'library_diagnostics',
  'library_diagnostics.summary.csv'
)

out_dir <- file.path(
  project_root,
  'output', 'exports', export_version,
  'qc', 'images'
)

fs::dir_create(out_dir)
out_png <- file.path(out_dir, figure_name)

if (!file.exists(summary_file)) {
  stop('Could not find summary file:\n', summary_file)
}

# -----------------------------
# 3) Read + validate input
# -----------------------------
summary_tbl <- readr::read_csv(summary_file, show_col_types = FALSE)

required_cols <- c(
  'sample',
  'total_reads',
  'unique_sites',
  'jackpot_top_frac',
  'jackpot_frac_reads',
  'n_jackpot_sites',
  'midlc_est',
  'depth_ratio_R_over_midlc'
)

missing_cols <- setdiff(required_cols, names(summary_tbl))
if (length(missing_cols) > 0) {
  stop(
    'The summary file is missing required columns: ',
    paste(missing_cols, collapse = ', ')
  )
}

# canonical sample order if present
sample_order <- c(
  'yH298-parent-pool1',
  'yH298-parent-pool2',
  'yH299-parent-pool3',
  'yH299-parent-pool4',
  'yH298-H2O2-treated-facs-pool1',
  'yH298-H2O2-treated-facs-pool2',
  'yH299-H2O2-treated-facs-pool3',
  'yH299-H2O2-treated-facs-pool4'
)

if (!is.null(selected_samples)) {
  missing_samples <- setdiff(selected_samples, summary_tbl$sample)
  if (length(missing_samples) > 0) {
    stop('These selected samples were not found:\n', paste(missing_samples, collapse = '\n'))
  }
  summary_tbl <- dplyr::filter(summary_tbl, sample %in% selected_samples)
}

if (nrow(summary_tbl) == 0) {
  stop('No rows remain after filtering.')
}

extra_samples <- setdiff(unique(summary_tbl$sample), sample_order)
final_order <- c(sample_order, sort(extra_samples))
final_order <- final_order[final_order %in% unique(summary_tbl$sample)]

summary_tbl <- summary_tbl %>%
  dplyr::mutate(
    sample = factor(sample, levels = final_order),
    depth_status = dplyr::case_when(
      is.na(depth_ratio_R_over_midlc) ~ 'NA',
      depth_ratio_R_over_midlc < 0.8 ~ '<0.8',
      depth_ratio_R_over_midlc <= 1.2 ~ '0.8-1.2',
      TRUE ~ '>1.2'
    )
  ) %>%
  dplyr::arrange(sample)

n_samples <- dplyr::n_distinct(summary_tbl$sample)

# -----------------------------
# 4) Plot helpers
# -----------------------------
base_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = 'bold', size = 13),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(face = 'bold'),
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'),
    axis.text.y = element_text(face = 'bold'),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = 'bold')
  )

count_label <- function(x) scales::comma(round(x))
pct_label <- function(x, digits = 1) paste0(format(round(100 * x, digits), nsmall = digits, trim = TRUE), '%')
ratio_label <- function(x, digits = 2) format(round(x, digits), nsmall = digits, trim = TRUE)

# -----------------------------
# 5) Build plots
# -----------------------------
p1 <- ggplot(summary_tbl, aes(x = sample, y = unique_sites)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = count_label(unique_sites)), vjust = -0.25, size = 3.2) +
  labs(
    title = 'Unique insertion sites',
    subtitle = 'Observed unique sites per library',
    x = NULL,
    y = 'Unique sites'
  ) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.02, 0.15))) +
  base_theme

p2 <- ggplot(summary_tbl, aes(x = sample, y = midlc_est)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = count_label(midlc_est)), vjust = -0.25, size = 3.2) +
  labs(
    title = 'MidLC estimate',
    subtitle = 'Estimated read depth at ~50% of max unique sites',
    x = NULL,
    y = 'MidLC estimate'
  ) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.02, 0.15))) +
  base_theme

p3 <- ggplot(summary_tbl, aes(x = sample, y = depth_ratio_R_over_midlc)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 1, linetype = 'dashed', linewidth = 0.6) +
  geom_text(aes(label = ratio_label(depth_ratio_R_over_midlc)), vjust = -0.25, size = 3.2) +
  labs(
    title = 'Depth ratio (R / MidLC)',
    subtitle = 'Sequencing depth relative to MidLC; dashed line = 1',
    x = NULL,
    y = 'Depth ratio'
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  base_theme

p4 <- ggplot(summary_tbl, aes(x = sample, y = jackpot_frac_reads)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = pct_label(jackpot_frac_reads, digits = 1)), vjust = -0.25, size = 3.2) +
  labs(
    title = 'Jackpot fraction of reads',
    subtitle = 'Fraction of total reads captured by jackpot sites',
    x = NULL,
    y = 'Fraction of reads'
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.15))) +
  base_theme

p5 <- ggplot(summary_tbl, aes(x = sample, y = n_jackpot_sites)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = count_label(n_jackpot_sites)), vjust = -0.25, size = 3.2) +
  labs(
    title = 'Number of jackpot sites',
    subtitle = 'Count of insertion sites classified as jackpots',
    x = NULL,
    y = 'Jackpot sites'
  ) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.02, 0.15))) +
  base_theme

p6 <- ggplot(summary_tbl, aes(x = sample, y = jackpot_top_frac)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = pct_label(jackpot_top_frac, digits = 1)), vjust = -0.25, size = 3.2) +
  labs(
    title = 'Top jackpot site fraction',
    subtitle = 'Fraction of reads represented by the single top jackpot site',
    x = NULL,
    y = 'Top-site fraction'
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.15))) +
  base_theme

final_plot <- (p1 + p2 + p3) / (p4 + p5 + p6) +
  patchwork::plot_annotation(
    title = 'Library complexity summary',
    subtitle = 'Companion summary to MidLC: complexity, saturation, and jackpot burden',
    theme = theme(
      plot.title = element_text(face = 'bold', size = 16),
      plot.subtitle = element_text(size = 11)
    )
  )

# -----------------------------
# 6) Save
# -----------------------------
width_out <- max(14, 1.1 * n_samples + 10)
height_out <- 9.5

ggsave(
  filename = out_png,
  plot = final_plot,
  width = width_out,
  height = height_out,
  dpi = 300,
  bg = 'white'
)

message('Library complexity summary figure saved to:\n', out_png)
