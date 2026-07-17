
#!/usr/bin/env Rscript

###########################################################################################
# YTAB sample normalization exploration
#
# Purpose:
#   After running ytab_SampleNormlization.txt, visually explore how different
#   MidLC normalization targets affect each sample.
#
# This script:
#   1) Reads LibraryDiagnostics normalization summaries from:
#        output/exports/<PROJECT_ID>/qc/sample_normalization/
#
#   2) Plots reads retained after normalization across targets.
#
#   3) Compares unnormalized SummaryTable feature tables against normalized
#      SummaryTable feature tables, if those normalized SummaryTable outputs exist.
#
# Important:
#   The correlation section requires normalized feature tables, not only normalized
#   hit files. If those files do not exist yet, this script will still make the
#   normalization-retention plots and will warn that normalized SummaryTable outputs
#   are missing.
#
# Expected unnormalized feature tables:
#   output/projects/<PROJECT_ID>/summary/<sample>/<sample>.feature_table.RDF_1.csv
#
# Expected normalized feature tables:
#   output/projects/<PROJECT_ID>/summary_normalized/T020/<sample>/<sample>_normalized.feature_table.RDF_1.csv
#   output/projects/<PROJECT_ID>/summary_normalized/T030/<sample>/<sample>_normalized.feature_table.RDF_1.csv
#   ...
#
# Outputs:
#   output/exports/<PROJECT_ID>/qc/sample_normalization/explore/tables/
#   output/exports/<PROJECT_ID>/qc/sample_normalization/explore/images/
#
# Usage:
#   Rscript scripts/ytab_plot_sample_normalization.R
#
# Optional:
#   PROJECT_ID=H2O2_screen_v1 Rscript scripts/ytab_plot_sample_normalization.R
#
#   SAMPLES="yH298-parent-pool1 yH298-parent-pool2 yH299-parent-pool3 yH299-parent-pool4" \
#   Rscript scripts/ytab_plot_sample_normalization.R
#
#   NORM_TARGETS="20 30 40 50 60 80 100" \
#   METRICS="reads hits neighborhood_index freedom_index" \
#   Rscript scripts/ytab_plot_sample_normalization.R
###########################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ---------------------------
# 0) Find repo root
# ---------------------------

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }

  return(normalizePath(getwd(), mustWork = TRUE))
}

find_repo_root <- function(start_path) {
  start_dir <- if (dir.exists(start_path)) start_path else dirname(start_path)
  current <- normalizePath(start_dir, mustWork = TRUE)

  repeat {
    marker <- file.path(current, "tnseq.yml")

    if (file.exists(marker)) {
      return(current)
    }

    parent <- dirname(current)

    if (identical(parent, current)) {
      stop(
        "Could not find YTAB repo root. ",
        "Run this script from inside the ytab repo, or keep it inside ytab/scripts/."
      )
    }

    current <- parent
  }
}

script_path <- get_script_path()
ytab_root <- find_repo_root(script_path)

# ---------------------------
# 1) Project settings
# ---------------------------

PROJECT_ID <- Sys.getenv("PROJECT_ID", unset = "H2O2_screen_v1")

# Default to parent/control samples because this normalization is mainly used
# before parent-combined classifier input.
default_samples <- c(
  "yH298-parent-pool1",
  "yH298-parent-pool2",
  "yH299-parent-pool3",
  "yH299-parent-pool4"
)

sample_env <- Sys.getenv("SAMPLES", unset = "")

if (nzchar(sample_env)) {
  samples <- str_split(sample_env, "\\s+")[[1]]
} else {
  samples <- default_samples
}

target_env <- Sys.getenv("NORM_TARGETS", unset = "20 30 40 50 60 80 100")
targets <- as.integer(str_split(target_env, "\\s+")[[1]])

metric_env <- Sys.getenv(
  "METRICS",
  unset = "reads hits insertion_index neighborhood_index freedom_index"
)

metrics <- str_split(metric_env, "\\s+")[[1]]

# ---------------------------
# 2) YTAB paths
# ---------------------------

unnorm_dir <- file.path(
  ytab_root,
  "output",
  "projects",
  PROJECT_ID,
  "summary"
)

norm_summary_dir <- file.path(
  ytab_root,
  "output",
  "exports",
  PROJECT_ID,
  "qc",
  "sample_normalization"
)

norm_feature_dir <- file.path(
  ytab_root,
  "output",
  "projects",
  PROJECT_ID,
  "summary_normalized"
)

explore_root <- file.path(
  norm_summary_dir,
  "explore"
)

plots_dir <- file.path(explore_root, "images")
tables_dir <- file.path(explore_root, "tables")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

message("YTAB root:           ", ytab_root)
message("PROJECT_ID:          ", PROJECT_ID)
message("Unnormalized dir:    ", unnorm_dir)
message("Normalization dir:   ", norm_summary_dir)
message("Norm feature dir:    ", norm_feature_dir)
message("Explore output dir:  ", explore_root)
message("Samples:             ", paste(samples, collapse = ", "))
message("Targets:             ", paste(targets, collapse = ", "))
message("Metrics:             ", paste(metrics, collapse = ", "))

# ---------------------------
# 3) Helpers
# ---------------------------

read_feature_table <- function(f, sample, T = NA_integer_) {
  readr::read_csv(f, skip = 1, show_col_types = FALSE) %>%
    mutate(
      sample = sample,
      normalize_target = T
    )
}

unnorm_file <- function(sample) {
  file.path(
    unnorm_dir,
    sample,
    paste0(sample, ".feature_table.RDF_1.csv")
  )
}

norm_file_candidates <- function(sample, T) {
  tag <- sprintf("T%03d", T)

  c(
    file.path(
      norm_feature_dir,
      tag,
      sample,
      paste0(sample, "_normalized.feature_table.RDF_1.csv")
    ),
    file.path(
      norm_feature_dir,
      tag,
      sample,
      paste0(sample, ".feature_table.RDF_1.csv")
    ),
    file.path(
      norm_feature_dir,
      tag,
      sample,
      paste0(sample, "_analysis.csv")
    )
  )
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

make_gene_vec <- function(df, metric = "hits") {
  if (!metric %in% colnames(df)) {
    stop("Metric not found in feature table: ", metric)
  }

  df %>%
    select(standard_name, value = all_of(metric)) %>%
    mutate(value = as.numeric(value)) %>%
    distinct(standard_name, .keep_all = TRUE)
}

plot_theme <- theme_bw() +
  theme(
    axis.title   = element_text(face = "bold", size = 14),
    axis.text    = element_text(face = "bold", size = 13, color = "black"),
    plot.title   = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12),
    strip.text.x = element_text(face = "bold", size = 13),
    strip.text.y = element_text(face = "bold", size = 13)
  )

# ---------------------------
# 4) Load normalization summaries
# ---------------------------

summary_files <- file.path(
  norm_summary_dir,
  paste0("library_diagnostics.summary.", sprintf("T%03d", targets), ".csv")
)

existing_summary_files <- summary_files[file.exists(summary_files)]

if (length(existing_summary_files) == 0) {
  stop(
    "No normalization summary files found in: ", norm_summary_dir, "\n",
    "Expected files like: library_diagnostics.summary.T020.csv"
  )
}

norm_tbl <- purrr::map_dfr(existing_summary_files, function(f) {
  tag <- str_match(basename(f), "T\\d+")[, 1]
  T_value <- as.integer(str_remove(tag, "^T"))

  readr::read_csv(f, show_col_types = FALSE) %>%
    mutate(
      normalize_target = if ("normalize_target" %in% colnames(.)) {
        as.integer(normalize_target)
      } else {
        T_value
      },
      target_tag = tag
    )
})

write_csv(
  norm_tbl,
  file.path(tables_dir, "sample_normalization_summary.all_targets.csv")
)

# ---------------------------
# 5) Plot reads retained after normalization
# ---------------------------
# The normalization target at depth T is:
#   Rtarget = T x MidLC
# where MidLC is the estimated mid-library complexity.
# Downsampling is applied only when total reads exceed the target.

norm_tbl_plot <- norm_tbl %>%
  filter(sample %in% samples) %>%
  arrange(sample, normalize_target)

needed_norm_cols <- c("sample", "normalize_target", "reads_after_norm")
missing_norm_cols <- setdiff(needed_norm_cols, colnames(norm_tbl_plot))

if (length(missing_norm_cols) > 0) {
  stop(
    "Normalization summary is missing required columns: ",
    paste(missing_norm_cols, collapse = ", ")
  )
}

p_norm_reads <- ggplot(
  norm_tbl_plot,
  aes(
    x = normalize_target,
    y = reads_after_norm,
    group = sample,
    shape = sample
  )
) +
  geom_line(
    linewidth = 0.9,
    color = "black",
    na.rm = TRUE
  ) +
  geom_point(
    size = 3,
    alpha = 0.95,
    color = "black",
    position = position_dodge(width = 1.2),
    na.rm = TRUE
  ) +
  scale_x_continuous(breaks = targets) +
  scale_y_continuous(
    labels = scales::label_number(scale = 1e-6, suffix = "M")
  ) +
  labs(
    title = "Reads retained after MidLC normalization",
    x = "Normalization target (T x MidLC)",
    y = "Reads retained after normalization",
    shape = "Library"
  ) +
  plot_theme

ggsave(
  file.path(plots_dir, "reads_retained_after_normalization.png"),
  p_norm_reads,
  width = 10,
  height = 7,
  dpi = 300
)

p_norm_reads_facet <- p_norm_reads +
  facet_wrap(~ sample, ncol = 1, scales = "free_y") +
  theme(legend.position = "none")

ggsave(
  file.path(plots_dir, "reads_retained_after_normalization.faceted.png"),
  p_norm_reads_facet,
  width = 10,
  height = 14,
  dpi = 300
)

# Optional retained fraction plot if total_reads exists
if ("total_reads" %in% colnames(norm_tbl_plot)) {
  norm_tbl_plot <- norm_tbl_plot %>%
    mutate(frac_reads_retained = reads_after_norm / total_reads)

  write_csv(
    norm_tbl_plot,
    file.path(tables_dir, "sample_normalization_summary.with_fraction_retained.csv")
  )

  p_frac <- ggplot(
    norm_tbl_plot,
    aes(
      x = normalize_target,
      y = frac_reads_retained,
      group = sample,
      shape = sample
    )
  ) +
    geom_line(linewidth = 0.9, color = "black", na.rm = TRUE) +
    geom_point(size = 3, alpha = 0.95, color = "black", na.rm = TRUE) +
    scale_x_continuous(breaks = targets) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    labs(
      title = "Fraction of reads retained after normalization",
      x = "Normalization target (T x MidLC)",
      y = "Reads retained / total reads",
      shape = "Library"
    ) +
    plot_theme

  ggsave(
    file.path(plots_dir, "fraction_reads_retained_after_normalization.png"),
    p_frac,
    width = 10,
    height = 7,
    dpi = 300
  )
}

# ---------------------------
# 6) Load unnormalized feature tables
# ---------------------------

unnorm_tbl <- purrr::map_dfr(samples, function(s) {
  f <- unnorm_file(s)

  if (!file.exists(f)) {
    warning("Missing unnormalized feature table: ", f)
    return(NULL)
  }

  read_feature_table(f, sample = s, T = 0)
})

if (nrow(unnorm_tbl) == 0) {
  warning("No unnormalized SummaryTable feature tables found. Skipping correlation plots.")
  quit(save = "no", status = 0)
}

# ---------------------------
# 7) Compare normalized vs unnormalized feature tables
# ---------------------------
# This requires that normalized hits have already been processed through SummaryTable.
# If normalized feature tables are missing, they are skipped.

all_corr <- list()

for (metric in metrics) {
  message("Computing correlations for metric: ", metric)

  if (!metric %in% colnames(unnorm_tbl)) {
    warning("Skipping metric not found in unnormalized tables: ", metric)
    next
  }

  corr_vs_unnorm <- purrr::map_dfr(samples, function(s) {
    base <- unnorm_tbl %>%
      filter(sample == s, normalize_target == 0) %>%
      make_gene_vec(metric)

    purrr::map_dfr(targets, function(T) {
      f <- first_existing(norm_file_candidates(s, T))

      if (is.na(f)) {
        warning(
          "Missing normalized feature table for sample=", s,
          ", target=T", sprintf("%03d", T),
          ". Skipping."
        )
        return(NULL)
      }

      cur_raw <- read_feature_table(f, sample = s, T = T)

      if (!metric %in% colnames(cur_raw)) {
        warning("Metric ", metric, " not found in normalized file: ", f)
        return(NULL)
      }

      cur <- cur_raw %>%
        make_gene_vec(metric)

      joined <- inner_join(
        base,
        cur,
        by = "standard_name",
        suffix = c("_unnorm", "_norm")
      )

      tibble(
        sample = s,
        metric = metric,
        normalize_target = T,
        n_genes = nrow(joined),
        spearman_rho = suppressWarnings(
          cor(joined$value_unnorm, joined$value_norm, method = "spearman")
        ),
        pearson_r = suppressWarnings(
          cor(joined$value_unnorm, joined$value_norm, method = "pearson")
        )
      )
    })
  })

  if (nrow(corr_vs_unnorm) == 0) {
    warning(
      "No normalized feature tables found for metric ", metric,
      ". Correlation plots not generated for this metric."
    )
    next
  }

  all_corr[[metric]] <- corr_vs_unnorm

  write_csv(
    corr_vs_unnorm,
    file.path(tables_dir, paste0("correlation_vs_unnormalized.", metric, ".csv"))
  )

  p_spearman <- ggplot(
    corr_vs_unnorm,
    aes(x = normalize_target, y = spearman_rho)
  ) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2.4, na.rm = TRUE) +
    facet_wrap(~ sample, ncol = 1) +
    scale_x_continuous(
      limits = c(min(targets), max(targets)),
      breaks = targets
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = paste0("Spearman correlation vs unnormalized (metric: ", metric, ")"),
      x = "Normalization target (T x MidLC)",
      y = "Spearman rho"
    ) +
    plot_theme

  ggsave(
    file.path(plots_dir, paste0("spearman_vs_unnormalized.", metric, ".png")),
    p_spearman,
    width = 10,
    height = 12,
    dpi = 300
  )

  p_pearson <- ggplot(
    corr_vs_unnorm,
    aes(x = normalize_target, y = pearson_r)
  ) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2.4, na.rm = TRUE) +
    facet_wrap(~ sample, ncol = 1) +
    scale_x_continuous(
      limits = c(min(targets), max(targets)),
      breaks = targets
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = paste0("Pearson correlation vs unnormalized (metric: ", metric, ")"),
      x = "Normalization target (T x MidLC)",
      y = "Pearson r"
    ) +
    plot_theme

  ggsave(
    file.path(plots_dir, paste0("pearson_vs_unnormalized.", metric, ".png")),
    p_pearson,
    width = 10,
    height = 12,
    dpi = 300
  )
}

if (length(all_corr) > 0) {
  all_corr_tbl <- bind_rows(all_corr)

  write_csv(
    all_corr_tbl,
    file.path(tables_dir, "correlation_vs_unnormalized.all_metrics.csv")
  )
} else {
  warning(
    "No normalized-vs-unnormalized correlations were generated. ",
    "Most likely, normalized SummaryTable feature tables are not present yet. ",
    "Run SummaryTable on the normalized hits for one or more targets first."
  )
}

message("Done.")
message("Tables: ", tables_dir)
message("Plots:  ", plots_dir)
