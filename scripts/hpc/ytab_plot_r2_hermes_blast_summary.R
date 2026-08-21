#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) default else args[[hit + 1L]]
}
summary_csv <- get_arg("--summary", "")
figures_dir <- get_arg("--figures-dir", "")
if (!nzchar(summary_csv) || !file.exists(summary_csv)) stop("missing --summary CSV", call. = FALSE)
if (!nzchar(figures_dir)) stop("missing --figures-dir", call. = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

x <- read_csv(summary_csv, show_col_types = FALSE)
required <- c("sample_id", "condition_or_role", "pool", "total_r2_reads",
              "hermes_positive_reads", "percent_r2_reads_hermes_positive")
missing <- setdiff(required, names(x))
if (length(missing)) stop("summary CSV missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
x <- x %>%
  mutate(
    percent_r2_reads_hermes_positive = as.numeric(percent_r2_reads_hermes_positive),
    total_r2_reads = as.numeric(total_r2_reads),
    role_label = case_when(
      str_detect(str_to_lower(condition_or_role), "treated") ~ "Zn-treated",
      str_detect(str_to_lower(condition_or_role), "mock|control") ~ "mock/control",
      TRUE ~ condition_or_role
    ),
    sample_label = if_else(nzchar(pool), paste0(pool, " — ", role_label), sample_id),
    sample_sorted = factor(sample_label, levels = sample_label[order(percent_r2_reads_hermes_positive)])
  )

theme_diag <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

p_sample <- ggplot(x, aes(x = sample_sorted, y = percent_r2_reads_hermes_positive, fill = role_label)) +
  geom_col(width = 0.76, color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", percent_r2_reads_hermes_positive)),
            vjust = -0.25, size = 3, check_overlap = TRUE) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Percent of Zn toxicity R2 reads with Hermes CDS BLAST hit",
    x = "Sample",
    y = "R2 reads with Hermes hit (%)",
    fill = "Role"
  ) +
  theme_diag

ggsave(file.path(figures_dir, "r2_hermes_percent_by_sample.png"), p_sample, width = 11, height = 6.5, dpi = 300)
ggsave(file.path(figures_dir, "r2_hermes_percent_by_sample.pdf"), p_sample, width = 11, height = 6.5)

p_role <- ggplot(x, aes(x = role_label, y = percent_r2_reads_hermes_positive, fill = role_label)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.72, color = "grey25") +
  geom_jitter(width = 0.11, size = 2.4, alpha = 0.85) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "Percent of Zn toxicity R2 reads with Hermes CDS BLAST hit by role",
    x = "Role",
    y = "R2 reads with Hermes hit (%)",
    fill = "Role"
  ) +
  theme_diag +
  theme(axis.text.x = element_text(angle = 0), legend.position = "none")

ggsave(file.path(figures_dir, "r2_hermes_percent_by_role.png"), p_role, width = 7, height = 5, dpi = 300)
ggsave(file.path(figures_dir, "r2_hermes_percent_by_role.pdf"), p_role, width = 7, height = 5)

p_depth <- ggplot(x, aes(x = total_r2_reads, y = percent_r2_reads_hermes_positive, color = role_label)) +
  geom_point(size = 2.8, alpha = 0.9) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(
    title = "R2 read depth versus Hermes-positive fraction",
    x = "Total R2 reads analyzed",
    y = "R2 reads with Hermes hit (%)",
    color = "Role"
  ) +
  theme_diag +
  theme(axis.text.x = element_text(angle = 0))

ggsave(file.path(figures_dir, "r2_hermes_depth_vs_percent.png"), p_depth, width = 7, height = 5, dpi = 300)
ggsave(file.path(figures_dir, "r2_hermes_depth_vs_percent.pdf"), p_depth, width = 7, height = 5)

cat("PASS\n")
