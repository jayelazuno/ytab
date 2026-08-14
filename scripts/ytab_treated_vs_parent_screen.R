#!/usr/bin/env Rscript

## title: General treated-vs-parent Tn-seq screen analysis for YTAB
## author: Joshua Ayelazuno
## Inputs are explicit raw per-sample SummaryTable files. CPM normalization is
## performed below; MidLC-normalized inputs are intentionally unsupported.

usage <- paste(
  "Usage: Rscript scripts/ytab_treated_vs_parent_screen.R",
  "--project-id ID --comparison-design FILE --output-dir DIR",
  "[--analysis-id ID] [--comparisons ID,ID]",
  "--sample-table SAMPLE=RAW_FEATURE_TABLE [--sample-table ...]",
  "[--pseudocount N] [--z-quantile N] [--top-label-n N]"
)

parse_args <- function(args) {
  result <- list(sample_table = character())
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) { cat(usage, "\n"); quit(status = 0) }
    if (i == length(args) || !startsWith(key, "--")) stop("Invalid arguments.\n", usage)
    value <- args[[i + 1L]]
    name <- gsub("-", "_", substring(key, 3L))
    if (name == "sample_table") result$sample_table <- c(result$sample_table, value) else result[[name]] <- value
    i <- i + 2L
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("project_id", "comparison_design", "output_dir")
missing_args <- required_args[!vapply(required_args, function(x) !is.null(args[[x]]) && nzchar(args[[x]]), logical(1))]
if (length(missing_args)) stop("Missing required arguments: ", paste(missing_args, collapse = ", "), "\n", usage)

required_pkgs <- c(
  "tidyverse",
  "slider",
  "minpack.lm",
  "ggrepel"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall with:\n",
    "install.packages(c(",
    paste(sprintf('\"%s\"', missing_pkgs), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(slider)
  library(minpack.lm)
  library(ggrepel)
})

# ---------------------------
# 0) Project configuration
# ---------------------------

PROJECT_ID <- args$project_id
ANALYSIS_ID <- if (is.null(args$analysis_id)) "treated_vs_parent_raw_cpm" else args$analysis_id
PSEUDOCOUNT <- as.numeric(if (is.null(args$pseudocount)) "0.5" else args$pseudocount)
Z_QUANTILE <- as.numeric(if (is.null(args$z_quantile)) "0.99" else args$z_quantile)
TOP_LABEL_N <- as.integer(if (is.null(args$top_label_n)) "10" else args$top_label_n)
out_root <- normalizePath(args$output_dir, mustWork = FALSE)

table_dir <- file.path(out_root, "tables")
plot_dir <- file.path(out_root, "plots")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(tibble::tibble(
  R_version = R.version.string,
  tidyverse = as.character(utils::packageVersion("tidyverse")),
  slider = as.character(utils::packageVersion("slider")),
  minpack_lm = as.character(utils::packageVersion("minpack.lm")),
  ggrepel = as.character(utils::packageVersion("ggrepel"))
), file.path(table_dir, "analysis_runtime_environment.csv"))

message("PROJECT_ID:   ", PROJECT_ID)
message("ANALYSIS_ID:  ", ANALYSIS_ID)
message("Input mode:   raw-summary")
message("Output root:  ", out_root)

# ---------------------------
# 1) Sample map
# ---------------------------

design <- readr::read_csv(args$comparison_design, show_col_types = FALSE, col_types = cols(.default = col_character()))
truthy <- function(x) tolower(trimws(as.character(x))) %in% c("1", "true", "yes", "y")
design <- design %>% filter(truthy(include))
selected_ids <- if (is.null(args$comparisons) || !nzchar(args$comparisons)) design$comparison_id else strsplit(args$comparisons, ",", fixed = TRUE)[[1]]
contrast_map <- design %>% filter(comparison_id %in% selected_ids) %>%
  transmute(contrast = comparison_id, background, pool, parent_sample, treated_sample)
if (!nrow(contrast_map)) stop("No selected comparisons were found in the design.")
sample_map <- bind_rows(
  design %>% transmute(sample = parent_sample, condition = parent_condition, background, pool),
  design %>% transmute(sample = treated_sample, condition = treated_condition, background, pool)
) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(pool = factor(pool), background = factor(background), condition = factor(condition, levels = c("parent", "treated")))

sample_paths <- args$sample_table
sample_names <- sub("=.*$", "", sample_paths)
sample_files <- sub("^[^=]*=", "", sample_paths)
if (length(sample_paths) == 0 || any(!nzchar(sample_names)) || any(!file.exists(sample_files))) stop("Every design sample requires an existing --sample-table SAMPLE=FILE input.")
sample_file_map <- stats::setNames(sample_files, sample_names)

# ---------------------------
# 2) Helpers
# ---------------------------

clean_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "_", x))
}

pick_col <- function(df, candidates, label) {
  nms <- names(df)
  clean_nms <- clean_name(nms)

  for (candidate in candidates) {
    hit <- which(clean_nms == clean_name(candidate))
    if (length(hit) > 0) {
      return(nms[[hit[[1]]]])
    }
  }

  message("Available columns:")
  print(nms)

  stop("Could not detect ", label, " column.")
}

read_feature_csv <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  clean_cols <- clean_name(names(df))

  if (!("reads" %in% clean_cols) && !("standard_name" %in% clean_cols)) {
    df2 <- readr::read_csv(path, skip = 1, show_col_types = FALSE)

    if (ncol(df2) > ncol(df)) {
      df <- df2
    }
  }

  df
}

find_feature_table <- function(sample) {
  path <- unname(sample_file_map[[sample]])
  if (is.null(path) || !file.exists(path)) stop("Missing explicit raw SummaryTable input for sample: ", sample)
  path
}

read_sample_feature_table <- function(sample) {
  f <- find_feature_table(sample)
  df <- read_feature_csv(f)

  feature_col <- pick_col(
    df,
    c(
      "standard_name",
      "feature_name",
      "feature",
      "gene",
      "gene_name",
      "locus_tag",
      "orf",
      "qng_id"
    ),
    "feature/gene id"
  )

  reads_col <- pick_col(
    df,
    c(
      "reads",
      "read_count",
      "feature_reads",
      "reads_per_feature",
      "total_reads",
      "sum_reads"
    ),
    "reads"
  )

  df %>%
    mutate(
      sample = sample,
      feature_id = as.character(.data[[feature_col]]),
      reads = suppressWarnings(as.numeric(.data[[reads_col]]))
    ) %>%
    filter(!is.na(feature_id), feature_id != "", !is.na(reads))
}

collapse_first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

# ---------------------------
# 3) Read all feature tables
# ---------------------------

feature_tbl <- purrr::map_dfr(sample_map$sample, read_sample_feature_table) %>%
  left_join(sample_map, by = "sample")

readr::write_csv(
  feature_tbl,
  file.path(table_dir, "feature_table.long.raw_reads.csv")
)

# Keep general annotation columns from feature tables.
# This avoids OSR-specific joins but preserves useful identifiers/descriptions if present.
annotation_cols <- setdiff(
  names(feature_tbl),
  c("sample", "condition", "background", "pool", "reads")
)

feature_annot <- feature_tbl %>%
  select(any_of(annotation_cols)) %>%
  group_by(feature_id) %>%
  summarise(
    across(everything(), collapse_first_nonempty),
    .groups = "drop"
  )

# ---------------------------
# 4) Build gene x sample matrix and CPM
# ---------------------------

reads_long <- feature_tbl %>%
  group_by(feature_id, sample, condition, background, pool) %>%
  summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop")

library_sizes <- reads_long %>%
  group_by(sample, condition, background, pool) %>%
  summarise(
    total_feature_reads = sum(reads, na.rm = TRUE),
    features_with_reads = sum(reads > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(background, pool, condition)

readr::write_csv(
  library_sizes,
  file.path(table_dir, "library_sizes.feature_reads.csv")
)

reads_cpm_long <- reads_long %>%
  left_join(library_sizes %>% select(sample, total_feature_reads), by = "sample") %>%
  mutate(
    cpm = if_else(total_feature_reads > 0, reads / total_feature_reads * 1e6, 0)
  )

readr::write_csv(
  reads_cpm_long,
  file.path(table_dir, "feature_reads_cpm.long.csv")
)

reads_wide <- reads_long %>%
  select(feature_id, sample, reads) %>%
  pivot_wider(names_from = sample, values_from = reads, values_fill = 0)

cpm_wide <- reads_cpm_long %>%
  select(feature_id, sample, cpm) %>%
  pivot_wider(names_from = sample, values_from = cpm, values_fill = 0)

# ---------------------------
# 5) Parent-parent background noise model
# ---------------------------

make_control_running <- function(cpm_wide, ctrl1, ctrl2, w = 40, minA = 0.1) {
  c1 <- cpm_wide[[ctrl1]]
  c2 <- cpm_wide[[ctrl2]]

  ctrl_tab <- tibble(
    feature_id = cpm_wide$feature_id,
    A_ctrl = (c1 + c2) / 2,
    M_ctrl = log2((c1 + PSEUDOCOUNT) / (c2 + PSEUDOCOUNT))
  ) %>%
    arrange(desc(A_ctrl)) %>%
    mutate(
      run_A = slider::slide_dbl(
        A_ctrl,
        mean,
        .before = floor(w / 2),
        .after = floor(w / 2),
        .complete = TRUE
      ),
      run_SD = slider::slide_dbl(
        M_ctrl,
        sd,
        .before = floor(w / 2),
        .after = floor(w / 2),
        .complete = TRUE
      )
    )

  fit_dat <- ctrl_tab %>%
    filter(!is.na(run_A), !is.na(run_SD), run_A >= minA, run_SD > 0)

  if (nrow(fit_dat) < 20) {
    stop("Too few points to fit control-control noise model for ", ctrl1, " vs ", ctrl2)
  }

  lm0 <- lm(log(run_SD) ~ log(run_A), data = fit_dat)

  m3_start <- as.numeric(coef(lm0)[2])
  m3_start <- max(min(m3_start, 0.3), -2.5)

  m1_start <- as.numeric(quantile(fit_dat$run_SD, 0.02, na.rm = TRUE))
  x_ref <- as.numeric(median(fit_dat$run_A, na.rm = TRUE))
  y_ref <- as.numeric(median(fit_dat$run_SD, na.rm = TRUE))
  m2_start <- max((y_ref - m1_start) / (x_ref ^ m3_start), 1e-6)

  start <- list(m1 = m1_start, m2 = m2_start, m3 = m3_start)

  fit <- minpack.lm::nlsLM(
    run_SD ~ m1 + m2 * (run_A ^ m3),
    data = fit_dat,
    start = start,
    lower = c(m1 = 0, m2 = 0, m3 = -5),
    upper = c(m1 = Inf, m2 = Inf, m3 = 0.5),
    control = minpack.lm::nls.lm.control(maxiter = 500)
  )

  list(
    ctrl_tab = ctrl_tab,
    fit_dat = fit_dat,
    fit = fit,
    start = start
  )
}

sd_from_fit <- function(fit_obj, A) {
  cf <- coef(fit_obj)
  m1 <- cf[["m1"]]
  m2 <- cf[["m2"]]
  m3 <- cf[["m3"]]

  pmax(m1 + m2 * (A ^ m3), 1e-8)
}

parent_controls <- sample_map %>% filter(condition == "parent") %>%
  group_by(background) %>% summarise(samples = list(unique(sample)), .groups = "drop")
if (any(lengths(parent_controls$samples) < 2)) stop("Each background requires at least two parent samples for the established parent-parent noise model.")
control_models <- setNames(vector("list", nrow(parent_controls)), as.character(parent_controls$background))
for (i in seq_len(nrow(parent_controls))) {
  controls <- parent_controls$samples[[i]][1:2]
  model <- make_control_running(cpm_wide, controls[[1]], controls[[2]])
  z_control <- model$ctrl_tab$M_ctrl / sd_from_fit(model$fit, model$ctrl_tab$A_ctrl)
  control_models[[as.character(parent_controls$background[[i]])]] <- list(
    model = model, threshold = as.numeric(quantile(abs(z_control), Z_QUANTILE, na.rm = TRUE)),
    z = z_control, controls = controls)
}

z_thresholds <- purrr::imap_dfr(control_models, function(x, background) tibble(
  background = background, z_quantile = Z_QUANTILE, z_threshold = x$threshold,
  parent_comparison = paste(x$controls, collapse = "_vs_")
))

readr::write_csv(
  z_thresholds,
  file.path(table_dir, "background_specific_z_thresholds.csv")
)

fit_coef <- purrr::imap_dfr(control_models, function(x, background) tibble(
  background = background, parameter = names(coef(x$model$fit)), value = as.numeric(coef(x$model$fit))))

readr::write_csv(
  fit_coef,
  file.path(table_dir, "background_noise_model_coefficients.csv")
)

# ---------------------------
# 6) Treated-vs-parent CPM log2FC and z scores
# ---------------------------

compute_contrast <- function(cpm_wide, contrast, background, pool, parent_sample, treated_sample) {
  parent_cpm <- cpm_wide[[parent_sample]]
  treated_cpm <- cpm_wide[[treated_sample]]

  A_exp <- (treated_cpm + parent_cpm) / 2
  M_exp <- log2((treated_cpm + PSEUDOCOUNT) / (parent_cpm + PSEUDOCOUNT))

  control <- control_models[[as.character(background)]]
  if (is.null(control)) stop("No parent-parent noise model for background: ", background)
  fit_obj <- control$model$fit
  z_thr <- control$threshold

  sd_local <- sd_from_fit(fit_obj, A_exp)
  z <- M_exp / sd_local

  tibble(
    feature_id = cpm_wide$feature_id,
    contrast = contrast,
    background = background,
    pool = pool,
    parent_sample = parent_sample,
    treated_sample = treated_sample,
    parent_cpm = parent_cpm,
    treated_cpm = treated_cpm,
    A_exp = A_exp,
    log2FC = M_exp,
    sd_local = sd_local,
    z = z,
    Zthr = z_thr,
    call = case_when(
      z >= z_thr ~ "enriched",
      z <= -z_thr ~ "depleted",
      TRUE ~ "none"
    )
  )
}

contrast_results <- purrr::pmap_dfr(
  contrast_map,
  function(contrast, background, pool, parent_sample, treated_sample) {
    compute_contrast(
      cpm_wide = cpm_wide,
      contrast = contrast,
      background = background,
      pool = pool,
      parent_sample = parent_sample,
      treated_sample = treated_sample
    )
  }
) %>%
  mutate(
    pool = factor(pool),
    background = factor(background),
    call = factor(call, levels = c("depleted", "none", "enriched"))
  ) %>%
  left_join(feature_annot, by = "feature_id")

readr::write_csv(
  contrast_results,
  file.path(table_dir, "treated_vs_parent.by_pool.log2fc_z.csv")
)

summary_results <- contrast_results %>%
  group_by(feature_id) %>%
  summarise(
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    median_log2FC = median(log2FC, na.rm = TRUE),
    sd_log2FC = sd(log2FC, na.rm = TRUE),
    mean_z = mean(z, na.rm = TRUE),
    max_abs_z = max(abs(z), na.rm = TRUE),
    n_pairs = sum(!is.na(log2FC)),
    n_enriched_pairs = sum(call == "enriched", na.rm = TRUE),
    n_depleted_pairs = sum(call == "depleted", na.rm = TRUE),
    pools_enriched = paste(sort(unique(as.character(pool[call == "enriched"]))), collapse = ";"),
    pools_depleted = paste(sort(unique(as.character(pool[call == "depleted"]))), collapse = ";"),
    final_call = case_when(
      n_enriched_pairs >= 2 & n_depleted_pairs == 0 ~ "consistently_enriched",
      n_depleted_pairs >= 2 & n_enriched_pairs == 0 ~ "consistently_depleted",
      n_enriched_pairs >= 1 & n_depleted_pairs >= 1 ~ "mixed",
      n_enriched_pairs == 1 ~ "single_pool_enriched",
      n_depleted_pairs == 1 ~ "single_pool_depleted",
      TRUE ~ "none"
    ),
    .groups = "drop"
  ) %>%
  left_join(feature_annot, by = "feature_id") %>%
  arrange(mean_log2FC)

readr::write_csv(
  summary_results,
  file.path(table_dir, "treated_vs_parent.summary_by_feature.csv")
)

top_depleted <- summary_results %>%
  arrange(mean_log2FC) %>%
  slice_head(n = 100)

top_enriched <- summary_results %>%
  arrange(desc(mean_log2FC) ) %>%
  slice_head(n = 100)

readr::write_csv(
  top_depleted,
  file.path(table_dir, "top100_depleted_in_treated.csv")
)

readr::write_csv(
  top_enriched,
  file.path(table_dir, "top100_enriched_in_treated.csv")
)

call_counts <- contrast_results %>%
  count(background, pool, contrast, call) %>%
  arrange(background, pool, call)

readr::write_csv(
  call_counts,
  file.path(table_dir, "call_counts.by_pool.csv")
)

# ---------------------------
# 7) Plot helpers
# ---------------------------

label_col <- if ("standard_name" %in% names(contrast_results)) {
  "standard_name"
} else {
  "feature_id"
}

contrast_results <- contrast_results %>%
  mutate(
    plot_label = coalesce(
      na_if(str_trim(as.character(.data[[label_col]])), ""),
      feature_id
    )
  )

sig_labels <- contrast_results %>%
  filter(call %in% c("enriched", "depleted")) %>%
  group_by(background, pool, call) %>%
  slice_max(order_by = abs(z), n = TOP_LABEL_N, with_ties = FALSE) %>%
  ungroup()

# ---------------------------
# 8) Plots
# ---------------------------

p_lib <- library_sizes %>%
  ggplot(aes(x = sample, y = total_feature_reads, fill = condition)) +
  geom_col() +
  facet_wrap(~ background, scales = "free_x") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Feature-level library sizes",
    subtitle = paste0(PROJECT_ID, "; raw summary / CPM"),
    x = "Sample",
    y = "Total feature reads"
  )

ggsave(
  file.path(plot_dir, "library_sizes_feature_reads.png"),
  p_lib,
  width = 9,
  height = 5,
  dpi = 300
)

ctrl_z_df <- purrr::imap_dfr(control_models, function(x, background) tibble(background = background, z = x$z))

p_ctrl_z <- ctrl_z_df %>%
  ggplot(aes(x = z)) +
  geom_histogram(bins = 100) +
  facet_wrap(~ background, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Parent-parent control z-score distributions",
    subtitle = paste0("Z threshold quantile = ", Z_QUANTILE),
    x = "Control-control z",
    y = "Number of genes/features"
  )

ggsave(
  file.path(plot_dir, "control_control_z_histogram.png"),
  p_ctrl_z,
  width = 8,
  height = 4.5,
  dpi = 300
)

combined_ma <- contrast_results %>%
  group_by(feature_id) %>%
  summarise(
    A_combined = mean(A_exp, na.rm = TRUE),
    M_combined = mean(log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(summary_results %>% select(feature_id, final_call), by = "feature_id") %>%
  mutate(
    combined_call = case_when(
      final_call %in% c("consistently_depleted", "single_pool_depleted") ~ "depleted",
      final_call %in% c("consistently_enriched", "single_pool_enriched") ~ "enriched",
      final_call == "mixed" ~ "mixed",
      TRUE ~ "none"
    ),
    combined_call = factor(combined_call, levels = c("depleted", "none", "enriched", "mixed"))
  )

p_ma_all <- combined_ma %>%
  ggplot(aes(x = A_combined + 1, y = M_combined)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(color = combined_call), alpha = 0.45, size = 1.1) +
  scale_color_manual(values = c(depleted = "blue", none = "grey70", enriched = "red", mixed = "black")) +
  scale_x_log10() +
  theme_bw(base_size = 12) +
  theme(
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "Combined treated versus parent MA plot",
    subtitle = paste0(PROJECT_ID, "; all selected pools combined; raw summary / CPM"),
    x = "A = mean average CPM across treated-parent comparisons + 1",
    y = "M = mean log2FC treated / parent",
    color = "Call"
  )

ggsave(
  file.path(plot_dir, "MA_treated_vs_parent_combined.png"),
  p_ma_all,
  width = 10,
  height = 7,
  dpi = 300
)

base_ma <- function(df) {
  ggplot(df, aes(x = A_exp + 1, y = log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(aes(color = call), alpha = 0.35, size = 1.2) +
    scale_color_manual(values = c(depleted = "blue", none = "grey70", enriched = "red")) +
    scale_x_log10() +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
      axis.title = element_text(face = "bold", size = 13),
      strip.text = element_text(face = "bold", size = 11),
      axis.text.y = element_text(color = "black", face = "bold"),
      legend.title = element_text(color = "black", face = "bold"),
      legend.text = element_text(color = "black")
    ) +
    labs(
      x = "A = average CPM (treated + parent) / 2 + 1",
      y = "M = log2FC treated / parent",
      color = "Call"
    )
}

for (bg in unique(as.character(contrast_results$background))) {
  plot <- base_ma(contrast_results %>% filter(as.character(background) == bg)) +
    facet_wrap(~ pool, nrow = 1) +
    ggrepel::geom_text_repel(data = sig_labels %>% filter(as.character(background) == bg),
      aes(label = plot_label), size = 3.5, max.overlaps = Inf, box.padding = 0.25,
      point.padding = 0.2, min.segment.length = 0, seed = 1) +
    ggtitle(paste(bg, "treated versus parent MA plots"))
  ggsave(file.path(plot_dir, paste0("MA_", clean_name(bg), "_labeled_top_hits.png")),
         plot, width = 10, height = 5.5, dpi = 300)
}

p_rank <- summary_results %>%
  arrange(mean_log2FC) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = mean_log2FC)) +
  geom_point(alpha = 0.45, size = 0.8) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(
    title = "Ranked treated-versus-parent log2FC",
    subtitle = paste0(PROJECT_ID, "; raw summary / CPM"),
    x = "Genes/features ranked by mean log2FC",
    y = "Mean log2FC"
  )

ggsave(
  file.path(plot_dir, "ranked_mean_log2FC.png"),
  p_rank,
  width = 9,
  height = 5,
  dpi = 300
)

p_hist <- summary_results %>%
  ggplot(aes(x = mean_log2FC)) +
  geom_histogram(bins = 80) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(
    title = "Distribution of mean treated-versus-parent log2FC",
    subtitle = paste0(PROJECT_ID, "; raw summary / CPM"),
    x = "Mean log2FC",
    y = "Number of genes/features"
  )

ggsave(
  file.path(plot_dir, "mean_log2FC_distribution.png"),
  p_hist,
  width = 7,
  height = 5,
  dpi = 300
)

top_heatmap_features <- summary_results %>%
  arrange(mean_log2FC) %>%
  slice_head(n = 25) %>%
  bind_rows(summary_results %>% arrange(desc(mean_log2FC)) %>% slice_head(n = 25)) %>%
  distinct(feature_id)

heatmap_df <- contrast_results %>%
  semi_join(top_heatmap_features, by = "feature_id") %>%
  mutate(
    plot_label = factor(plot_label, levels = rev(unique(top_heatmap_features$feature_id)))
  )

# Use feature_id order if label mapping is not unique.
heatmap_df <- contrast_results %>%
  semi_join(top_heatmap_features, by = "feature_id") %>%
  mutate(feature_id = factor(feature_id, levels = rev(top_heatmap_features$feature_id)))

p_heat <- heatmap_df %>%
  ggplot(aes(x = contrast, y = feature_id, fill = log2FC)) +
  geom_tile() +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    title = "Top depleted and enriched features",
    subtitle = "Log2FC in treated versus parent; raw summary / CPM",
    x = "Matched pool contrast",
    y = "Gene/feature",
    fill = "log2FC"
  )

ggsave(
  file.path(plot_dir, "top_features_log2FC_heatmap.png"),
  p_heat,
  width = 10,
  height = 10,
  dpi = 300
)

message("Done.")
message("Tables: ", table_dir)
message("Plots:  ", plot_dir)
