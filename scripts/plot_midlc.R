#!/usr/bin/env Rscript


# usage 
# Rscript scripts/plot_midlc.R --project_root=.

#Rscript scripts/plot_midlc_qc.R \
 # --project_root=. \
  #--samples=yH298-parent-pool1,yH298-parent-pool2

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
    figure_name    = "midlc_overview.png"
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

standardize_midlc_cols <- function(df) {
  sampled_candidates <- names(df)[tolower(names(df)) %in% c(
    "sampled",
    "reads sampled",
    "reads_sampled",
    "sampled reads",
    "sampled_reads"
  )]

  trial_cols <- grep("^Unique Sites Trial", names(df), value = TRUE)

  # Wide-format MidLC table, e.g.:
  #   Reads Sampled | Unique Sites Trial1 | Unique Sites Trial2 | Unique Sites Trial3
  if (length(sampled_candidates) >= 1 && length(trial_cols) >= 1) {
    sampled_col <- sampled_candidates[1]

    df <- df %>%
      rename(sampled_reads = all_of(sampled_col)) %>%
      pivot_longer(
        cols = all_of(trial_cols),
        names_to = "trial",
        values_to = "unique_sites"
      ) %>%
      mutate(trial = stringr::str_replace(trial, "^Unique Sites\\s*", ""))

    return(df)
  }

  rename_pairs <- c(
    Sampled = "sampled_reads",
    `Reads Sampled` = "sampled_reads",
    depth = "sampled_reads",
    sampled_depth = "sampled_reads",
    reads = "sampled_reads",
    unique = "unique_sites",
    unique_insertions = "unique_sites",
    unique_hits = "unique_sites",
    sites = "unique_sites",
    replicate = "trial"
  )

  for (old_name in names(rename_pairs)) {
    new_name <- rename_pairs[[old_name]]
    if (old_name %in% names(df) && !new_name %in% names(df)) {
      names(df)[names(df) == old_name] <- new_name
    }
  }

  if (!"trial" %in% names(df) && all(c("sampled_reads", "unique_sites") %in% names(df))) {
    df$trial <- "Trial1"
  }

  df
}


read_midlc_files <- function(diagnostics_dir, selected_samples = NULL) {
  all_files <- list.files(
    diagnostics_dir,
    pattern = "\\.midlc\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(all_files) == 0) {
    stop("No .midlc.csv files were found under: ", diagnostics_dir, call. = FALSE)
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
    file_tbl <- file_tbl %>% arrange(sample)
  }

  midlc_tbl <- map_dfr(seq_len(nrow(file_tbl)), function(i) {
    path_i <- file_tbl$path[i]
    sample_i <- file_tbl$sample[i]

    read_csv(path_i, show_col_types = FALSE) %>%
      standardize_midlc_cols() %>%
      mutate(sample = sample_i, .before = 1)
  })

  list(files = file_tbl, data = midlc_tbl)
}

read_midlc_summary <- function(diagnostics_dir) {
  summary_file <- file.path(diagnostics_dir, "library_diagnostics.summary.csv")
  if (!file.exists(summary_file)) return(NULL)

  summ <- read_csv(summary_file, show_col_types = FALSE)
  if (!all(c("sample", "midlc_est") %in% names(summ))) return(NULL)
  summ
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  project_root <- normalizePath(args$project_root, winslash = "/", mustWork = FALSE)
  diagnostics_dir <- file.path(project_root, "output", "smoketests", "library_diagnostics")
  export_img_dir <- file.path(project_root, "output", "exports", args$export_version, "qc", "images")
  dir.create(export_img_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(diagnostics_dir)) {
    stop("Library diagnostics directory was not found: ", diagnostics_dir, call. = FALSE)
  }

  midlc_obj <- read_midlc_files(diagnostics_dir, args$samples)
  midlc_tbl <- midlc_obj$data
  file_tbl <- midlc_obj$files

  required_cols <- c("sampled_reads", "unique_sites")
  missing_cols <- setdiff(required_cols, names(midlc_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "The MidLC files are missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  summary_tbl <- read_midlc_summary(diagnostics_dir)
  if (!is.null(summary_tbl)) {
    summary_tbl <- summary_tbl %>% filter(sample %in% file_tbl$sample)
  }

  p <- ggplot(midlc_tbl, aes(x = sampled_reads, y = unique_sites, group = interaction(sample, trial))) +
    geom_line(alpha = 0.28, linewidth = 0.5, color = "grey60") +
    stat_summary(
      aes(group = sample),
      fun = mean,
      geom = "line",
      linewidth = 1.1,
      color = "black"
    ) +
    facet_wrap(~ sample, scales = "free_y") +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Library complexity (MidLC)",
      subtitle = "Thin lines show subsampling trials; bold line shows the mean curve",
      x = "Sampled reads",
      y = "Unique insertion sites"
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

  if (!is.null(summary_tbl)) {
    p <- p +
      geom_vline(
        data = summary_tbl,
        aes(xintercept = midlc_est),
        linetype = "dashed",
        linewidth = 0.7,
        color = "firebrick",
        inherit.aes = FALSE
      )
  }

  out_file <- file.path(export_img_dir, args$figure_name)
  ggsave(
    filename = out_file,
    plot = p,
    width = max(10, 3.4 * nrow(file_tbl)),
    height = max(6, 3 + ceiling(nrow(file_tbl) / 2) * 2.4),
    dpi = 300,
    bg = "white"
  )

  message("Saved MidLC figure: ", out_file)
}

main()
