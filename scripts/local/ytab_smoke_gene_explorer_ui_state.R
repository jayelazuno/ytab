#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: ytab_smoke_gene_explorer_ui_state.R <project.yaml>", call. = FALSE)
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1L]]),
                      winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/", mustWork = TRUE)

library(shiny)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left

source(file.path(repo_root, "app/shiny/R/ui_components.R"))
source(file.path(repo_root, "app/shiny/R/plot_customization_helpers.R"))
source(file.path(repo_root, "app/shiny/R/plot_display_helpers.R"))
source(file.path(repo_root, "app/shiny/R/gene_domain_explorer_state.R"))

samples <- data.frame(
  sample = c(
    "yH298-parent-pool1", "yH298-parent-pool2",
    "yH299-parent-pool3", "yH299-parent-pool4",
    "yH298-H2O2-treated-facs-pool1", "yH298-H2O2-treated-facs-pool2",
    "yH299-H2O2-treated-facs-pool3", "yH299-H2O2-treated-facs-pool4"
  ),
  stringsAsFactors = FALSE
)
tracks <- data.frame(sample = samples$sample, track_name = samples$sample,
                     stringsAsFactors = FALSE)
ordered <- gene_domain_order_tracks(tracks, samples)

parents <- gene_domain_preset_track_rows("parents", ordered)
treated <- gene_domain_preset_track_rows("treated", ordered)
all_tracks <- gene_domain_preset_track_rows("all", ordered)
pool1 <- gene_domain_preset_track_rows("pool1_pair", ordered)
matched <- gene_domain_preset_track_rows("matched_pairs", ordered)

stopifnot(
  nrow(parents) == 4L,
  all(parents$role == "parent"),
  nrow(treated) == 4L,
  all(treated$role == "treated"),
  nrow(all_tracks) == 8L,
  identical(as.character(all_tracks$sample[1:4]), as.character(parents$sample)),
  nrow(pool1) == 2L,
  identical(as.character(pool1$role), c("parent", "treated")),
  nrow(matched) == 8L,
  identical(as.character(matched$role), rep(c("parent", "treated"), 4L)),
  all(c("Parents only", "Treated only", "Matched pairs", "Custom") %in%
        names(gene_domain_preset_choices(ordered)))
)

empty_html <- htmltools::renderTags(ytab_empty_state(
  "Search for a gene, choose tracks, then click Generate figure.",
  "The figure preview and downloads appear here after a successful run."
))$html
failure_html <- htmltools::renderTags(ytab_result_card(
  "Gene Explorer could not generate a figure.",
  status = "warning",
  tags$p("Reason: test failure"),
  tags$details(class = "ytab-technical-details", tags$summary("Show technical details"))
))$html

stopifnot(
  grepl("Search for a gene, choose tracks, then click Generate figure.", empty_html, fixed = TRUE),
  grepl("Gene Explorer could not generate a figure.", failure_html, fixed = TRUE),
  grepl("Show technical details", failure_html, fixed = TRUE),
  !grepl("character vector argument expected", empty_html, fixed = TRUE),
  !grepl("character vector argument expected", failure_html, fixed = TRUE)
)

cat("PASS\n")
