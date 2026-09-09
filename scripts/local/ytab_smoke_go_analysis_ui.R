#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]]),
                      winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/", mustWork = TRUE)

fail <- function(message) stop(message, call. = FALSE)
must_exist <- function(path) {
  full <- file.path(root, path)
  if (!file.exists(full)) fail(paste("Missing file:", path))
  full
}
must_contain <- function(path, pattern, label = pattern) {
  text <- paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
  if (!grepl(pattern, text, fixed = TRUE)) fail(paste("Missing", label, "in", path))
}
must_not_contain <- function(path, pattern, label = pattern) {
  text <- paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
  if (grepl(pattern, text, fixed = TRUE)) fail(paste("Unexpected", label, "in", path))
}

library(shiny)
library(bslib)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left

for (helper in c(
  "app/shiny/R/ui_helpers.R",
  "app/shiny/R/ui_components.R",
  "app/shiny/R/plot_customization_helpers.R",
  "app/shiny/R/plot_display_helpers.R",
  "app/shiny/R/table_display_helpers.R",
  "app/shiny/R/glabrata_annotation_lookup.R",
  "app/shiny/R/go_analysis.R",
  "app/shiny/R/navigation.R"
)) {
  full <- must_exist(helper)
  invisible(parse(file = full))
  source(full, local = TRUE)
}

invisible(parse(file = must_exist("app/shiny/app.R")))
invisible(parse(file = must_exist("scripts/local/ytab_topgo_rrvgo_worker.R")))

for (fn in c(
  "go_analysis_ui", "go_analysis_server", "go_analysis_plot_top_terms",
  "go_analysis_plot_semantic_heatmap", "go_analysis_gene_set_table",
  "go_analysis_discover_results", "go_analysis_go_display_table",
  "go_analysis_semantic_display_table", "go_analysis_display_datatable"
)) {
  if (!exists(fn, mode = "function")) fail(paste("Missing GO helper:", fn))
}

ui <- paste(as.character(go_analysis_ui()), collapse = "\n")
for (needle in c("Gene set", "GO enrichment", "Semantic reduction",
                 "Current Fitness selected top hits", "Pasted genes",
                 "Upload gene list CSV/TXT", "Comparison view", "Hit direction",
                 "Number of genes", "Run GO analysis",
                 "Download selected GO genes CSV", "Download topGO results CSV",
                 "Download semantic reduction CSV", "All ontologies")) {
  if (!grepl(needle, ui, fixed = TRUE)) fail(paste("GO UI missing:", needle))
}
if (grepl("nav_panel(\"Downloads\"", ui, fixed = TRUE))
  fail("GO downloads should be embedded in result tabs, not a separate Downloads subtab")

state <- navigation_state_defaults("go_analysis")
stopifnot(state$active_go_analysis_subtab == "gene_set")
state <- record_navigation_selection(state, "go_analysis_tabs", "semantic_reduction")
stopifnot(state$active_go_analysis_subtab == "semantic_reduction")
stopifnot(navigation_nested_id_for_top("go_analysis") == "go_analysis_tabs")

must_contain("app/shiny/app.R", "go_analysis.R", "GO source file")
must_contain("app/shiny/app.R", "go_analysis_server(input, output, session", "GO server call")
must_contain("app/shiny/R/ui_workspace.R", "GO Analysis", "GO workspace tab")
must_contain("app/shiny/www/ytab_nav_menus.js", "go_analysis", "GO navigation menu")
must_contain("app/shiny/R/job_manager.R", "\"go_analysis\"", "GO progress job support")
must_contain("app/shiny/R/go_analysis.R", "Concordance filter", "GO selected-hit concordance label")
must_contain("app/shiny/R/go_analysis.R", "go_fitness_min_support", "GO selected-hit concordance control")
must_contain("app/shiny/R/go_analysis.R", "fitness_ma_top_hits_table", "GO uses Fitness selected-hit logic")
must_contain("app/shiny/R/go_analysis.R", "go_analysis_selected_hit_display_table", "GO selected-hit metrics table")
must_contain("app/shiny/R/go_analysis.R", "go_analysis_go_display_table", "GO enrichment app-facing display table")
must_contain("app/shiny/R/go_analysis.R", "go_analysis_semantic_display_table", "GO semantic app-facing display table")
must_contain("app/shiny/R/go_analysis.R", "Classic Fisher P", "readable topGO p-value column")
must_contain("app/shiny/R/go_analysis.R", "weight01 Fisher P", "readable topGO weight01 p-value column")
must_contain("app/shiny/R/go_analysis.R", "Representative theme", "readable rrvgo representative column")
must_contain("app/shiny/R/go_analysis.R", "Best enrichment score", "readable rrvgo score column")
must_contain("app/shiny/www/ytab_release_ui.css", "ytab-go-result-table", "GO table wrapping CSS")
must_contain("app/shiny/R/go_analysis.R", "go_results_tick <- reactiveVal(0L)", "GO result refresh tick")
must_contain("app/shiny/R/go_analysis.R", "go_results_tick()\n    active()", "GO results depend on explicit refresh tick")
must_contain("app/shiny/R/go_analysis.R", "latest_go_result <- reactive({\n    go_analysis_latest_result(go_results())", "GO latest result is render-stable")
must_contain("app/shiny/R/go_analysis.R", "observe({\n    if (!identical(jobs$stage(), \"go_analysis\")) return(NULL)", "GO job polling separated from rendered results")
must_contain("app/shiny/R/go_analysis.R", "jobs$poll_job()", "GO background job polling")
must_contain("app/shiny/R/go_analysis.R", "go_results_tick(go_results_tick() + 1L)", "GO results refresh after job completion")
must_not_contain("app/shiny/R/go_analysis.R", "latest_go_result <- reactive({\n    if (identical(jobs$stage(), \"go_analysis\") && jobs$job_is_running())", "GO result reactive polling loop")
must_not_contain("app/shiny/R/go_analysis.R", "go_results <- reactive({\n    invalidateLater(1500, session)", "unconditional GO output-directory refresh")
must_contain("app/shiny/R/go_analysis.R", "ggplot2::geom_col(width = 0.74", "dynamic ggplot GO enrichment bars")
must_contain("app/shiny/R/go_analysis.R", "linewidth = 0.95", "GO enrichment bar borders")
must_contain("app/shiny/R/go_analysis.R", "ggplot2::geom_tile(color = \"grey45\"", "dynamic semantic-reduction heatmap tiles")
must_contain("app/shiny/R/go_analysis.R", "plot.margin = ggplot2::margin", "GO plot anti-clipping margins")
must_contain("app/shiny/R/job_progress.R", "progress_is_stale(p,60)&&!job_manager$job_is_running()", "compact job banner stale guard")
must_not_contain("app/shiny/R/go_analysis.R", "nav_panel(\"Downloads\"", "separate GO Downloads tab")
must_not_contain("app/shiny/R/go_analysis.R", "ytab_static_image", "static GO image helper")
must_not_contain("app/shiny/R/go_analysis.R", "imageOutput", "static GO image output")
must_not_contain("app/shiny/R/go_analysis.R", "Bars are rendered from the current topGO result table", "implementation note in GO enrichment UI")

sample_hits <- data.frame(
  feature_id = c("GWK60_A00001", "GWK60_A00002"),
  cagl_display_id = c("CAGL0A00011g", "CAGL0A00022g"),
  gene_display_name = c("ABC1", "DEF2"),
  selected_rank = 1:2,
  hit_direction = c("depleted", "enriched"),
  stringsAsFactors = FALSE
)
gene_set <- go_analysis_gene_set_table(sample_hits, source = "fitness")
stopifnot(nrow(gene_set) == 2L, "Rank" %in% names(gene_set), "Gene name" %in% names(gene_set))

sample_topgo <- data.frame(
  GO.ID = c("GO:1", "GO:2"),
  Term = c("very long biological process term for display wrapping", "second term"),
  Annotated = c(10, 20),
  Significant = c(3, 1),
  Expected = c(0.2, 0.8),
  classic_numeric = c(0.0012, 0.2),
  weight01_numeric = c(0.0012, 0.2),
  ontology_code = c("BP", "BP"),
  ontology_label = c("Biological Process", "Biological Process"),
  target_genes_with_GO = c(15, 15),
  universe_genes_with_GO = c(1000, 1000),
  check.names = FALSE
)
go_display <- go_analysis_go_display_table(sample_topgo, ontology = "BP", p_cutoff = 0.05, top_n = 1)
stopifnot(
  nrow(go_display) == 2L,
  identical(go_display$`Plot term`[[1L]], "yes"),
  all(c("GO ID", "Classic Fisher P", "weight01 Fisher P", "-log10(weight01 P)") %in% names(go_display)),
  !"weight01_numeric" %in% names(go_display),
  !"ontology_code" %in% names(go_display)
)

sample_semantic <- data.frame(
  go = c("GO:1", "GO:2"),
  cluster = c(1, 1),
  parent = c("GO:1", "GO:1"),
  score = c(2.3456, 1.2345),
  size = c(2, 2),
  term = c("member process", "second member process"),
  parentTerm = c("representative process", "representative process"),
  termUniqueness = c(0.8, 0.7),
  termUniquenessWithinCluster = c(0.5, 0.6),
  termDispensability = c(0, 0.4),
  ontology_code = c("BP", "BP"),
  ontology_label = c("Biological Process", "Biological Process"),
  semantic_similarity_threshold = c(0.7, 0.7),
  semantic_similarity_method = c("Rel", "Rel"),
  check.names = FALSE
)
semantic_display <- go_analysis_semantic_display_table(sample_semantic, ontology = "BP")
stopifnot(
  nrow(semantic_display) == 2L,
  all(c("Representative GO ID", "Representative theme", "Member GO term",
        "Best enrichment score", "Semantic threshold") %in% names(semantic_display)),
  !"parentTerm" %in% names(semantic_display),
  !"termUniquenessWithinCluster" %in% names(semantic_display)
)

cat("PASS\n")
