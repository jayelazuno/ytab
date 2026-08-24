#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)][1])
repo_root <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)
app_file <- file.path(repo_root, "app", "shiny", "app.R")
ui_file <- file.path(repo_root, "app", "shiny", "R", "ui_gene_domain_explorer.R")
server_file <- file.path(repo_root, "app", "shiny", "R", "gene_domain_explorer_server.R")
state_file <- file.path(repo_root, "app", "shiny", "R", "gene_domain_explorer_state.R")
qc_plot_file <- file.path(repo_root, "app", "shiny", "R", "qc_plot_utils.R")
table_helpers_file <- file.path(repo_root, "app", "shiny", "R", "table_display_helpers.R")

fail <- function(message) {
  cat("FAIL:", message, "\n", file = stderr())
  quit(status = 1)
}

tryCatch(invisible(parse(file = app_file)), error = function(e) fail(conditionMessage(e)))
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(qc_plot_file)
source(state_file)

all_text <- paste(vapply(c(app_file, ui_file, server_file, state_file, qc_plot_file, table_helpers_file), function(path) {
  if (!file.exists(path)) fail(paste("missing file", path))
  paste(readLines(path, warn = FALSE), collapse = "\n")
}, ""), collapse = "\n")

required <- c(
  "Gene & Domain Insertion Explorer",
  "gene_domain_query",
  "gene_domain_search",
  "gene_domain_candidates",
  "gene_domain_track_source",
  "gene_domain_track_preset",
  "gene_domain_preset_message",
  "gene_domain_samples",
  "gene_domain_flank_bp",
  "gene_domain_show_domains",
  "gene_domain_show_direction",
  "gene_domain_label_mode",
  "gene_domain_show_site_counts",
  "gene_domain_generate",
  "gene_domain_download_png",
  "gene_domain_download_table"
)

missing <- required[!vapply(required, function(pattern)
  grepl(pattern, all_text, fixed = TRUE), logical(1))]
if (length(missing)) fail(paste("missing UI controls:", paste(missing, collapse = ", ")))

project <- list(species = "glabrata", repo_root = repo_root)
if (!identical(qc_plot_chromosome_display("CP048230.1", project), "Chr A"))
  fail("CP048230.1 does not map to Chr A")
if (!identical(qc_plot_chromosome_display("ChrM_C_glabrata_CBS138", project), "Chr M"))
  fail("ChrM_C_glabrata_CBS138 does not map to Chr M")
display_map_calls <- gregexpr("qc_plot_chromosome_display\\(data\\$chromosome", all_text, perl = TRUE)[[1]]
if (identical(display_map_calls[[1]], -1L) || length(display_map_calls) < 2L)
  fail("Gene Explorer candidate/insertion tables do not both display-map chromosome values")
if (!grepl("qc_plot_chromosome_display\\(row\\$chromosome", all_text))
  fail("Gene Explorer selected-gene summary does not display-map chromosome values")
if (!grepl("\"Chromosome\"", all_text, fixed = TRUE))
  fail("Gene Explorer does not use a readable Chromosome display column")
if (!grepl("write.csv\\(data, file, row.names = FALSE\\)", all_text))
  fail("Gene Explorer insertion table download is not written from display-mapped data")
if (grepl("Internal feature ID", all_text, fixed = TRUE))
  fail("Gene Explorer visible results still expose Internal feature ID")
if (!grepl("ytab_gene_details_datatable\\(data", all_text))
  fail("Gene Explorer search results do not use the shared clickable gene-details table")
if (!grepl("ytab_glabrata_gene_detail_columns\\(data\\)", all_text))
  fail("Gene Explorer search results do not carry hidden gene-detail columns")
if (!grepl("selection = \"single\"", all_text, fixed = TRUE))
  fail("Gene Explorer search results do not preserve single-row selection")

sample_sheet <- data.frame(
  sample = c(
    "parent-pool1-treated",
    "parent-pool1-mock",
    "parent-pool2-treated",
    "parent-pool2-mock",
    "unknown-track"
  ),
  condition = c("treated", "mock", "treated", "mock", ""),
  treatment = c("H2O2", "mock", "1_5mM-Zn", "mock", ""),
  pool = c("1", "1", "2", "2", ""),
  stringsAsFactors = FALSE
)
tracks <- data.frame(
  sample = sample_sheet$sample,
  track_name = sample_sheet$sample,
  stringsAsFactors = FALSE
)
ordered <- gene_domain_order_tracks(tracks, sample_sheet)
choices <- gene_domain_preset_choices(ordered)
if (!all(c("All tracks", "Parents only", "Treated only", "Matched pairs", "Custom") %in% names(choices)))
  fail("track preset choices are incomplete")
parents <- gene_domain_preset_track_rows("parents", ordered)
treated <- gene_domain_preset_track_rows("treated", ordered)
matched <- gene_domain_preset_track_rows("matched_pairs", ordered)
custom <- gene_domain_preset_track_rows("custom", ordered)
if (nrow(parents) != 2L || any(parents$role != "parent")) fail("Parents only preset did not select parent/control tracks")
if (nrow(treated) != 2L || any(treated$role != "treated")) fail("Treated only preset did not select treated tracks")
if (nrow(matched) != 4L || !identical(as.character(matched$role), c("parent", "treated", "parent", "treated")))
  fail("Matched pairs preset did not return parent/treated pairs in pool order")
if (nrow(custom) != nrow(ordered)) fail("Custom preset did not keep all tracks available for manual selection")
if (gene_domain_infer_track_role("parent-pool1-treated", sample_sheet) != "treated")
  fail("treated role did not override parent token in a treated sample name")
for (pattern in c("custom_track_selection", "clear_generated_figure", "No tracks match this preset for the current project."))
  if (!grepl(pattern, all_text, fixed = TRUE)) fail(paste("missing preset synchronization behavior:", pattern))

forbidden_imports <- paste0(
  c("codex/", "from ", "import ", "from ", "import ", "from ", "import ", "from ", "import "),
  c("DomainFigures", "DomainFigures", "DomainFigures", "Organisms", "Organisms",
    "SortedCollection", "SortedCollection", "RangeSet", "RangeSet")
)
hits <- forbidden_imports[vapply(forbidden_imports, function(pattern)
  grepl(pattern, all_text, fixed = TRUE), logical(1))]
if (length(hits)) fail(paste("legacy import reference detected:", paste(hits, collapse = ", ")))

hard_coded_genes <- c("GWK60_A00033", "EPA24", "yH298", "yH299")
new_tab_text <- paste(vapply(c(ui_file, server_file, state_file), function(path)
  paste(readLines(path, warn = FALSE), collapse = "\n"), ""), collapse = "\n")
hits <- hard_coded_genes[vapply(hard_coded_genes, function(pattern)
  grepl(pattern, new_tab_text, fixed = TRUE), logical(1))]
if (length(hits)) fail(paste("hard-coded H2O2-specific identifier detected:", paste(hits, collapse = ", ")))

cat("PASS\n")
