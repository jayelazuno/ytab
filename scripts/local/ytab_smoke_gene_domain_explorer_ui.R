#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)][1])
repo_root <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)
app_file <- file.path(repo_root, "app", "shiny", "app.R")
ui_file <- file.path(repo_root, "app", "shiny", "R", "ui_gene_domain_explorer.R")
server_file <- file.path(repo_root, "app", "shiny", "R", "gene_domain_explorer_server.R")
state_file <- file.path(repo_root, "app", "shiny", "R", "gene_domain_explorer_state.R")

fail <- function(message) {
  cat("FAIL:", message, "\n", file = stderr())
  quit(status = 1)
}

tryCatch(invisible(parse(file = app_file)), error = function(e) fail(conditionMessage(e)))

all_text <- paste(vapply(c(app_file, ui_file, server_file, state_file), function(path) {
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
