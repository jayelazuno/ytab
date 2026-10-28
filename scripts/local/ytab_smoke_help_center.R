#!/usr/bin/env Rscript
# Help center: every indexed methods page, the tutorial and the release notes render; every image they use exists;
# search finds tabs, pages, command options and datasets; every dataset's details render. Usage: Rscript scripts/local/ytab_smoke_help_center.R REPO_ROOT
args <- commandArgs(trailingOnly = TRUE); root <- normalizePath(args[[1L]])
suppressPackageStartupMessages(library(shiny))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
source(file.path(root, "app/shiny/R/datasets_center.R")); source(file.path(root, "app/shiny/R/help_center.R"))
index <- utils::read.delim(file.path(root, "docs/methods/index.tsv"), stringsAsFactors = FALSE)
stopifnot(all(file.exists(file.path(root, index$file))))
pages <- help_pages(root)
stopifnot(nrow(pages) == nrow(index) + 2L, all(c("tutorial", "RELEASE_NOTES") %in% pages$id))
for (i in seq_len(nrow(pages))) {
  html <- as.character(help_render(root, pages$file[[i]]))
  stopifnot(nchar(html) > 200)
  src <- regmatches(html, gregexpr('src="ytab-docs/[^"]+"', html))[[1]]
  files <- file.path(root, "docs", sub('^src="ytab-docs/', "", sub('"$', "", src)))
  if (length(files) && !all(file.exists(files))) stop("Missing image in ", pages$file[[i]], ": ", paste(files[!file.exists(files)], collapse = ", "))
}
hits <- help_search(root, "colonies", pages); stopifnot("screen_coverage" %in% hits$id[hits$kind == "tab"])
hits <- help_search(root, "beta-geometric", pages); stopifnot(any(hits$kind == "page"))
stopifnot(!nrow(help_search(root, "", pages)))
hits <- help_search(root, "--project-config", pages); stopifnot("cli_reference" %in% hits$id)
sets <- datasets_catalog(root); stopifnot(nrow(sets$catalog) > 0, nrow(sets$libraries) > 0)
for (i in seq_len(nrow(sets$catalog))) stopifnot(nchar(as.character(datasets_detail_ui(as.list(sets$catalog[i, ]), sets$libraries))) > 200)
run <- sets$libraries$sample[[1]]; hits <- help_search(root, run, pages, datasets = sets)
stopifnot(any(hits$kind == "dataset" & hits$id == sets$libraries$study[[1]]))
studies <- unique(sets$catalog$study); stopifnot(nrow(datasets_table_data(sets$catalog, sets$libraries)) == length(studies))
cat("PASS\n")
