#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1L]] else default
}
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(
  sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1L]]),
  winslash = "/"
)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
rnacross_dir <- get_arg("--rnacross-dir", paste0("codex/", "RNAcross"))
out_dir <- get_arg("--out-dir", "resources/comparative")
script <- file.path(root, "scripts/local/ytab_import_rnacross_resources.R")
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  shQuote(c(script, "--rnacross-dir", rnacross_dir, "--out-dir", out_dir)),
  stdout = TRUE, stderr = TRUE
)
code <- attr(status, "status") %||% 0L
if (!identical(as.integer(code), 0L)) {
  writeLines(status)
  stop("RNAcross import failed.", call. = FALSE)
}
expected <- file.path(root, out_dir, c(
  "orthology/gene_lookup.csv",
  "orthology/orthogroups_long.csv",
  "orthology/qng_gwk_cagl_threeway_map_complete.csv",
  "manifest/comparative_species_manifest.yaml",
  "manifest/resource_import_manifest.json"
))
stopifnot(all(file.exists(expected)))
lookup <- read.csv(expected[[1]], stringsAsFactors = FALSE, check.names = FALSE)
orthology <- read.csv(expected[[2]], stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(lookup) > 0L,
          all(c("gene_id", "species", "gene_name", "og_id") %in% names(lookup)),
          all(c("orthogroup_id", "species", "gene_id", "gene_name",
                "source_column") %in% names(orthology)))
cat("PASS\n")
