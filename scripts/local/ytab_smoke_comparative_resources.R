#!/usr/bin/env Rscript
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(
  sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1L]]),
  winslash = "/"
)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(root, "app/shiny/R/comparative_resources.R"))
species <- load_comparative_species_manifest(root)
lookup <- load_comparative_gene_lookup(root)
orthology <- load_comparative_orthogroups(root)
available <- comparative_resources_available(root)
stopifnot(nrow(species) == 5L,
          all(c("glabrata", "albicans", "cerevisiae", "lactis", "pombe") %in%
                species$species),
          identical(species$enabled[species$species == "pombe"], FALSE),
          is.data.frame(lookup),
          is.data.frame(orthology),
          is.data.frame(available))
paths <- paste(c(available$path, capture.output(str(species))), collapse = "\n")
stopifnot(!grepl("/Users/", paths, fixed = TRUE),
          !grepl(paste0("codex/", paste0("RNA", "cross")), paths, fixed = TRUE))
if (nrow(lookup) && nrow(orthology)) {
  mapped <- map_gene_across_species(lookup$gene_id[[1]], orthology, lookup)
  stopifnot(is.data.frame(mapped))
}
cat("PASS\n")
