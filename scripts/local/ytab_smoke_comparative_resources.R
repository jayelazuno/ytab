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
          # S. pombe was enabled in the manifest when its projects were added; this assertion was left behind
          # asserting the opposite, so the check has been failing since. It is a working species here now.
          identical(species$enabled[species$species == "pombe"], TRUE),
          is.data.frame(lookup),
          is.data.frame(orthology),
          is.data.frame(available))
# The check is that no path outside this checkout is baked in, not that none is absolute: these are resolved at
# run time against the repo root, so on any machine with the repo under /Users every one of them contains it.
stopifnot(all(startsWith(available$path[nzchar(available$path)], root)))
paths <- paste(c(sub(root, "", available$path, fixed = TRUE), capture.output(str(species))), collapse = "\n")
stopifnot(!grepl("/Users/", paths, fixed = TRUE),
          !grepl(paste0("codex/", paste0("RNA", "cross")), paths, fixed = TRUE))
if (nrow(lookup) && nrow(orthology)) {
  mapped <- map_gene_across_species(lookup$gene_id[[1]], orthology, lookup)
  stopifnot(is.data.frame(mapped))
}
cat("PASS\n")
