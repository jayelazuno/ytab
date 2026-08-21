#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
self <- normalizePath(sub(file_arg, "", args[grep(file_arg, args)][1]),
                      winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/",
                      mustWork = TRUE)

`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(root, "app/shiny/R/glabrata_annotation_lookup.R"))
source(file.path(root, "app/shiny/R/fitness_plot_utils.R"))

fail <- function(message) {
  cat("FAIL:", message, "\n", file = stderr())
  quit(status = 1)
}

lookup_path <- ytab_glabrata_annotation_lookup_file(root)
if (!file.exists(lookup_path)) fail("new glabrata annotation lookup file is missing")
lookup <- ytab_load_glabrata_annotation_lookup(root)
if (!nrow(lookup)) fail("lookup did not load")
required <- ytab_glabrata_annotation_required_columns()
if (!all(required %in% names(lookup))) fail("required lookup columns missing")
for (key in c("gwk60_key", "cagl_key", "qng_key", "scer_gene_id_key",
              "scer_gene_name_key", "cgla_common_name_key"))
  if (!key %in% names(lookup)) fail(paste("missing normalized key", key))
if ("cgla_gene_name_from_deseq_key" %in% names(lookup))
  fail("DESeq gene-name key should not be created")

hit <- lookup[nzchar(lookup$gwk60_id_clean) & nzchar(lookup$cagl_id), , drop = FALSE][1, ]
input <- data.frame(feature_id = c(hit$gwk60_id_clean[[1]], "GWK60_UNMATCHED"),
                    value = c(1, 2), stringsAsFactors = FALSE)
joined <- ytab_join_glabrata_display(input, root, "feature_id")
if (nrow(joined) != nrow(input)) fail("join changed row count")
if (!identical(joined$feature_id, input$feature_id)) fail("join changed internal feature IDs")
if (joined$cagl_display_id[[1]] != hit$cagl_id[[1]]) fail("CAGL display did not use cagl_id")
if (joined$cagl_display_id[[2]] != "GWK60_UNMATCHED") fail("CAGL display fallback did not use original feature ID")
expected_gene <- if (nzchar(hit$scer_gene_name[[1]])) hit$scer_gene_name[[1]] else hit$cagl_id[[1]]
if (joined$gene_display_name[[1]] != expected_gene) fail("gene display fallback is incorrect")
if (joined$gene_display_name[[2]] != "GWK60_UNMATCHED") fail("unmatched gene display fallback is incorrect")
if (!"cg_to_sc_relationship_display" %in% names(joined)) fail("relationship display column missing")
if (!"gwk60_id_clean" %in% names(joined)) fail("raw lookup columns were not preserved internally")

plot_label_input <- data.frame(feature_id = hit$gwk60_id_clean[[1]],
                               label = hit$gwk60_id_clean[[1]],
                               stringsAsFactors = FALSE)
plot_label_data <- fitness_ma_apply_display_labels(plot_label_input, root)
if (plot_label_data$label[[1]] != expected_gene)
  fail("Fitness plot labels do not use the glabrata display mapping")

fallback <- data.frame(feature_id = c("a", "b", "c"), stringsAsFactors = FALSE)
fallback$cagl_id <- c("CAGL_TEST", "", "")
fallback$scer_gene_name <- c("SCER_TEST", "", "")
fallback <- ytab_join_glabrata_display(fallback[, "feature_id", drop = FALSE], root, "feature_id")
if (!all(c("cagl_display_id", "gene_display_name") %in% names(fallback)))
  fail("fallback display columns missing")

standardized <- ytab_standardize_glabrata_annotation_names(joined)
if ("cgla_gene_name_from_deseq" %in% names(standardized))
  fail("DESeq gene name exposed after standardization")
for (name in c("CAGL ID", "Gene name", "Cg-to-Sc relationship", "Original feature ID"))
  if (!name %in% names(standardized)) fail(paste("standardized user-facing name missing:", name))

if (!file.exists(file.path(root, "resources/comparative/orthology/gene_lookup.csv")))
  fail("gene_lookup.csv was replaced or removed")
if (!file.exists(file.path(root, "resources/comparative/orthology/orthogroups_long.csv")))
  fail("orthogroups_long.csv was replaced or removed")

helper_text <- paste(readLines(file.path(root, "app/shiny/R/glabrata_annotation_lookup.R"),
                               warn = FALSE), collapse = "\n")
if (grepl("qng_gwk_cagl_threeway_map_complete", helper_text, fixed = TRUE))
  fail("new helper requires the old three-way map")
if (grepl("C. glabrata DESeq gene name", helper_text, fixed = TRUE))
  fail("DESeq gene name label exposed by helper")

cat("PASS\n")
