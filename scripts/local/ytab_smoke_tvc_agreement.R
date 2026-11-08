#!/usr/bin/env Rscript
# Treated vs Control > Compare methods: `ytab tvc-agreement` on Gale 2020's z-score, resampling and Gumbel results writes
# the agreement tables, and every plot of the view builds from them. Usage: Rscript scripts/local/ytab_smoke_tvc_agreement.R REPO_ROOT
args <- commandArgs(trailingOnly = TRUE); root <- normalizePath(if (length(args)) args[[1L]] else ".")
if (grepl("[.]yaml$", root)) root <- normalizePath(file.path(dirname(root), "../../../.."))
suppressPackageStartupMessages({library(shiny); library(bslib)})
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
for (name in c("ui_helpers.R", "ui_components.R", "plot_customization_helpers.R", "plot_display_helpers.R", "table_display_helpers.R",
               "analysis_runs.R", "project_discovery.R", "ui_fitness.R", "treated_vs_control.R")) source(file.path(root, "app/shiny/R", name))
g <- file.path(root, "output/projects/Gale_2020"); out <- tempfile("agreement_")
python <- file.path(dirname(dirname(R.home())), "bin", "python"); if (!file.exists(python)) python <- Sys.which("python3")
status <- system2(python, c(file.path(root, "scripts/local/ytab_tvc_agreement.py"), "--zscore", file.path(g, "treated_vs_parent/treated_vs_parent_BGC_22Oct2026_1048"),
                            "--resampling", file.path(g, "treated_vs_control/resampling_BGC"), "--gumbel", file.path(g, "treated_vs_control/gumbel_20261104_165910"), "--output", out),
                  stdout = FALSE, stderr = FALSE)
stopifnot(status == 0, all(file.exists(file.path(out, c("genes.tsv", "pairs.tsv", "sets.tsv", "manifest.json")))))
pairs <- ytab_read_analysis_tsv(file.path(out, "pairs.tsv")); sets <- ytab_read_analysis_tsv(file.path(out, "sets.tsv")); genes <- ytab_read_analysis_tsv(file.path(out, "genes.tsv"))
stopifnot(nrow(pairs) == 3L, all(is.na(pairs$enriched_a_recovered_by_b[pairs$method_b == "gumbel"])), !any(grepl("jaccard", names(pairs))),
          all(pairs$depleted_a_recovered_by_b > 0 & pairs$depleted_a_recovered_by_b <= 1),
          all(abs(pairs$depleted_a_recovered_by_b - pairs$depleted_n_both / pairs$depleted_n_a) < 1e-9),
          pairs$spearman_effect[pairs$method_a == "zscore" & pairs$method_b == "resampling"] > .9, !any(sets$direction == "enriched" & grepl("gumbel", sets$methods)))
s <- ytab_plot_customization_values(list(), "tvcag")
for (d in c("depleted", "enriched")) invisible(ytab_ggplotly(tvc_agreement_matrix(pairs, d, s)))
for (w in c("effect", "score")) invisible(ytab_ggplotly(tvc_agreement_scatter(genes, "zscore", "gumbel", w, s)))
invisible(ytab_ggplotly(tvc_agreement_upset(sets[sets$direction == "depleted", ], c("zscore", "resampling", "gumbel"), s)))
v <- tvc_agreement_venn(sets[sets$direction == "depleted", ], c("zscore", "resampling", "gumbel"), s); invisible(ytab_ggplotly(v))
# the seven Venn regions add up to every gene called by at least one of the three methods
stopifnot(sum(as.numeric(gsub(",", "", v$layers[[2]]$data$label))) == sum(sets$genes[sets$direction == "depleted"]))
invisible(ytab_ggplotly(tvc_agreement_venn(sets[sets$direction == "enriched", ], c("zscore", "resampling"), s)))
html <- htmltools::renderTags(tvc_results_tab_ui())$html; stopifnot(grepl("Compare methods", html, fixed = TRUE), grepl("tvc_ag_inputs", html, fixed = TRUE))
cat("PASS\n")
