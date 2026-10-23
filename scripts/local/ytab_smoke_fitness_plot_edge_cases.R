#!/usr/bin/env Rscript
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]]))
root <- normalizePath(file.path(dirname(self), "../.."))
`%||%` <- function(x,y) if(is.null(x)||!length(x)) y else x
source(file.path(root, "app/shiny/R/fitness_plot_utils.R"))
directory <- tempfile("ytab-fitness-edge-")
dir.create(file.path(directory,"tables"), recursive=TRUE)
result <- list(output_dir=directory)
save_table <- function(x, filename) write.csv(x, file.path(directory,"tables",filename), row.names=FALSE)
# Contract columns (ytab.fitness.contract): the plot abundance comes from the pipeline's summary.
save_table(data.frame(feature_id=c("gene-a","gene-b"),mean_log2FC=c(-2,-3),call_log2FC=c(-2,-3),final_call="none",mean_abundance=c(10,20)),
  "treated_vs_parent.summary_by_feature.csv")
pool <- data.frame(feature_id=c("gene-a","gene-b"),contrast="control-treatment-1",pool="1",
  log2FC=c(-2,-3),z=c(-2,-3),call="none",parent_cpm=c(10,20),treated_cpm=c(2,3),A_exp=c(6,11.5))
save_table(pool,"treated_vs_parent.by_pool.log2fc_z.csv")
# A result without the pipeline's top-hit evidence (made before it existed): no hits, and a warning saying why.
data <- fitness_ma_rank_data(result)
stopifnot(nrow(data)==2, !any(data$support_pass), all(nzchar(data$ranking_warning)))
# With the evidence (ytab.fitness.zscore.hit_evidence), the app only reads it and applies the viewer's choices.
evidence <- data.frame(feature_id=c("gene-a","gene-b"),standard_name=c("A","B"),common_name=c("A","B"),call="none",log2FC=c(-2,-3),z=c(-2,-3),
  mean_abundance=c(6,11.5),n_pairs_valid=1L,n_pairs_depleted_support=1L,n_pairs_enriched_support=0L,n_pairs_depleted_called=0L,
  n_pairs_enriched_called=0L,pools_depleted_support="pool1",pools_enriched_support="",passes_abundance=TRUE,
  z_strength_depleted=c(2,3),rank_depleted=c(2,1),z_strength_enriched=0,rank_enriched=NA,z_strength_both=c(2,3),rank_both=c(2,1))
save_table(evidence,"top_hit_evidence.csv"); save_table(cbind(contrast="control-treatment-1",evidence),"top_hit_evidence.by_pair.csv")
for (direction in c("depleted","enriched","both","none")) {
  data <- fitness_ma_rank_data(result,direction=direction)
  stopifnot(nrow(data)==2, all(data$enriched_support_n==0), all(data$depleted_stored_call_n==0), all(data$enriched_stored_call_n==0))
  hits <- fitness_ma_highlight_rows(data,direction,10L,0L)
  stopifnot(identical(as.character(data$feature_id[hits]), if (direction %in% c("depleted","both")) c("gene-b","gene-a") else character()))
}
stopifnot(!length(fitness_ma_highlight_rows(fitness_ma_rank_data(result,direction="depleted",min_support=2L),"depleted",10L,2L)))
auto <- fitness_ma_rank_data(result,"individual")
explicit <- fitness_ma_rank_data(result,"individual","control-treatment-1")
stopifnot(identical(auto,explicit))
cat("PASS: missing evidence, pipeline order per direction, support filter, automatic pair selection\n")
