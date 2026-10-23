#!/usr/bin/env Rscript
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]]))
root <- normalizePath(file.path(dirname(self), "../.."))
`%||%` <- function(x,y) if(is.null(x)||!length(x)) y else x
for (file in c("project_discovery.R", "qc_plot_utils.R", "qc_summary_library_plots.R"))
  source(file.path(root,"app/shiny/R",file))
directory <- tempfile("ytab-summary-labels-")
for (sample in c("sample-a","sample-b")) {
  dir.create(file.path(directory,"summary",sample),recursive=TRUE)
  writeLines(c("# summary", "standard_name,hits", "gene-1,1", "gene-2,0"),
    file.path(directory,"summary",sample,paste0(sample,".feature_table.RDF_1.csv")))
}
project <- list(project_root=directory,samples=data.frame(sample=c("sample-a","sample-b"),
  condition=c("control","treated")))
labels <- qc_summary_combined_feature_group_labels(project)
stopifnot(identical(labels$treated,"Treated combined"),
  nrow(qc_summary_combined_features_hit_data(project))==2L,
  identical(qc_summary_combined_feature_default_group(project),"both"))
project$samples$treatment_label <- c("","Drug-X")
stopifnot(identical(qc_summary_combined_feature_group_labels(project)$treated,"Drug-X-treated combined"))
project$samples <- project$samples[1,,drop=FALSE]
stopifnot(identical(qc_summary_combined_feature_default_group(project),"control"))
empty <- list(project_root=tempfile(),samples=data.frame())
stopifnot(length(qc_summary_combined_feature_group_choices(empty))==0L,
  identical(qc_summary_combined_feature_default_group(empty),"control"))
cat("PASS: generic and named treatments, single-group and empty projects\n")
