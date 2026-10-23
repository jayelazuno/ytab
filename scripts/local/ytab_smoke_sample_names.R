#!/usr/bin/env Rscript
all <- commandArgs(FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", all, value=TRUE)[[1]]))
root <- normalizePath(file.path(dirname(self), "../.."))
library(shiny)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
for (file in c("project_discovery.R", "sample_selector.R", "qc_plot_utils.R"))
  source(file.path(root, "app/shiny/R", file))

# Arbitrary IDs: no species, study or naming convention required.
samples <- data.frame(sample=c("instrument-001", "run-XYZ", "sample_C"),
  display_name=c("Untreated culture 1", "Drug culture 1", ""),
  include=c("True", "True", "False"), condition=c("control", "treated", "parent"),
  fastq_1=c("reads/a.fastq.gz", "reads/b.fastq.gz", "reads/c.fastq.gz"))
labels <- ytab_sample_labels(samples$sample, samples)
stopifnot(identical(labels, c("Untreated culture 1", "Drug culture 1", "sample_C")),
  identical(ytab_sample_labels(c("unknown", "run-XYZ"), samples), c("unknown", "Drug culture 1")),
  identical(ytab_sample_labels(samples$sample, samples[setdiff(names(samples), "display_name")]), samples$sample))
status <- data.frame(sample=rev(samples$sample), mapping=c("failed", "success", "pending"))
table <- compact_sample_table(samples, status, "run-XYZ", include_selected=TRUE)
stopifnot(identical(table$Sample, labels), identical(table$Selected, c("No", "Yes", "No")),
  identical(table$Mapping, c("pending", "success", "failed")))
server <- function(input, output, session) {
  session$userData$chosen <- sample_selector_server("test", reactive(samples), status_data=reactive(status))
}
shiny::testServer(server, {
  session$flushReact()
  stopifnot(setequal(session$userData$chosen(), samples$sample[1:2]))
  session$setInputs(`test-samples`="run-XYZ")
  session$flushReact()
  stopifnot(identical(session$userData$chosen(), "run-XYZ"))
})

directory <- tempfile("ytab-sample-names-")
dir.create(directory)
sheet <- file.path(directory, "sample_sheet.csv")
write.csv(samples, sheet, row.names=FALSE)
before <- read.csv(sheet, colClasses="character", check.names=FALSE)
checksum <- unname(tools::md5sum(sheet))
expect_error <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error=identity), "error"))
expect_error(ytab_save_sample_labels(sheet, samples$sample, rep("duplicate", 3), checksum))
stopifnot(identical(unname(tools::md5sum(sheet)), checksum))
ytab_save_sample_labels(sheet, samples$sample, c("New control", "New treatment", ""), checksum)
after <- read.csv(sheet, colClasses="character", check.names=FALSE)
stopifnot(identical(before[setdiff(names(before), "display_name")], after[setdiff(names(after), "display_name")]),
  identical(ytab_sample_labels(samples$sample, after), c("New control", "New treatment", "sample_C")))
expect_error(ytab_save_sample_labels(sheet, samples$sample, labels, checksum))
expect_error(ytab_save_sample_labels(sheet, rev(samples$sample), labels, unname(tools::md5sum(sheet))))
for (file in c("app/shiny/app.R", "app/shiny/R/ui_workspace.R", "app/shiny/R/fitness_design_state.R",
               "app/shiny/R/gene_domain_explorer_state.R", "app/shiny/R/qc_library_diagnostics_plots.R"))
  invisible(parse(file.path(root, file)))
cat("PASS: generic aliases, selector IDs, status joins, metadata preservation, duplicate/stale-edit protection and app syntax\n")
