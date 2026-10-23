#!/usr/bin/env Rscript
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]]))
root <- normalizePath(file.path(dirname(self), "../.."))
`%||%` <- function(x,y) if(is.null(x)||!length(x)) y else x
source(file.path(root,"app/shiny/R/qc_plot_utils.R"))
source(file.path(root,"app/shiny/R/plot_customization_helpers.R"))
linear <- qc_plot_count_display(c(0,100000,600000))
stopifnot(identical(linear$values,c(0,100000,600000)),linear$divisor==1000,
  identical(qc_plot_axis_labels(linear$values/linear$divisor),c("0","100","600")))
logged <- qc_plot_count_display(c(0,9,99),"log10")
stopifnot(identical(logged$values,c(0,1,2)),logged$divisor==1,
  grepl("count + 1",logged$suffix,fixed=TRUE))
stopifnot(identical(qc_plot_axis_labels(c(0,200000,400000,600000)),
  c("0","200,000","400,000","600,000")),
  !any(grepl("[eE][+-]",qc_plot_axis_labels(c(1e-8,1e8)))))
old <- options(scipen=0)
ytab_with_plot_display_options(list(),"test", {
  stopifnot(getOption("scipen")==999,
    !grepl("[eE][+-]",format(600000)),
    !grepl("[eE][+-]",format(0.00000001)))
})
ytab_with_plot_display_options(list(test_count_scale="log10"),"test", {
  stopifnot(identical(getOption("ytab.plot.count_scale"),"log10"))
})
stopifnot(getOption("scipen")==0)
options(old)
cat("PASS: comma-separated count ticks, decimal notation and shared plot/export options\n")
