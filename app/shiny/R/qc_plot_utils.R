qc_plot_empty <- function(message) {
  old <- qc_plot_par()
  on.exit(par(old), add = TRUE)
  plot.new()
  text(0.5, 0.5, message, cex = 1.05, col = "#617285")
}

qc_plot_numeric <- function(x) suppressWarnings(as.numeric(gsub("%$", "", as.character(x))))

qc_plot_column <- function(data, candidates, n = nrow(data)) {
  hit <- intersect(candidates, names(data))
  if (length(hit)) data[[hit[[1]]]] else rep(NA, n)
}

qc_plot_sample_order <- function(samples) {
  canonical <- c(
    "yH298-parent-pool1", "yH298-parent-pool2",
    "yH299-parent-pool3", "yH299-parent-pool4",
    "yH298-H2O2-treated-facs-pool1", "yH298-H2O2-treated-facs-pool2",
    "yH299-H2O2-treated-facs-pool3", "yH299-H2O2-treated-facs-pool4"
  )
  c(intersect(canonical, samples), sort(setdiff(samples, canonical)))
}

qc_plot_condition <- function(sample) {
  ifelse(grepl("parent", sample, ignore.case = TRUE), "parent",
         ifelse(grepl("treated|H2O2|facs", sample, ignore.case = TRUE),
                "H2O2-treated-facs", "other"))
}

qc_plot_pool <- function(sample) {
  value <- sub(".*pool([0-9]+).*", "\\1", sample)
  ifelse(value == sample, "", value)
}

qc_plot_card <- function(title, output_id, height = "300px") {
  tags$article(class = "ytab-plot-card", tags$h4(title), plotOutput(output_id, height = height))
}

qc_plot_par <- function(mar = c(5, 5, 3, 1), ...) {
  par(
    mar = mar,
    font.axis = 2,
    font.lab = 2,
    font.main = 2,
    cex.axis = 1.2,
    cex.lab = 1.25,
    cex.main = 1.25,
    lwd = 1.35,
    ...
  )
}

qc_plot_fill <- "grey70"
qc_plot_fill_light <- "grey85"
qc_plot_border <- "black"
qc_plot_lwd <- 1.4
