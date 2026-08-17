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

qc_plot_text_scale <- function() {
  switch(getOption("ytab.plot.text_size", "medium"),
         small = 0.9, large = 1.18, 1)
}

qc_plot_label_las <- function(default = 2L) {
  angle <- as.character(getOption("ytab.plot.label_angle", "90"))
  if (identical(angle, "0")) return(1L)
  if (identical(angle, "90")) return(2L)
  default
}

qc_plot_label_angle <- function(default = 0L) {
  angle <- suppressWarnings(as.integer(getOption("ytab.plot.label_angle", default)))
  if (is.na(angle)) default else angle
}

qc_plot_label_mode <- function() {
  as.character(getOption("ytab.plot.label_mode", "full"))
}

qc_plot_display_labels <- function(labels) {
  labels <- as.character(labels)
  mode <- qc_plot_label_mode()
  if (identical(mode, "hide")) return(rep("", length(labels)))
  if (!identical(mode, "compact")) return(labels)
  compact <- labels
  compact <- sub("^yH[0-9]+-", "", compact)
  compact <- gsub("H2O2-treated-facs", "H2O2", compact, fixed = TRUE)
  compact <- gsub("parent", "parent", compact, fixed = TRUE)
  compact <- gsub("-pool", " p", compact, fixed = TRUE)
  compact
}

qc_plot_grid_enabled <- function() {
  identical(as.character(getOption("ytab.plot.grid", "show")), "show")
}

qc_plot_bar_horizontal <- function(labels = character()) {
  mode <- as.character(getOption("ytab.plot.bar_orientation", "auto"))
  if (identical(mode, "horizontal")) return(TRUE)
  if (identical(mode, "vertical")) return(FALSE)
  if (!length(labels)) return(FALSE)
  max(nchar(as.character(labels)), na.rm = TRUE) >= 18L
}

qc_plot_show_value_labels <- function() {
  isTRUE(getOption("ytab.plot.show_value_labels", TRUE))
}

qc_plot_par <- function(mar = c(5, 5, 3, 1), ...) {
  scale <- qc_plot_text_scale()
  par(
    mar = mar,
    font.axis = 2,
    font.lab = 2,
    font.main = 2,
    cex.axis = 1.2 * scale,
    cex.lab = 1.25 * scale,
    cex.main = 1.25 * scale,
    lwd = 1.35,
    ...
  )
}

qc_plot_fill <- "grey70"
qc_plot_fill_light <- "grey85"
qc_plot_border <- "black"
qc_plot_lwd <- 1.4
