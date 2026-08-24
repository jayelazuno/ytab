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

qc_plot_repo_root <- function(project = NULL) {
  root <- as.character((project %||% list())$repo_root %||% "")
  if (nzchar(root)) return(root)
  project_root <- as.character((project %||% list())$project_root %||% "")
  if (nzchar(project_root)) {
    return(normalizePath(file.path(project_root, "../../.."),
                         winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

qc_plot_glabrata_chromosome_map <- function(project = NULL) {
  fallback <- data.frame(
    contig = paste0("CP048", 230:242, ".1"),
    display = paste("Chr", LETTERS[1:13]),
    order = seq_len(13L),
    stringsAsFactors = FALSE
  )
  species <- tolower(as.character((project %||% list())$species %||% ""))
  if (nzchar(species) && !identical(species, "glabrata")) return(data.frame())
  repo_root <- qc_plot_repo_root(project)
  path <- file.path(repo_root, "resources", "species", "glabrata", "reference_genome", "chr_to_cp.tsv")
  if (!file.exists(path)) return(fallback)
  map <- tryCatch(read.delim(path, header = FALSE, stringsAsFactors = FALSE,
                             col.names = c("chromosome_name", "contig")),
                  error = function(e) data.frame())
  if (!nrow(map) || !all(c("chromosome_name", "contig") %in% names(map))) return(fallback)
  map$letter <- sub("^Chr([A-Z]).*$", "\\1", as.character(map$chromosome_name))
  map$display <- ifelse(grepl("^[A-Z]$", map$letter), paste("Chr", map$letter),
                        as.character(map$chromosome_name))
  map$order <- match(map$letter, LETTERS)
  map <- map[is.finite(map$order), c("contig", "display", "order"), drop = FALSE]
  if (!nrow(map)) return(fallback)
  map[order(map$order), , drop = FALSE]
}

qc_plot_chromosome_display <- function(chrom, project = NULL) {
  chrom <- as.character(chrom)
  map <- qc_plot_glabrata_chromosome_map(project)
  if (!nrow(map)) return(chrom)
  aliases <- data.frame(
    key = c(map$contig, map$display, gsub(" ", "", map$display),
            paste0(gsub(" ", "", map$display), "_C_glabrata_CBS138")),
    display = rep(map$display, times = 4L),
    order = rep(map$order, times = 4L),
    stringsAsFactors = FALSE
  )
  idx <- match(chrom, aliases$key)
  out <- ifelse(is.na(idx), chrom, aliases$display[idx])
  out
}

qc_plot_chromosome_order <- function(chrom, project = NULL) {
  chrom <- as.character(chrom)
  map <- qc_plot_glabrata_chromosome_map(project)
  if (!nrow(map)) return(chrom)
  display <- qc_plot_chromosome_display(chrom, project)
  known <- map$display[map$display %in% unique(display)]
  c(known, sort(setdiff(unique(display), known)))
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
  text_scale <- switch(getOption("ytab.plot.text_size", "medium"),
                       small = 0.95, large = 1.55, 1.2)
  text_scale * qc_plot_style_scale("text")
}

qc_plot_style <- function() {
  style <- as.character(getOption("ytab.plot.style", "clean"))
  if (style %in% c("clean", "compact", "presentation")) style else "clean"
}

qc_plot_style_scale <- function(component = c("text", "margin", "line", "point", "key")) {
  component <- match.arg(component)
  style <- qc_plot_style()
  values <- switch(
    component,
    text = c(clean = 1.00, compact = 0.82, presentation = 1.24),
    margin = c(clean = 1.00, compact = 0.72, presentation = 1.16),
    line = c(clean = 1.00, compact = 0.78, presentation = 1.35),
    point = c(clean = 1.00, compact = 0.84, presentation = 1.22),
    key = c(clean = 1.00, compact = 0.82, presentation = 1.18)
  )
  unname(values[[style]] %||% values[["clean"]])
}

qc_plot_style_margins <- function(mar) {
  mar <- as.numeric(mar)
  scale <- qc_plot_style_scale("margin")
  if (identical(qc_plot_style(), "clean")) return(mar)
  minimum <- c(3.2, 3.2, 2.2, 0.8)
  pmax(minimum, mar * scale)
}

qc_plot_text_sizes <- function() {
  scale <- qc_plot_text_scale()
  list(
    axis = 1.05 * scale,
    sample = 1.05 * scale,
    lab = 1.18 * scale,
    main = 1.28 * scale,
    key = 0.98 * scale,
    value = 0.95 * scale
  )
}

qc_plot_label_las <- function(default = 2L) {
  angle <- as.character(getOption("ytab.plot.label_angle", "30"))
  if (identical(angle, "0")) return(1L)
  if (identical(angle, "90")) return(2L)
  default
}

qc_plot_label_angle <- function(default = 30L) {
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

qc_plot_wrapped_labels <- function(labels, width = NULL) {
  labels <- as.character(labels)
  if (is.null(width)) {
    mode <- qc_plot_label_mode()
    width <- if (identical(mode, "compact")) 16L else 24L
  }
  vapply(labels, function(label) {
    if (!nzchar(label)) return("")
    paste(strwrap(label, width = width), collapse = "\n")
  }, character(1))
}

qc_plot_grid_enabled <- function() {
  identical(as.character(getOption("ytab.plot.grid", "show")), "show")
}

qc_plot_bar_horizontal <- function(labels = character()) {
  mode <- as.character(getOption("ytab.plot.bar_orientation", "vertical"))
  if (identical(mode, "horizontal")) return(TRUE)
  if (identical(mode, "vertical")) return(FALSE)
  if (!length(labels)) return(FALSE)
  max(nchar(as.character(labels)), na.rm = TRUE) >= 18L
}

qc_plot_show_value_labels <- function() {
  isTRUE(getOption("ytab.plot.show_value_labels", TRUE))
}

qc_plot_label_cex <- function(multiplier = 1) {
  multiplier * qc_plot_text_sizes()$sample
}

qc_plot_value_cex <- function(multiplier = 1) {
  multiplier * qc_plot_text_sizes()$value
}

qc_plot_key_cex <- function(multiplier = 1) {
  multiplier * qc_plot_style_scale("key") * qc_plot_text_sizes()$key
}

qc_plot_point_cex <- function(multiplier = 1) {
  multiplier * qc_plot_style_scale("point")
}

qc_plot_label_margin_lines <- function(labels, angle = 0L, orientation = c("horizontal", "vertical")) {
  orientation <- match.arg(orientation)
  labels <- as.character(labels)
  visible <- labels[nzchar(labels)]
  if (!length(visible)) return(if (orientation == "horizontal") 7 else 5)
  scale <- qc_plot_text_scale()
  longest <- max(nchar(visible), na.rm = TRUE)
  if (orientation == "horizontal") {
    return(max(22, min(60, longest * 0.82 * scale)))
  }
  if (angle == 90L) return(max(15, min(42, longest * 0.48 * scale)))
  if (angle == 45L) return(max(12, min(34, longest * 0.36 * scale)))
  if (angle == 30L) return(max(10, min(30, longest * 0.30 * scale)))
  max(8, min(22, ceiling(longest / 12) * 2.3 * scale + 4))
}

qc_plot_bar_margins <- function(labels, horizontal = FALSE, top = 3.5, right = 2.4) {
  labels <- qc_plot_display_labels(labels)
  if (isTRUE(horizontal)) {
    c(6.2, qc_plot_label_margin_lines(labels, orientation = "horizontal"), top, right)
  } else {
    angle <- qc_plot_label_angle()
    c(qc_plot_label_margin_lines(labels, angle = angle, orientation = "vertical"),
      5.5, top, right)
  }
}

qc_plot_draw_vertical_labels <- function(at, labels, angle = qc_plot_label_angle(),
                                         cex = qc_plot_label_cex(1)) {
  labels <- qc_plot_display_labels(labels)
  axis(1, at = at, labels = FALSE)
  if (!length(labels) || !any(nzchar(labels))) return(invisible())
  text_labels <- if (identical(as.integer(angle), 0L)) qc_plot_wrapped_labels(labels, width = 18L) else labels
  usr <- par("usr")
  line_height <- strheight("M", cex = cex)
  y <- usr[[3L]] - line_height * if (identical(as.integer(angle), 0L)) 1.2 else 1.8
  adj <- if (identical(as.integer(angle), 90L)) c(1, 0.5) else if (angle > 0) c(1, 1) else c(0.5, 1)
  text(at, y, labels = text_labels, srt = angle, adj = adj, xpd = NA, cex = cex, font = 2)
}

qc_plot_draw_horizontal_labels <- function(at, labels, cex = qc_plot_label_cex(1)) {
  labels <- qc_plot_display_labels(labels)
  axis(2, at = at, labels = FALSE, las = 1)
  if (!length(labels) || !any(nzchar(labels))) return(invisible())
  axis(2, at = at, labels = qc_plot_wrapped_labels(labels, width = 32L),
       las = 1, tick = FALSE, cex.axis = cex, font = 2)
}

qc_plot_add_grid <- function(horizontal = FALSE) {
  if (!qc_plot_grid_enabled()) return(invisible())
  if (isTRUE(horizontal)) {
    abline(v = pretty(par("usr")[1:2]), col = "#e3e9ee", lty = 3, lwd = 1)
  } else {
    abline(h = pretty(par("usr")[3:4]), col = "#e3e9ee", lty = 3, lwd = 1)
  }
}

qc_plot_palette <- function(n) {
  if (n <= 0L) return(character())
  grDevices::gray.colors(n, start = 0.2, end = 0.7)
}

qc_plot_begin_key_layout <- function(key_height = 0.16) {
  layout(matrix(c(1, 2), ncol = 1), heights = c(1, key_height))
}

qc_plot_reset_layout <- function() {
  layout(1)
}

qc_plot_metric_key_row <- function(labels, fill = NULL, col = NULL, border = NULL,
                                   pch = NULL, lty = NULL, lwd = NULL, ncol = NULL,
                                   title = "Metric key") {
  labels <- as.character(labels)
  labels <- labels[nzchar(labels)]
  if (!length(labels)) return(invisible())
  if (is.null(col)) col <- rep("black", length(labels))
  if (is.null(lwd)) lwd <- rep(qc_plot_lwd, length(labels))
  par(mar = c(0.2, 0.5, 0.2, 0.5))
  plot.new()
  ncol <- ncol %||% min(3L, length(labels))
  legend(
    "center",
    legend = labels,
    fill = fill,
    col = col,
    border = border,
    pch = pch,
    lty = lty,
    lwd = lwd,
    ncol = ncol,
    title = title,
    bty = "n",
    text.font = 2,
    cex = qc_plot_key_cex(0.92),
    x.intersp = 0.8,
    y.intersp = 0.95
  )
}

qc_plot_par <- function(mar = c(5, 5, 3, 1), ...) {
  sizes <- qc_plot_text_sizes()
  par(
    mar = qc_plot_style_margins(mar),
    font.axis = 2,
    font.lab = 2,
    font.main = 2,
    cex.axis = sizes$axis,
    cex.lab = sizes$lab,
    cex.main = sizes$main,
    lwd = 1.35 * qc_plot_style_scale("line"),
    ...
  )
}

qc_plot_fill <- "grey70"
qc_plot_fill_light <- "grey85"
qc_plot_border <- "black"
.qc_plot_lwd_base <- 1.4
if (exists("qc_plot_lwd", inherits = FALSE))
  rm(qc_plot_lwd, envir = environment())
makeActiveBinding(
  "qc_plot_lwd",
  function(value) {
    if (missing(value)) return(.qc_plot_lwd_base * qc_plot_style_scale("line"))
    .qc_plot_lwd_base <<- as.numeric(value)
    invisible(.qc_plot_lwd_base)
  },
  environment()
)
