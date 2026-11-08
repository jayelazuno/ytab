ytab_plot_output <- function(output_id, height = "medium") {
  plotOutput(output_id, width = "100%", height = ytab_plot_height_px(height))
}

ytab_plot_width_class <- function(choice) {
  choice <- choice %||% "standard"
  if (identical(choice, "wide")) return("ytab-plot-width-wide")
  if (identical(choice, "full")) return("ytab-plot-width-full")
  "ytab-plot-width-standard"
}

ytab_plot_frame <- function(plot_ui, width = "standard", type = "app-rendered") {
  tags$div(
    class = paste("ytab-plot-frame", paste0("ytab-plot-frame-", type),
                  ytab_plot_width_class(width)),
    plot_ui
  )
}

ytab_static_image_ui <- function(src, alt = "", max_height = "760px", show_filename = FALSE) {
  tagList(
    tags$div(
      class = "ytab-static-image",
      tags$img(src = src, alt = alt, loading = "lazy",
               style = paste0("max-height:", max_height, ";"))
    ),
    if (isTRUE(show_filename)) tags$small(class = "text-muted", alt)
  )
}

ytab_static_image_output_card <- function(title, output_id, description = NULL,
                                          height = "large", downloads = NULL,
                                          details = NULL) {
  tags$article(
    class = "ytab-static-image-card ytab-release-card",
    tags$div(
      class = "ytab-plot-card-header",
      tags$div(
        tags$h4(title),
        if (!is.null(description)) tags$p(class = "text-muted", description)
      ),
      if (!is.null(downloads)) tags$div(class = "ytab-actions", downloads)
    ),
    tags$div(
      class = "ytab-static-image ytab-static-image-output",
      imageOutput(output_id, width = "100%", height = ytab_plot_height_px(height))
    ),
    details
  )
}

ytab_generated_file_gallery <- function(files, title_col = "title", url_col = "served_url",
                                        filename_col = "filename") {
  if (!is.data.frame(files) || !nrow(files)) {
    return(ytab_empty_state("No generated files are available."))
  }
  tags$div(
    class = "ytab-generated-file-gallery",
    lapply(seq_len(nrow(files)), function(i) {
      title <- as.character(files[[title_col]][[i]] %||% files[[filename_col]][[i]])
      url <- as.character(files[[url_col]][[i]] %||% "")
      filename <- as.character(files[[filename_col]][[i]] %||% title)
      tags$article(
        class = "ytab-static-image-card ytab-release-card",
        tags$div(class = "ytab-plot-card-header",
                 tags$h4(title)),
        if (nzchar(url)) ytab_static_image_ui(url, filename, "520px"),
        tags$details(tags$summary("Show filename"), tags$code(filename))
      )
    })
  )
}

# Libraries on two principal components (QC library correlation and Normalization; scores and variance from the
# pipelines' pca.tsv / pca_variance.tsv). x: one row per library with the PC columns, `label` and the colour/shape
# columns; variance: pc, variance_explained (for the axis titles). Library names are kept as the plot's "tvc_labels"
# (draggable annotations through tvc_plotly(); ytab_static_labels() draws them for downloads).
# Shared plot helpers used by both the Treated vs Control results and the QC/Normalization PCA plots, so they
# live here rather than in one module's file.
# Labels rows `rows` of `d` (points at columns x, y) on plot `p`: kept as the plot's "tvc_labels" attribute, which
# tvc_plotly() draws as draggable, de-overlapped annotations, or drawn with ggrepel for downloads (interactive = FALSE).
cmp_label_plot <- function(p, d, rows, x, y, s, interactive = FALSE) {
  p <- p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .08)) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.05, .12)))
  l <- d[rows, , drop = FALSE]; if (!nrow(l)) return(p)
  l$tag <- ifelse(!is.na(l$Name) & nzchar(l$Name) & l$Name != l$ORF, l$Name, l$ORF); l$.lx <- as.numeric(l[[x]]); l$.ly <- as.numeric(l[[y]])
  l <- l[is.finite(l$.lx) & is.finite(l$.ly), , drop = FALSE]; if (!nrow(l)) return(p)
  attr(p, "tvc_labels") <- l[, c(".lx", ".ly", "tag")]
  if (isTRUE(interactive)) return(p)
  p + ggrepel::geom_text_repel(data = l, ggplot2::aes(x = .lx, y = .ly, label = tag), inherit.aes = FALSE, colour = "black",
    size = s$text_size / 3.6, min.segment.length = 0, segment.colour = "grey40", box.padding = .5, max.overlaps = Inf, seed = 1)
}

tvc_theme <- function(s) {
  size <- s$text_size + switch(s$style %||% "clean", compact = -2, presentation = 3, 0)
  ggplot2::theme_minimal(base_size = size) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = if (identical(s$grid, "show")) ggplot2::element_line(colour = "grey90") else ggplot2::element_blank())
}


ytab_pca_plot <- function(x, axes, variance, s, colour = "condition", shape = "", labels = TRUE) {
  pct <- function(pc) { v <- suppressWarnings(as.numeric(variance$variance_explained[variance$pc == pc]))
    if (length(v) && is.finite(v[[1L]])) sprintf("%s (%.0f%% of variance)", pc, 100 * v[[1L]]) else pc }
  x$.x <- as.numeric(x[[axes[[1L]]]]); x$.y <- as.numeric(x[[axes[[2L]]]])
  value <- function(col) { v <- if (nzchar(col) && col %in% names(x)) as.character(x[[col]]) else rep("", nrow(x)); ifelse(is.na(v) | !nzchar(v), "(none)", v) }
  x$.colour <- value(colour); x$.shape <- value(shape); use_shape <- nzchar(shape) && length(unique(x$.shape)) <= 6L
  x$.tip <- sprintf("%s\n%s: %s%s\n%s %.2f, %s %.2f", x$label, colour, x$.colour, if (use_shape) sprintf("\n%s: %s", shape, x$.shape) else "",
                    axes[[1L]], x$.x, axes[[2L]], x$.y)
  aes <- if (use_shape) ggplot2::aes(x = .x, y = .y, colour = .colour, shape = .shape, text = .tip) else ggplot2::aes(x = .x, y = .y, colour = .colour, text = .tip)
  p <- ggplot2::ggplot(x, aes) + ggplot2::geom_hline(yintercept = 0, colour = "grey85") + ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_point(size = (s$point_size %||% 2) + 1.5, alpha = max(.6, s$point_opacity %||% .8)) +
    ggplot2::labs(x = pct(axes[[1L]]), y = pct(axes[[2L]]), colour = tools::toTitleCase(colour), shape = tools::toTitleCase(shape)) +
    tvc_theme(s) + ggplot2::theme(legend.position = "bottom")
  d <- data.frame(ORF = x$label, Name = x$label, .x = x$.x, .y = x$.y)
  if (isTRUE(labels)) cmp_label_plot(p, d, seq_len(nrow(d)), ".x", ".y", s, interactive = TRUE) else p
}

# A plot's "tvc_labels" drawn with ggrepel, for downloads of plots shown through tvc_plotly().
ytab_static_labels <- function(p, s) {
  l <- attr(p, "tvc_labels"); if (is.null(l) || !nrow(l)) return(p)
  p + ggrepel::geom_text_repel(data = l, ggplot2::aes(x = .lx, y = .ly, label = tag), inherit.aes = FALSE, colour = "black",
    size = s$text_size / 3.6, min.segment.length = 0, segment.colour = "grey40", box.padding = .5, max.overlaps = Inf, seed = 1)
}
