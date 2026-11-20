# The gene/domain insertion figure, drawn by the app from the payload the explorer already returns.
#
# The explorer's CLI also writes a PNG of this figure. That image is not shown: it cannot be restyled, resized,
# hovered or exported to PDF, and it bypasses the house plot style. Every value it draws is already in the JSON
# the app receives -- `insertions` (one row per insertion, with position, track and region_class), `tracks`
# (the lane order and per-track counts), `gene` and `region` (the gene model and the window) and `domains` --
# so the app builds the plot itself and the PNG stays on disk as the pipeline's own record.

# Lane order follows the track summaries, so the figure's rows match the track list the user chose.
gene_domain_track_levels <- function(manifest) {
  tracks <- manifest$track_summaries
  if (is.null(tracks) || !length(tracks)) return(character())
  vapply(tracks, function(t) as.character(t$track_name %||% t$sample %||% ""), "")
}

# "<track name>  n=<insertions in gene>", the label the figure has always carried.
gene_domain_track_labels <- function(manifest) {
  tracks <- manifest$track_summaries
  if (is.null(tracks) || !length(tracks)) return(character())
  stats::setNames(
    vapply(tracks, function(t) sprintf("%s  n=%s", as.character(t$track_name %||% t$sample %||% ""),
                                       as.character(t$insertions_in_gene %||% 0L)), ""),
    gene_domain_track_levels(manifest))
}

gene_domain_insertion_frame <- function(manifest) {
  path <- as.character(manifest$table_path %||% "")
  if (!length(path) || !nzchar(path) || !file.exists(path)) return(data.frame())
  d <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d) || !"position" %in% names(d)) return(data.frame())
  d$position <- suppressWarnings(as.numeric(d$position))
  d$track <- as.character(if ("track_name" %in% names(d)) d$track_name else d$sample)
  blank <- is.na(d$track) | !nzchar(d$track)
  if (any(blank) && "sample" %in% names(d)) d$track[blank] <- as.character(d$sample)[blank]
  d$region_class <- if ("region_class" %in% names(d)) as.character(d$region_class) else NA_character_
  d$reads <- suppressWarnings(as.numeric(if ("read_count" %in% names(d)) d$read_count else NA))
  d[is.finite(d$position) & !is.na(d$track), c("position", "track", "region_class", "reads"), drop = FALSE]
}

# The insertion lanes and the gene model, as one figure.
#   payload   the explorer's parsed JSON manifest
#   s         plot customization values
gene_domain_figure_ggplot <- function(manifest, s, base) {
  d <- gene_domain_insertion_frame(manifest)
  levels_ <- gene_domain_track_levels(manifest)
  if (!nrow(d) || !length(levels_)) return(NULL)
  labels <- gene_domain_track_labels(manifest)
  d <- d[d$track %in% levels_, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  # Lanes read top to bottom in the order the track list gives, which reversing the factor achieves on a
  # discrete y axis.
  d$track <- factor(d$track, levels = rev(levels_))
  d$lane <- as.integer(d$track)

  x0 <- suppressWarnings(as.numeric(manifest$displayed_region_start %||% manifest$start %||% min(d$position)))
  x1 <- suppressWarnings(as.numeric(manifest$displayed_region_end %||% manifest$end %||% max(d$position)))
  gene_start <- suppressWarnings(as.numeric(manifest$gene_start %||% NA))
  gene_end <- suppressWarnings(as.numeric(manifest$gene_end %||% NA))

  # In the gene body the insertion is evidence about the gene; in the flank it is context. The figure has always
  # drawn that difference as ink, and it is the one thing a reader must not have to work out.
  d$where <- ifelse(d$region_class %in% "gene_body", "in gene", "flanking")
  d$hover <- sprintf("%s<br>%s bp<br>%s%s", d$track, format(d$position, big.mark = ","), d$where,
                     ifelse(is.finite(d$reads), sprintf("<br>%s reads", format(d$reads, big.mark = ",")), ""))

  p <- ggplot2::ggplot(d, ggplot2::aes(x = position, y = lane, colour = where, text = hover))
  if (is.finite(gene_start) && is.finite(gene_end)) {
    # The gene's extent behind the lanes, so a tick's position relative to the gene is readable without
    # looking down at the model.
    p <- p + ggplot2::annotate("rect", xmin = gene_start, xmax = gene_end,
                               ymin = .4, ymax = length(levels_) + .6, fill = "#2f6fb5", alpha = .07)
  }
  p <- p +
    ggplot2::geom_hline(yintercept = seq_along(levels_), colour = "grey90", linewidth = .3) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lane - .34, ymax = lane + .34), linewidth = .55) +
    ggplot2::scale_colour_manual(values = c(`in gene` = "#1a1a1a", flanking = "#b9c0c7"), drop = FALSE) +
    ggplot2::scale_y_continuous(breaks = seq_along(levels_), labels = unname(labels[rev(levels_)]),
                                expand = ggplot2::expansion(add = .6)) +
    ggplot2::coord_cartesian(xlim = c(x0, x1)) +
    ggplot2::labs(x = sprintf("%s coordinate (bp)",
                              as.character(manifest$chromosome_display %||% manifest$chromosome %||% "")),
                  y = "Insertion tracks", colour = NULL) + base
  p
}

# The gene model and its domains, drawn under the lanes as a second panel.
gene_domain_model_ggplot <- function(manifest, s, base, gene = list()) {
  gene_start <- suppressWarnings(as.numeric(manifest$gene_start %||% NA))
  gene_end <- suppressWarnings(as.numeric(manifest$gene_end %||% NA))
  if (!is.finite(gene_start) || !is.finite(gene_end)) return(NULL)
  x0 <- suppressWarnings(as.numeric(manifest$displayed_region_start %||% manifest$start %||% gene_start))
  x1 <- suppressWarnings(as.numeric(manifest$displayed_region_end %||% manifest$end %||% gene_end))
  strand <- as.character(gene$strand %||% "")
  title <- trimws(paste(as.character(manifest$resolved_gene_name %||% gene$common_name %||% ""),
                        sprintf("(%s)", as.character(manifest$resolved_gene_id %||% gene$standard_name %||% ""))))

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = gene_start, xmax = gene_end, ymin = .72, ymax = 1.28,
                      fill = "#2f6fb5", colour = "black", linewidth = .3)
  # The arrow states the reading direction, which a bare rectangle cannot.
  if (nzchar(strand)) {
    span <- (gene_end - gene_start) * .32
    from <- if (identical(strand, "-")) gene_start + span else gene_end - span
    to <- if (identical(strand, "-")) gene_start + span * .05 else gene_end - span * .05
    p <- p + ggplot2::annotate("segment", x = from, xend = to, y = 1, yend = 1,
                               arrow = grid::arrow(length = grid::unit(.12, "cm"), type = "closed"),
                               linewidth = .6, colour = "black")
  }
  domains <- manifest$domains
  if (!is.null(domains) && length(domains)) {
    dd <- do.call(rbind, lapply(domains, function(x) {
      v <- suppressWarnings(as.numeric(unlist(x)[1:2]))
      if (any(!is.finite(v))) return(NULL)
      data.frame(start = min(v), end = max(v))
    }))
    if (!is.null(dd) && nrow(dd))
      p <- p + ggplot2::annotate("rect", xmin = dd$start, xmax = dd$end, ymin = .58, ymax = 1.42,
                                 fill = "#c8553d", alpha = .75, colour = NA)
  }
  p + ggplot2::annotate("text", x = (gene_start + gene_end) / 2, y = 1.62, label = title,
                        fontface = "bold", size = (s$text_size %||% 12) / 3.4) +
    ggplot2::coord_cartesian(xlim = c(x0, x1), ylim = c(.4, 1.85)) +
    ggplot2::labs(x = NULL, y = NULL) + base +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
}
