# Normalization module: site-level normalization of a project's insertion counts with any available method
# (`ytab normalize-project`), MidLC common saturation included.

normalization_method_labels <- c(betageom = "Beta-geometric correction (BGC)", ttr = "Trimmed total reads (TTR)",
  totreads = "Total reads", nzmean = "Nonzero mean", aBGC = "Adaptive BGC", quantile = "Zero-preserving quantile",
  zinfnb = "Zero-inflated negative binomial", midlc = "MidLC common saturation", nonorm = "Unnormalized")
# Short names for this tab's plots; the checkboxes use the full names.
normalization_method_short <- c(betageom = "BGC", ttr = "TTR", totreads = "Total reads", nzmean = "Nonzero mean", aBGC = "aBGC",
  quantile = "Quantile", zinfnb = "ZINB", midlc = "MidLC", nonorm = "Unnormalized")
# Every panel already shows unnormalized counts beside the method, so Unnormalized is not a choice of its own.
norm_shown_methods <- setdiff(names(normalization_method_short), "nonorm")

normalization_site_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Normalization",
    selectizeInput("norm_methods", "Show methods", stats::setNames(norm_shown_methods, normalization_method_labels[norm_shown_methods]),
                   selected = "betageom", multiple = TRUE,
                   options = list(placeholder = "Choose methods to show\u2026", plugins = list("remove_button"))),
    helpText("Each run normalizes with every method; choose the ones to look at. BGC is suggested for condition comparisons (see Methods)."),
    tags$details(class = "ytab-more-options", tags$summary("Settings"),
      numericInput("norm_trim", "Trim per tail (TTR and BGC scaling)", .05, min = 0, max = .45, step = .01),
      numericInput("norm_target", "Target scale", 100, min = 1, step = 1)),
    actionButton("norm_run", "Run normalization", class = "btn-primary"),
    ytab_reference_projects_ui("norm", display_modes = FALSE),
    tags$hr(),
    uiOutput("norm_result_selector"),
    radioButtons("norm_plot", "View", c("Each normalization, before and after" = "native", "Replicate agreement (methods side by side)" = "pairs")),
    conditionalPanel("input.norm_plot == 'pairs'",
      radioButtons("norm_pair_plot", "Show", c("Shift and spread at shared sites" = "sites", "Gene-level shift and rank agreement" = "gene_summary",
                                               "Gene log2 ratios, each pair" = "gene_ratios")),
      ytab_plot_downloads_ui("norm", c(pairs = "Download pair diagnostics (TSV)", gene_agreement = "Download gene-level agreement (TSV)"))),
    downloadButton("norm_download_stats", "Download per-sample statistics (TSV)")),
  main = tagList(
    # The window height sizes the panels so that all chosen methods fit on screen.
    tags$script(HTML("(function(){var t;function send(){if(window.Shiny&&Shiny.setInputValue)Shiny.setInputValue('norm_viewport_h',window.innerHeight);}
      $(document).on('shiny:connected',send);$(window).on('resize',function(){clearTimeout(t);t=setTimeout(send,300);});})();")),
    uiOutput("norm_status"),
    conditionalPanel("input.norm_plot == 'native'", uiOutput("norm_native_panels")),
    conditionalPanel("input.norm_plot == 'pairs'",
      ytab_plot_card("Replicate agreement", uiOutput("norm_plot_frame"), description = uiOutput("norm_plot_caption"))),
    ytab_table_card("Per-sample statistics", DT::DTOutput("norm_stats_table")),
    ytab_methods_drawer(
      tags$p("Counts at each insertion site are normalized per sample, within each library (one background and pool), so that libraries sequenced to different depths can be compared. Every method starts from the same raw counts; zeros stay zero."),
      tags$p(tags$b("Which method? "), "In the pre-registered Hermes benchmarks, beta-geometric correction (BGC, then TTR scaling) with recentering kept false discoveries lowest: 7.0% versus 12.4% with TTR in the two-condition test. The multi-condition and genetic-interaction tests were validated with BGC. BGC compresses the upper tail, so it can attenuate a real enrichment carried by a few very abundant sites; compare with TTR when that matters. Unnormalized counts are not suitable when libraries differ in depth."),
      tags$p(tags$b("Reading the views. "), "Each method is shown with the plot that shows what it does, before (unnormalized) and after: a scaling method makes one statistic equal across samples (the shape of the counts is unchanged); BGC and adaptive BGC reshape the counts toward a geometric distribution (quantile–quantile plot, ranked counts); quantile normalization makes the count distributions identical; MidLC thins libraries to a common saturation. Replicate agreement compares methods directly: matched controls from the same pool should differ little after normalization (median log2 ratio near 0, narrow spread). The per-sample statistics are the standard Tn-seq library statistics (see Reference)."),
      tags$p(tags$b("MidLC common saturation. "), "MidLC is the read depth at which a library shows half of its insertion sites; depth divided by MidLC says how far past that point it was sequenced. Each library is randomly thinned (reads removed, never added) to the smallest of these ratios among the libraries normalized together, here the libraries of one background and pool. It removes differences in sequencing saturation, not in insertion preference. A strongly selected library has few sites with many reads, looks highly saturated and is thinned the most, so MidLC suits unselected libraries (such as parents for essentiality) better than treated–control comparisons."),
      uiOutput("norm_method_guide"),
      ytab_reference("TTR, BGC and the other site methods are adapted from the normalization methods of TRANSIT (DeJesus et al. 2015, PLoS Comput Biol 11:e1004401); the quantile–quantile plot against a geometric distribution, the ranked counts and the per-sample statistics follow TRANSIT's quality-control views."),
      ytab_reference("MidLC depth correction follows Gale et al. 2020 (G3 10:3859–3870, doi:10.1534/g3.120.401595), with the target set to the smallest saturation among the libraries normalized together."),
      tags$p("Command: ", tags$code("ytab normalize-project"), ". The normalized counts for condition comparisons are written to the run's ", tags$code("experiment/methods"), " folder."))))

normalization_module_ui <- function() navset_tab(id = "normalization_tabs",
  nav_panel("Normalization", value = "site", normalization_site_ui()))

norm_plot_data <- function(run_dir) {
  read <- function(...) ytab_read_analysis_tsv(file.path(run_dir, ...))
  list(scale = read("plots", "effective_scale.tsv"), distributions = read("plots", "distributions.tsv"),
       pairs = read("experiment", "comparison_diagnostics.tsv"), stats = read("plots", "native_stats.tsv"),
       scaling = read("plots", "native_scaling.tsv"), qq = read("plots", "native_qq.tsv"), ranked = read("plots", "native_ranked.tsv"),
       saturation = read("plots", "native_saturation.tsv"), midlc = read("plots", "native_midlc.tsv"),
       abgc = read("experiment", "methods", "aBGC", "factors.tsv"), gene_pairs = read("plots", "pair_gene_agreement.tsv"))
}
# Per-gene log2 ratios of every pair (large): read only for the view that shows them.
norm_pair_genes <- function(run_dir) ytab_read_analysis_tsv(file.path(run_dir, "plots", "pair_genes.tsv.gz"))

norm_method_label <- function(x) {
  x <- as.character(x); out <- unname(normalization_method_labels[x]); out[is.na(out)] <- x[is.na(out)]
  factor(out, levels = unique(c(normalization_method_labels, out)))
}
norm_method_short_label <- function(x) {
  x <- as.character(x); out <- unname(normalization_method_short[x]); out[is.na(out)] <- x[is.na(out)]
  factor(out, levels = unique(c(normalization_method_short, out)))
}

# Plot data restricted to the ticked methods; unnormalized counts are always kept as the reference.
norm_filter <- function(data, methods) lapply(data, function(d) if (is.data.frame(d) && "method" %in% names(d)) d[d$method %in% c(methods, "nonorm"), , drop = FALSE] else d)

norm_scaling_methods <- c("ttr", "totreads", "nzmean", "zinfnb")
# Compact numbers on every axis (1K, 10K, 1M) so large counts never crowd the axis title, whatever the depth.
norm_number <- scales::label_number(scale_cut = scales::cut_short_scale())
norm_statistic_short <- c(ttr = "Density × trimmed mean", totreads = "Total reads", nzmean = "Nonzero mean", zinfnb = "Fitted NB mean")
norm_native_caption <- function(method, data = NULL) {
  if (method %in% norm_scaling_methods) return("A scaling method multiplies every count of a sample by one factor, chosen so that one statistic becomes equal across the samples of a library. Left: that statistic before; right: after. The shape of each sample's counts is unchanged.")
  switch(method,
    betageom = "BGC fits a beta-geometric distribution to each sample's counts over all its sites, empty sites included, and moves every count to the same quantile of a geometric distribution; this compresses the jackpot tail. The result is then scaled by TTR. Top: quantiles of the nonzero counts against a geometric distribution fitted to them (on the dashed line = geometric). Bottom: counts in ranked order. Left: unnormalized; right: after.",
    aBGC = {
      chosen <- if (is.data.frame(data$abgc) && nrow(data$abgc) && "abgc_selected_trim" %in% names(data$abgc))
        paste(sprintf("%s %.0f%%", data$abgc$display_name, 100 * as.numeric(data$abgc$abgc_selected_trim)), collapse = "; ") else ""
      paste0("Adaptive BGC does the same mapping as BGC, but it fits only the nonzero counts and does not fix how much of the tail is ignored in the fit: it tries trimming 0–30% of the counts and keeps, for each sample, the trim whose result is closest to a geometric distribution. Use it when the plain BGC fit is pulled by extreme counts.",
             if (nzchar(chosen)) paste0(" Trim chosen: ", chosen, ".") else "", " Top: quantiles against a fitted geometric; bottom: ranked counts. Left: unnormalized; right: after.")
    },
    quantile = "Quantile normalization gives every sample of a library the same distribution of nonzero counts. Counts at each percentile of each sample's insertion sites, unnormalized (left) and after (right).",
    midlc = "Expected insertion sites against reads as each library is thinned (exact rarefaction from its own counts). Open circles: MidLC, the depth showing half of the library's sites; filled circles: the library as sequenced; triangles: after thinning to the smallest depth/MidLC ratio of the library group. Libraries are only thinned, never added to.",
    nonorm = "The raw counts: quantiles of the nonzero counts against a geometric distribution fitted to them (strong upward curvature = jackpots, skew) and counts in ranked order.",
    "")
}

# One colour per sample, the same in every panel; the samples are listed once, in a legend above the panels.
norm_sample_colours <- function(data) {
  samples <- unique(c(data$stats$display_name[data$stats$method == "nonorm"], data$stats$display_name, data$scaling$display_name, data$saturation$display_name))
  samples <- samples[!is.na(samples)]
  if (!length(samples)) return(character())
  stats::setNames(scales::hue_pal()(length(samples)), samples)
}
norm_sample_legend <- function(colours) {
  if (!length(colours)) return(NULL)
  tags$div(class = "ytab-norm-legend", style = "display: flex; flex-wrap: wrap; gap: 4px 16px; margin: 0 0 8px; font-size: .9em;",
    lapply(names(colours), function(n) tags$span(style = "display: inline-flex; align-items: center; gap: 6px;",
      tags$span(style = sprintf("width: 14px; height: 4px; border-radius: 2px; background: %s; display: inline-block;", colours[[n]])), n)))
}

norm_native_ggplot <- function(data, method) {
  base <- ggplot2::theme_minimal(base_size = 12); label <- as.character(norm_method_short_label(method))
  colours <- ggplot2::scale_colour_manual(values = norm_sample_colours(data), guide = "none")
  stage <- function(m) factor(ifelse(m == "nonorm", "Unnormalized", label), levels = unique(c("Unnormalized", label)))
  tip <- function(...) sprintf(...)
  if (method %in% norm_scaling_methods) {
    d <- data$scaling[data$scaling$method == method, , drop = FALSE]; req(nrow(d) > 0)
    long <- rbind(data.frame(d[, c("display_name", "role", "statistic", "factor")], panel = "Unnormalized", value = d$before),
                  data.frame(d[, c("display_name", "role", "statistic", "factor")], panel = label, value = d$after))
    long$panel <- factor(long$panel, levels = unique(long$panel))
    return(ggplot2::ggplot(long, ggplot2::aes(display_name, value, colour = display_name,
        text = tip("%s (%s)<br>%s<br>%.4g<br>factor %.4g", display_name, role, panel, value, factor))) +
      ggplot2::geom_point(size = 3.5) + ggplot2::scale_y_continuous(labels = norm_number) + ggplot2::facet_wrap(~panel, scales = "free_y") + colours + base +
      ggplot2::labs(x = "Sample", y = norm_statistic_short[[method]], colour = "Sample") +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()))   # samples identified by colour (legend)
  }
  if (method %in% c("betageom", "aBGC", "nonorm")) {
    # One small plot per panel, each with its own axes: quantiles against a fitted geometric, and the ranked counts,
    # unnormalized and after (only unnormalized for the unnormalized view).
    qq <- data$qq[data$qq$method %in% unique(c("nonorm", method)), , drop = FALSE]
    rk <- data$ranked[data$ranked$method %in% unique(c("nonorm", method)), , drop = FALSE]; req(nrow(qq) > 0, nrow(rk) > 0)
    stages <- unique(c("nonorm", method)); title <- function(kind, m) paste0(kind, ": ", if (m == "nonorm") "Unnormalized" else label)
    qq_plot <- function(m) { d <- qq[qq$method == m, , drop = FALSE]; r <- range(c(d$expected, d$observed), na.rm = TRUE)
      ggplot2::ggplot(d, ggplot2::aes(expected, observed, colour = display_name, group = display_name,
          text = tip("%s<br>quantile %.4f<br>geometric %.3g<br>observed %.3g", display_name, probability, expected, observed))) +
        ggplot2::geom_line(data = data.frame(x = r, y = r), ggplot2::aes(x, y), inherit.aes = FALSE, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_line() + ggplot2::geom_point(size = .9) +
        ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 4), labels = norm_number) + ggplot2::scale_y_log10(breaks = scales::breaks_log(n = 5), labels = norm_number) +
        colours + base + ggplot2::labs(title = title("Quantiles vs geometric", m), x = "Geometric quantile (log)", y = "Count (log)", colour = "Sample") }
    rk_plot <- function(m) { d <- rk[rk$method == m, , drop = FALSE]
      ggplot2::ggplot(d, ggplot2::aes(rank, value, colour = display_name, group = display_name, text = tip("%s<br>rank %d<br>count %.3g", display_name, rank, value))) +
        ggplot2::geom_line() + ggplot2::geom_point(size = .9) +
        ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 4), labels = norm_number) + ggplot2::scale_y_log10(breaks = scales::breaks_log(n = 5), labels = norm_number) +
        colours + base + ggplot2::labs(title = title("Ranked counts", m), x = "Rank (log)", y = "Count (log)", colour = "Sample") }
    parts <- c(lapply(stages, qq_plot), lapply(stages, rk_plot))   # row 1: quantiles; row 2: ranked
    attr(parts, "grid") <- c(rows = 2L, cols = length(stages))
    return(parts)
  }
  if (identical(method, "quantile")) {
    d <- data$distributions[data$distributions$method %in% c("nonorm", "quantile"), , drop = FALSE]; req(nrow(d) > 0)
    d$panel <- stage(d$method)
    return(ggplot2::ggplot(d, ggplot2::aes(percentile, value, colour = display_name, group = display_name,
        text = tip("%s<br>%s<br>percentile %d: %.3g", display_name, panel, percentile, value))) +
      ggplot2::geom_line() + ggplot2::scale_y_log10(labels = norm_number) + ggplot2::facet_wrap(~panel, scales = "free_y") + colours + base +
      ggplot2::labs(x = "Percentile of the sample's insertion sites", y = "Count (log)", colour = "Sample"))
  }
  if (identical(method, "midlc")) {
    curve <- data$saturation; pts <- data$midlc; req(nrow(curve) > 0, nrow(pts) > 0)
    marks <- rbind(data.frame(pts[, c("display_name", "role")], what = "MidLC (half of the sites)", reads = pts$midlc_reads, sites = pts$midlc_sites),
                   data.frame(pts[, c("display_name", "role")], what = "As sequenced", reads = pts$raw_reads, sites = pts$raw_sites),
                   data.frame(pts[, c("display_name", "role")], what = "After thinning", reads = pts$thinned_reads, sites = pts$thinned_sites))
    marks$text <- tip("%s<br>%s<br>%s reads, %s sites", marks$display_name, marks$what, format(round(marks$reads), big.mark = ","), format(round(marks$sites), big.mark = ","))
    marks$what <- factor(marks$what, levels = c("MidLC (half of the sites)", "As sequenced", "After thinning"))
    return(ggplot2::ggplot(curve, ggplot2::aes(reads, sites, colour = display_name, group = display_name,
        text = tip("%s<br>%s reads: %s sites", display_name, format(round(reads), big.mark = ","), format(round(sites), big.mark = ",")))) +
      ggplot2::geom_line() + ggplot2::geom_point(data = marks, ggplot2::aes(reads, sites, colour = display_name, shape = what, text = text), size = 3, inherit.aes = FALSE, show.legend = FALSE) +
      ggplot2::scale_shape_manual(values = c(1, 16, 17)) + ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 4), labels = norm_number) +
      ggplot2::scale_y_continuous(labels = norm_number) + colours + base +
      ggplot2::labs(x = "Reads (log scale)", y = "Insertion sites", colour = "Sample", shape = NULL))
  }
  NULL
}

# One plotly figure from a method's plot(s): two plots are stacked (side by side for unnormalized counts), each keeping its
# own axis titles; one legend. Vertical axis titles get extra room from the tick labels.
norm_native_plotly <- function(parts, method, height = 500) {
  finish <- function(out) {   # one shared legend above the panels; plotly sizes the axis margins
    out$x$layout$showlegend <- FALSE
    for (a in names(out$x$layout)[grepl("^[xy]axis[0-9]*$", names(out$x$layout))]) {
      out$x$layout[[a]]$automargin <- TRUE
      if (is.list(out$x$layout[[a]]$title)) out$x$layout[[a]]$title$standoff <- 6
    }
    out
  }
  if (inherits(parts, "ggplot")) {
    out <- plotly::plotly_build(ytab_ggplotly(parts))
    out$x$layout$annotations <- lapply(out$x$layout$annotations, function(a) {   # facet plots: vertical axis title off the tick labels
      if (identical(as.numeric(a$textangle %||% 0), -90)) a$xshift <- (a$xshift %||% 0) - 14
      a })
    out$x$layout$margin$l <- (out$x$layout$margin$l %||% 40) + 14
    return(finish(out))
  }
  grid <- attr(parts, "grid") %||% c(rows = length(parts), cols = 1L)
  figs <- lapply(seq_along(parts), function(i) { f <- ytab_ggplotly(parts[[i]] + ggplot2::labs(title = NULL))
    if (i > 1L) f$x$data <- lapply(f$x$data, function(tr) { tr$showlegend <- FALSE; tr })
    f })
  # Room between rows for tick labels, the x title and the next panel title (~100 px whatever the panel height).
  gap_y <- if (grid[["rows"]] > 1L) min(.45, 100 / max(height, 1)) else 0
  out <- plotly::subplot(figs, nrows = grid[["rows"]], titleX = TRUE, titleY = TRUE, margin = c(.05, .05, gap_y / 2, gap_y / 2))
  out <- plotly::plotly_build(out)
  # every column shares the y title of the first one, so only the first column shows it
  y_axes <- names(out$x$layout)[grepl("^yaxis[0-9]*$", names(out$x$layout))]
  y_axes <- y_axes[order(suppressWarnings(as.integer(sub("yaxis", "", y_axes))), na.last = FALSE)]
  for (i in seq_along(y_axes)) if ((i - 1L) %% grid[["cols"]] > 0L) out$x$layout[[y_axes[[i]]]]$title <- list(text = "")
  # panel titles above each small plot (subplot drops ggplot titles)
  titles <- vapply(parts, function(p) p$labels$title %||% "", "")
  axes <- names(out$x$layout)[grepl("^xaxis[0-9]*$", names(out$x$layout))]
  axes <- axes[order(suppressWarnings(as.integer(sub("xaxis", "", axes))), na.last = FALSE)]
  out$x$layout$annotations <- c(out$x$layout$annotations, lapply(seq_along(titles), function(i) {
    xa <- out$x$layout[[axes[[i]]]]; ya <- out$x$layout[[sub("^x", "y", xa$anchor %||% "y")]] %||% out$x$layout[[paste0("yaxis", if (i > 1L) i else "")]]
    list(text = titles[[i]], x = mean(xa$domain), y = ya$domain[[2L]], xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom",
         showarrow = FALSE, font = list(size = 12)) }))
  finish(out)
}

# Downloads carry the sample legend (in the app it is shown once above the panels): below a single plot, once below a grid.
norm_with_legend <- function(p) p + ggplot2::guides(colour = ggplot2::guide_legend(title = "Sample", ncol = 2)) + ggplot2::theme(legend.position = "bottom")

norm_save_native <- function(file, parts, method, width, height) {
  if (inherits(parts, "ggplot")) return(ytab_save_ggplot(file, norm_with_legend(parts), width = width, height = height))
  grid <- attr(parts, "grid") %||% c(rows = length(parts), cols = 1L); rows <- grid[["rows"]]; cols <- grid[["cols"]]
  table <- ggplot2::ggplotGrob(ytab_style_ggplot(norm_with_legend(parts[[1L]])))
  legend <- table$grobs[[which(grepl("guide-box", table$layout$name) & !vapply(table$grobs, inherits, logical(1), "zeroGrob"))[[1L]]]]
  legend_height <- grid::convertHeight(grid::grobHeight(legend), "in", valueOnly = TRUE) + .2
  if (grepl("\\.pdf$", file, ignore.case = TRUE)) ytab_pdf_device()(file, width = width, height = height + legend_height) else grDevices::png(file, width = width, height = height + legend_height, units = "in", res = 600)
  on.exit(grDevices::dev.off())
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(rows + 1L, cols, heights = grid::unit(c(rep(height / rows, rows), legend_height), "in"))))
  for (i in seq_along(parts)) print(ytab_style_ggplot(parts[[i]]), vp = grid::viewport(layout.pos.row = (i - 1L) %/% cols + 1L, layout.pos.col = (i - 1L) %% cols + 1L))
  grid::pushViewport(grid::viewport(layout.pos.row = rows + 1L, layout.pos.col = seq_len(cols))); grid::grid.draw(legend); grid::popViewport()
  invisible(file)
}

# Replicate agreement: the one view where methods are compared side by side. sites: shift and spread of log2 ratios at
# the sites both samples contain; gene_summary: per pair, the median gene log2 ratio (shift) and the Spearman
# correlation of gene read sums (agreement no scaling can change); gene_ratios: every gene's log2 ratio, per pair.
norm_pair_kind <- function(kind) ifelse(kind == "control_control", "control vs control", "treated vs control")
norm_ggplot <- function(data, choice = "sites", genes = NULL) {
  base <- ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  if (identical(choice, "gene_ratios")) {
    g <- genes[genes$method %in% unique(data$gene_pairs$method), , drop = FALSE]; validate(need(nrow(g) > 0, "No gene-level ratios in this run."))
    g$method_id <- g$method; g$method <- norm_method_short_label(g$method)
    g$pair <- factor(sprintf("%s (%s)", g$comparison_id, norm_pair_kind(g$kind)))
    fills <- stats::setNames(YTAB_NORM_COLOURS[unique(g$method_id)], as.character(norm_method_short_label(unique(g$method_id))))
    return(ggplot2::ggplot(g, ggplot2::aes(method, log2_ratio, fill = method)) + ggplot2::geom_hline(yintercept = 0, linewidth = .5) +
      ggplot2::geom_violin(scale = "width", linewidth = .4, colour = "black") +
      ggplot2::geom_boxplot(width = .15, outlier.shape = NA, fill = "white", linewidth = .4) +
      ggplot2::facet_wrap(~pair, ncol = min(3L, nlevels(g$pair))) + ggplot2::scale_fill_manual(values = fills, guide = "none") +
      ggplot2::coord_cartesian(ylim = c(-5, 5)) + ggplot2::labs(x = NULL, y = "Gene log2 ratio, sample A / B") + base +
      ggplot2::theme(strip.text = ggplot2::element_text(face = "bold")))
  }
  if (identical(choice, "gene_summary")) {
    d <- data$gene_pairs; validate(need(nrow(d) > 0, "This run has no gene-level agreement; run normalization again to add it."))
    long <- rbind(data.frame(d[, c("method", "comparison_id", "kind")], metric = "Median gene log2 ratio (shift)", value = d$median_gene_log2_ratio),
                  data.frame(d[, c("method", "comparison_id", "kind")], metric = "Rank agreement of genes (Spearman)", value = d$spearman_genes))
    zero <- data.frame(metric = "Median gene log2 ratio (shift)")
  } else {
    d <- data$pairs
    validate(need(nrow(d) > 0 && "median_log2_ratio_shared" %in% names(d), "This run normalized a list of samples (for multi-condition tests) without matched pairs, so it has no replicate agreement. A run made from the comparison design (Run normalization) has it; earlier design runs are under Earlier runs."))
    d <- d[!is.na(d$median_log2_ratio_shared), , drop = FALSE]
    long <- rbind(data.frame(d[, c("method", "comparison_id", "kind")], metric = "Median log2 ratio (shift)", value = d$median_log2_ratio_shared),
                  data.frame(d[, c("method", "comparison_id", "kind")], metric = "Spread of log2 ratios (IQR)", value = d$log2_ratio_iqr_shared))
    zero <- data.frame(metric = "Median log2 ratio (shift)")
  }
  # Unnormalized first, then the methods; one line per pair joins its values, so the change each method makes reads left to right.
  long$method <- norm_method_short_label(long$method); long$method <- factor(long$method, unique(c("Unnormalized", levels(long$method))))
  long$method <- droplevels(long$method); long$kind <- norm_pair_kind(long$kind)
  ggplot2::ggplot(long, ggplot2::aes(method, value, colour = kind, group = comparison_id, text = sprintf("%s<br>%s<br>%s: %.3f", comparison_id, method, metric, value))) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), data = zero, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line(linewidth = .7, alpha = .7) + ggplot2::geom_point(size = 3) + ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c(`control vs control` = YTAB_FOCAL, `treated vs control` = YTAB_ACCENT)) + base +
    ggplot2::labs(x = NULL, y = NULL, colour = "Pair")
}

normalization_server <- function(input, output, session, repo_root, active, active_project_path, jobs, python_bin,
                                 projects = reactive(data.frame()), reference_state = NULL) {
  if (is.null(reference_state)) reference_state <- ytab_reference_state()
  ytab_reference_server("norm", input, output, session, repo_root, active, projects, reference_state, display_modes = FALSE)
  reference_normalization <- reactive(
    ytab_reference_normalization(ytab_reference_metas(reference_state(), repo_root), norm_plot_data))
  runs_dir <- reactive(file.path(active()$project_root, "normalization"))
  runs <- reactivePoll(3000, session,
    checkFunc = function() { x <- ytab_analysis_runs(tryCatch(runs_dir(), error = function(e) "")); paste(x$path, x$modified, collapse = "|") },
    valueFunc = function() { x <- ytab_analysis_runs(tryCatch(runs_dir(), error = function(e) "")); x[file.exists(file.path(x$path, "plots", "effective_scale.tsv")), , drop = FALSE] })
  design_ready <- reactive(file.exists(file.path(active()$project_root, "config", "comparison_design.csv")))

  # The latest run is shown; earlier runs (other settings or an older design) stay reachable under "Earlier runs".
  run_label <- function(path, modified) {
    m <- ytab_read_manifest(path) %||% list(); methods <- setdiff(unlist(m$methods), "nonorm")
    sprintf("%s · %d methods · trim %g%%, target %g", format(modified, "%d %b %Y %H:%M"), length(methods),
            100 * as.numeric(m$trim %||% NA), as.numeric(m$target %||% NA))
  }
  output$norm_result_selector <- renderUI({
    x <- runs()
    if (!nrow(x)) return(tags$p(class = "text-muted", "No normalization has been run for this project yet."))
    m <- ytab_read_manifest(x$path[[1L]]) %||% list()
    design <- file.path(active()$project_root, "config", "comparison_design.csv")
    stale <- file.exists(design) && file.mtime(design) > x$modified[[1L]]
    tagList(
      tags$p(class = "text-muted", sprintf("Latest normalization: %s, %s samples, trim %g%%, target %g.", format(x$modified[[1L]], "%d %b %Y %H:%M"),
                                          m$samples %||% "?", 100 * as.numeric(m$trim %||% NA), as.numeric(m$target %||% NA))),
      if (stale) tags$p(class = "text-warning", "The comparison design changed after this run; run normalization again to include the change."),
      if (nrow(x) > 1L) tags$details(class = "ytab-more-options", tags$summary(sprintf("Earlier runs (%d)", nrow(x) - 1L)),
        selectInput("norm_result", NULL, stats::setNames(x$path, vapply(seq_len(nrow(x)), function(i)
          paste0(run_label(x$path[[i]], x$modified[[i]]), if (i == 1L) " (latest)" else ""), "")))))
  })

  observeEvent(input$norm_run, {
    if (jobs$job_is_running()) { showNotification("Wait for the current job to finish.", type = "warning"); return() }
    if (!design_ready()) { showNotification("Create the comparison design (Condition Comparisons → Design) first.", type = "error"); return() }
    args <- c(file.path(repo_root, "scripts", "local", "ytab_normalize_project.py"), "--project-config", active_project_path(),
              "--trim", as.character(input$norm_trim %||% .05),
              "--target", as.character(input$norm_target %||% 100), "--output", ytab_new_run_dir(runs_dir(), "norm"))
    tryCatch({
      jobs$start_job(python_bin(), args, active_project_path(), "site_normalization", wd = repo_root, selected_items = norm_shown_methods)
      showNotification("Normalization started. The result appears here when it finishes.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error"))
  }, ignoreInit = TRUE)

  selected_run <- reactive({ req(nrow(runs()) > 0); path <- input$norm_result %||% runs()$path[[1L]]; req(path %in% runs()$path); path })
  data <- reactive({
    run <- selected_run(); d <- norm_plot_data(run)
    if (nrow(d$scale) && !nrow(d$stats)) {   # runs made before the per-method views: add their summaries once
      system2(python_bin(), c(file.path(repo_root, "scripts", "local", "ytab_normalize_project.py"), "--refresh-summaries", run), stdout = FALSE, stderr = FALSE)
      d <- norm_plot_data(run)
    }
    # Reference samples join as samples of their own, so every view and the legend pick them up unchanged.
    ytab_reference_normalization_bind(d, reference_normalization())
  })
  shown <- reactive(norm_filter(data(), input$norm_methods %||% character()))
  missing_methods <- reactive(setdiff(input$norm_methods %||% character(), unique(data()$scale$method)))   # nonorm is never a choice
  manifest <- reactive(ytab_read_manifest(selected_run()))

  output$norm_status <- renderUI({
    if (!design_ready()) return(ytab_empty_state("Comparison design needed",
      "Normalization uses the samples in the project's comparison design (each treated sample with its matched control). Create it in Condition Comparisons → Design."))
    if (!nrow(runs())) return(ytab_empty_state("No normalization yet", "Choose methods and run normalization. BGC and TTR are preselected."))
    missing <- missing_methods()
    if (length(missing)) tags$p(class = "text-warning", sprintf("Not in this run: %s. Run normalization again to add every method.",
                                                                paste(as.character(norm_method_short_label(missing)), collapse = ", ")))
  })
  pair_choice <- reactive(input$norm_pair_plot %||% "sites")
  pair_genes <- reactive(norm_filter(list(g = norm_pair_genes(selected_run())), input$norm_methods %||% character())$g)
  pair_plot <- reactive({
    d <- shown(); validate(need(length(input$norm_methods), "Choose at least one method to show."))
    if (identical(pair_choice(), "sites")) req(nrow(d$pairs) > 0)
    norm_ggplot(d, pair_choice(), if (identical(pair_choice(), "gene_ratios")) pair_genes())
  })
  output$norm_plot_caption <- renderUI(switch(pair_choice(),
    gene_summary = "Each point is one matched pair of samples. Normalized reads are summed over each gene's coding sequence (genes with at least 3 catalogue sites and reads in both samples). Shift: the median gene log2 ratio, near 0 when the two samples are on one scale. Rank agreement: the Spearman correlation of the gene sums, which no scaling can change; a method that reshapes counts (BGC) can raise it.",
    gene_ratios = "Every gene's log2 ratio between the two samples of each pair, under each method (violin: distribution; box: median and quartiles; axis cut at ±5). Controls from the same pool should centre on 0 with a narrow distribution.",
    "Each point is one matched pair of samples, compared at the insertion sites both contain. After a good normalization, controls from the same pool have a shift near 0 and a narrow spread; treated–control pairs also contain the treatment's effect."))
  output$norm_plot_frame <- renderUI({ req(nrow(runs()) > 0); plotly::plotlyOutput("norm_plotly", height = "520px") })
  output$norm_plotly <- plotly::renderPlotly(ytab_ggplotly(pair_plot()))

  # One card per ticked method, each with its own plot and axes.
  # Panels fit the window: one column for one method, two for up to four, three beyond; the window height is shared by the rows.
  native_layout <- reactive({
    methods <- intersect(norm_shown_methods, setdiff(input$norm_methods %||% character(), missing_methods()))
    n <- length(methods); columns <- if (n <= 1L) 1L else if (n <= 4L) 2L else 3L; rows <- max(1L, ceiling(n / columns))
    viewport <- suppressWarnings(as.numeric(input$norm_viewport_h %||% 900)); if (!length(viewport) || !is.finite(viewport)) viewport <- 900
    # every panel shrinks so that all of them fit in the window (no fixed minimum size)
    legend_rows <- ceiling(length(norm_sample_colours(data())) / 4)   # the shared sample legend wraps at about four names per row
    list(methods = methods, columns = columns, height = as.integer(max(150, (viewport - 170 - 24 * legend_rows) / rows - 95)))
  })
  output$norm_native_panels <- renderUI({
    req(nrow(runs()) > 0); layout <- native_layout(); methods <- layout$methods; n <- length(methods)
    if (!n) return(ytab_empty_state("No method ticked", "Tick a method on the left to see what it does to each sample."))
    tagList(norm_sample_legend(norm_sample_colours(data())),
    tags$div(class = "ytab-norm-panels", style = sprintf("display: grid; grid-template-columns: repeat(%d, minmax(0, 1fr)); gap: 12px;", layout$columns),
      lapply(methods, function(m) ytab_plot_card(normalization_method_labels[[m]],
        plotly::plotlyOutput(paste0("norm_native_", m), height = paste0(layout$height, "px")),
        description = tags$details(open = if (n == 1L) NA else NULL, tags$summary("What this shows"), norm_native_caption(m, data())),
        controls = tags$div(class = "ytab-actions", downloadButton(paste0("norm_native_png_", m), "PNG", class = "btn-sm"),
                            downloadButton(paste0("norm_native_pdf_", m), "PDF", class = "btn-sm"))))))
  })
  for (method in norm_shown_methods) local({
    m <- method
    output[[paste0("norm_native_", m)]] <- plotly::renderPlotly({
      layout <- native_layout(); req(m %in% layout$methods)   # redraw at the new size whenever the panel layout changes
      p <- norm_native_ggplot(data(), m); req(p); out <- norm_native_plotly(p, m, layout$height)
      out$x$layout$height <- layout$height; out$x$layout$autosize <- TRUE
      out
    })
    for (format in c("png", "pdf")) local({
      ext <- format
      output[[paste0("norm_native_", ext, "_", m)]] <- downloadHandler(
        filename = function() paste0(active()$project_id, "_normalization_", m, "_", basename(selected_run()), ".", ext),
        content = function(file) norm_save_native(file, norm_native_ggplot(data(), m), m, width = 11, height = if (m %in% c("betageom", "aBGC")) 9 else 6))
    })
  })

  output$norm_stats_table <- DT::renderDT({
    d <- shown()$stats; req(nrow(d) > 0)
    d$method <- as.character(norm_method_short_label(d$method))
    for (k in c("density", "mean", "nzmean", "nzmedian", "max", "total", "skew", "kurtosis")) d[[k]] <- signif(d[[k]], 4)
    d <- d[, c("display_name", "role", "method", "density", "mean", "nzmean", "nzmedian", "max", "total", "skew", "kurtosis")]
    names(d) <- c("Sample", "Role", "Method", "Density", "Mean", "Nonzero mean", "Nonzero median", "Max", "Total", "Skew", "Kurtosis")
    ytab_datatable(d, page_length = 10L)
  })
  output$norm_method_guide <- renderUI({
    catalogue <- tryCatch(jsonlite::fromJSON(paste(system2(python_bin(), c(file.path(repo_root, "scripts", "local", "ytab_normalize.py"), "--list-methods"),
                                                           stdout = TRUE, stderr = FALSE), collapse = "\n")), error = function(e) NULL)
    if (!is.data.frame(catalogue)) return(NULL)
    catalogue <- catalogue[catalogue$status == "available", c("name", "description", "may_suit", "caution"), drop = FALSE]
    names(catalogue) <- c("Method", "What it does", "Use when", "Limitations")
    tags$div(tags$h5("Available methods"), DT::datatable(catalogue, rownames = FALSE, options = list(dom = "t", pageLength = 20)))
  })

  for (format in c("png", "pdf")) local({
    ext <- format
    output[[paste0("norm_download_", ext)]] <- downloadHandler(filename = function() paste0(active()$project_id, "_normalization_agreement_", basename(selected_run()), ".", ext),
      content = function(file) ytab_save_ggplot(file, pair_plot(), width = 11, height = 7))
  })
  output$norm_download_stats <- downloadHandler(filename = function() paste0(active()$project_id, "_normalization_statistics_", basename(selected_run()), ".tsv"),
    content = function(file) file.copy(file.path(selected_run(), "plots", "native_stats.tsv"), file))
  output$norm_download_pairs <- downloadHandler(filename = function() paste0(active()$project_id, "_normalization_pairs_", basename(selected_run()), ".tsv"),
    content = function(file) file.copy(file.path(selected_run(), "experiment", "comparison_diagnostics.tsv"), file))
  output$norm_download_gene_agreement <- downloadHandler(filename = function() paste0(active()$project_id, "_normalization_gene_agreement_", basename(selected_run()), ".tsv"),
    content = function(file) file.copy(file.path(selected_run(), "plots", "pair_gene_agreement.tsv"), file))
  invisible(list(runs = runs, selected_run = selected_run, data = data, manifest = manifest))
}
