# Treated vs Control: three ways to call genes depleted or enriched in treated libraries relative to their controls.
# The local-abundance z-score (default) uses the Fitness run and results views (ytab fitness); resampling and the
# control-conditioned Gumbel run on a normalization run (ytab treated-vs-control) and have their own views here.

tvc_method_choices <- c("Local-abundance z-score" = "zscore",
                        "Resampling" = "resampling",
                        "Control-conditioned Gumbel (depletion only)" = "gumbel")   # validation status: Methods drawer

# Normalizations each method accepts: the Gumbel needs a scalar factor (its count model uses it as the exposure);
# resampling takes the scalar methods and the beta-geometric corrections.
tvc_allowed_methods <- list(resampling = c("betageom", "aBGC", "ttr", "totreads", "nzmean"), gumbel = c("ttr", "totreads", "nzmean"))

# Normalization choices grouped by normalization run, latest first; values are "<run>|<method>". `allowed`: the
# methods a test accepts (NULL: every method of the run, plus Unnormalized when the run has it).
tvc_norm_choices <- function(runs, allowed = NULL) {
  choices <- list()
  for (i in seq_len(nrow(runs))) {
    methods <- c(cmp_run_methods(runs$path[[i]]), if (file.exists(file.path(runs$path[[i]], "experiment", "methods", "nonorm", "normalized_counts.npz"))) "nonorm")
    if (!is.null(allowed)) methods <- intersect(allowed, methods)
    if (!length(methods)) next
    group <- sprintf("%s run of %s", if (i == 1L) "Latest" else "Earlier", format(runs$modified[[i]], "%d %b %Y %H:%M"))
    choices[[group]] <- stats::setNames(paste0(runs$path[[i]], "|", methods), as.character(norm_method_label(methods)))
  }
  choices
}
# The kept choice if still offered, else `preferred` of the latest run that has it, else the first choice.
tvc_norm_default <- function(choices, preferred, current = NULL) {
  values <- unlist(choices, use.names = FALSE)
  if (!is.null(current) && current %in% values) return(current)
  hit <- grep(paste0("\\|", preferred, "$"), values, value = TRUE)
  if (length(hit)) hit[[1L]] else values[[1L]]
}

tvc_methods_text <- function() ytab_methods_drawer(
  tags$p(tags$b("Local-abundance z-score (default). "), "The replicate libraries of each condition are averaged; each gene's log2 fold change between the treated and control means is divided by the noise expected at the gene's abundance. The noise model is fitted to the two controls of each background (running SD of control–control log2 ratios against abundance, SD = m1 + m2·A^m3), and the call threshold is the 99th percentile of control–control |z|. Default input: BGC from the latest normalization run (gene totals of its normalized site counts, expressed per million of the mean library total so the 0.5 pseudocount means the same as with counts per million); total reads (counts per million) and any other normalization method can be chosen. Calls: depleted, enriched or none (one call per gene from the replicate means); the per-pool z-scores are kept as a diagnostic."),
  tags$p(tags$b("Resampling. "), "A permutation test on the gene's normalized site counts in controls versus treated libraries (permuted within each pool when there are several), two-sided, Benjamini–Hochberg FDR. With BGC the treated libraries are recentred (trimmed sum) so that typical genes show no change. ", tags$b("Validation: "), "passed its pre-registered calibration gates with BGC (false discovery rate 6% on simulated Gale-like screens; 4.5% of null genes at p ≤ 0.05). Use it for hit lists with error control; the z-score ranks genes better but its calls are not error-controlled."),
  tags$p(tags$b("Control-conditioned Gumbel. "), "Finds runs of consecutive sites depleted in the treated libraries, compared with the longest runs expected from the controls' count model (exact tail probability; Gumbel fit reported as a diagnostic). Depletion only; needs at least two controls per background; scalar normalizations only (TTR, total reads, nonzero mean). ", tags$b("Validation: "), "passed its pre-registered calibration gates with TTR (false discovery rate 1%, no calls on null screens); depletion only."),
  ytab_reference("The z-score follows Gale et al. 2020 (G3 10:3859–3870). The resampling test is adapted from TRANSIT's resampling (DeJesus et al. 2015, PLoS Comput Biol 11:e1004401); the Gumbel test is adapted from TRANSIT's Gumbel method (DeJesus et al. 2013, Bioinformatics 29:695–703)."),
  tags$p("Commands: ", tags$code("ytab fitness"), " (z-score), ", tags$code("ytab treated-vs-control --method resampling|gumbel"), "."))

tvc_run_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel(uiOutput("tvc_run_title", inline = TRUE),
    uiOutput("tvc_inputs"),
    uiOutput("tvc_comparison_selector"),
    numericInput("tvc_fdr", "FDR", .05, min = .001, max = .25, step = .01),
    actionButton("tvc_run", "Run", class = "btn-primary")),
  main = tvc_methods_text())

# One results view for the three methods: the same controls, display options, plots and gene table, whichever method
# made the result. Only the plots that test something method-specific differ (the z-score's control-control null, the
# p-values of resampling and the Gumbel, the Gumbel's runs).
tvc_results_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Results",
    uiOutput("tvc_result_selector"),
    uiOutput("tvc_group_selector"),
    uiOutput("tvc_plot_selector"),
    uiOutput("tvc_gene_selector"),
    conditionalPanel("input.tvc_plot == 'volcano' || input.tvc_plot == 'ma' || input.tvc_plot == 'scatter'",
      conditionalPanel("input.tvc_results_method != 'gumbel'",
        selectInput("tvc_direction", "Calls shown", c("Depleted and enriched" = "both", "Depleted" = "depleted", "Enriched" = "enriched"))),
      numericInput("tvc_label_n", "Label top hits", 10, min = 0, max = 50, step = 1),
      textAreaInput("tvc_label_genes", "Also label", "", rows = 2, placeholder = "Gene names or IDs, separated by commas or new lines")),
    ytab_plot_downloads_ui("tvc", c(genes = "Download genes (TSV)")),
    ytab_plot_customization_controls("tvc", include_points = TRUE, include_labels = FALSE, open = FALSE)),
  main = tagList(uiOutput("tvc_summary"),
    ytab_plot_card("Treated vs control", uiOutput("tvc_plot_frame"), description = uiOutput("tvc_plot_caption")),
    ytab_table_card("Genes", DT::DTOutput("tvc_table")), tvc_methods_text()))

tvc_run_tab_ui <- function() tagList(
  radioButtons("tvc_method", "Method", tvc_method_choices, "zscore"),
  conditionalPanel("input.tvc_method == 'zscore'", fitness_run_ui()),
  conditionalPanel("input.tvc_method != 'zscore'", tvc_run_ui()))

tvc_results_tab_ui <- function() tagList(
  radioButtons("tvc_results_method", "Method", c(tvc_method_choices, "Compare methods" = "agreement"), "zscore"),
  conditionalPanel("input.tvc_results_method != 'agreement'", tvc_results_ui()),
  conditionalPanel("input.tvc_results_method == 'agreement'", tvc_agreement_ui()))

tvc_plot_choices <- function(method) switch(method %||% "resampling",
  zscore = c("Volcano" = "volcano", "MA plot" = "ma", "Condition vs control" = "scatter", "Control-control z (call threshold)" = "nullz", "Gene values" = "sites"),
  gumbel = c("Volcano" = "volcano", "Observed vs expected run" = "runs", "p-value histogram" = "pvalues", "Gene sites" = "sites"),
  c("Volcano" = "volcano", "MA plot" = "ma", "Condition vs control" = "scatter", "p-value histogram" = "pvalues", "Gene sites" = "sites"))

# Tested genes (and Gumbel family members); the z-score table has no in_family column.
tvc_tested <- function(genes) genes[genes$status %in% "tested" | (if ("in_family" %in% names(genes)) genes$in_family %in% c("True", TRUE) else FALSE), , drop = FALSE]

tvc_call_colours <- c(YTAB_DIRECTION_COLOURS, `depleted run` = YTAB_DIRECTION_COLOURS[["depleted"]], `other direction` = "#ececec")
tvc_call_levels <- function(method) if (identical(method, "gumbel")) c("depleted run", "not significant") else c("depleted", "enriched", "not significant", "other direction")

# Gene labels with leader lines in the interactive plot: plotly annotations (draggable, see tvc_plotly) whose text is
# placed by tvc_repel_labels so labels do not overlap each other or the labelled points.
# Label positions in data units. Works on axes scaled to [0, 1]; a label's box is estimated from its characters for a
# plot of about 800 x 450 px (12 px text). Labels start beside their point, toward the nearer outer edge, then pairs that
# overlap are pushed apart along the axis where they overlap least, and labels are pushed off the labelled points.
tvc_repel_labels <- function(x, y, tag, xr, yr, iterations = 300L) {
  n <- length(x); u0 <- (x - xr[[1L]]) / diff(xr); v0 <- (y - yr[[1L]]) / diff(yr)
  w <- (nchar(tag) * 7 + 10) / 800; h <- 20 / 450
  u <- u0 + ifelse(u0 > .5, 1, -1) * (w / 2 + .03); v <- v0 + .06
  for (it in seq_len(iterations)) {
    moved <- FALSE
    for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
      ox <- (w[[i]] + w[[j]]) / 2 - abs(u[[i]] - u[[j]]); oy <- (h + h) / 2 - abs(v[[i]] - v[[j]])
      if (ox > 0 && oy > 0) { moved <- TRUE
        if (oy < ox) { d <- if (v[[i]] >= v[[j]]) 1 else -1; v[[i]] <- v[[i]] + d * oy / 2 + d * 1e-3; v[[j]] <- v[[j]] - d * oy / 2 - d * 1e-3 }
        else { d <- if (u[[i]] >= u[[j]]) 1 else -1; u[[i]] <- u[[i]] + d * ox / 2 + d * 1e-3; u[[j]] <- u[[j]] - d * ox / 2 - d * 1e-3 } }
    }
    for (i in seq_len(n)) for (k in seq_len(n)) {   # keep a label's box off every labelled point
      if (abs(u[[i]] - u0[[k]]) < w[[i]] / 2 + .005 && abs(v[[i]] - v0[[k]]) < h / 2 + .005) { moved <- TRUE
        v[[i]] <- v0[[k]] + (if (v[[i]] >= v0[[k]]) 1 else -1) * (h / 2 + .006) }
    }
    u <- pmin(pmax(u, w / 2), 1 - w / 2); v <- pmin(pmax(v, h / 2), 1 - h / 2)
    if (!moved) break
  }
  list(x = xr[[1L]] + u * diff(xr), y = yr[[1L]] + v * diff(yr))
}
tvc_label_annotations <- function(l, xr = range(l$.lx), yr = range(l$.ly)) {
  if (is.null(l) || !nrow(l)) return(list())
  if (!all(is.finite(c(xr, yr))) || diff(xr) <= 0 || diff(yr) <= 0) { xr <- range(l$.lx) + c(-1, 1); yr <- range(l$.ly) + c(-1, 1) }
  at <- tvc_repel_labels(l$.lx, l$.ly, l$tag, xr, yr)
  lapply(seq_len(nrow(l)), function(k) list(x = l$.lx[[k]], y = l$.ly[[k]], ax = at$x[[k]], ay = at$y[[k]], axref = "x", ayref = "y",
    text = l$tag[[k]], xref = "x", yref = "y", showarrow = TRUE, arrowhead = 0, arrowwidth = 1, arrowcolor = "#555555", standoff = 2,
    font = list(size = 12, color = "black"), bgcolor = "rgba(255,255,255,0.8)"))
}
tvc_plotly <- function(p) {
  built <- ytab_ggplotly(p); axis <- function(a) suppressWarnings(as.numeric(built$x$layout[[a]]$range))
  notes <- tvc_label_annotations(attr(p, "tvc_labels"), axis("xaxis"), axis("yaxis"))
  if (length(notes)) built$x$layout$annotations <- c(built$x$layout$annotations %||% list(), notes)
  plotly::config(built, edits = list(annotationPosition = TRUE, annotationTail = TRUE))
}

# The volcano's y axis and the order of top hits: |z| for the z-score, which has no p-values (a gene is called when |z|
# passes the control-control threshold), -log10 p for resampling and the Gumbel.
tvc_significance <- function(d, method) if (identical(method, "zscore")) abs(as.numeric(d$z)) else -log10(pmax(as.numeric(d$p_value), 1e-300))

# The z-score's table per background x treated condition, under the column names of a resampling genes.tsv (renamed,
# not recomputed): mean_abundance is the noise model's abundance A, z_threshold the background's call threshold.
tvc_zscore_genes <- function(output_dir) {
  path <- file.path(output_dir, "tables", "treated_vs_parent.combined.log2fc_z.csv")
  if (!file.exists(path)) return(data.frame())
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  data.frame(group = paste(x$background, x$treated_condition, sep = "__"), background = x$background, condition = x$treated_condition,
             ORF = x$feature_id, Name = if ("common_name" %in% names(x)) x$common_name else NA_character_,
             controls = x$controls, treated = x$treated, mean_control = x$parent_mean, mean_treatment = x$treated_mean,
             mean_abundance = x$A_exp, log2FC = x$log2FC, z = x$z, z_threshold = x$Zthr,
             call = ifelse(x$call %in% c("depleted", "enriched"), x$call, "not significant"),
             status = ifelse(is.finite(x$z), "tested", "not tested"), stringsAsFactors = FALSE)
}

# Rows of `d` to label on the volcano and MA plot: the `n` most significant called genes, and any gene named in `extra`.
tvc_label_rows <- function(d, method, n, extra) cmp_label_rows(d, d$call, tvc_significance(d, method), n, extra)

# Shared by Treated vs Control, Multi-condition and Genetic interactions. Rows of `d` to label: the `n` called genes with
# the largest `sig` (calls not in `uncalled`), and any gene named in `extra` (IDs or names, comma/space separated).
cmp_uncalled <- c("not significant", "other direction", "No interaction", "other class")
cmp_label_rows <- function(d, call, sig, n, extra) {
  n <- suppressWarnings(as.integer(n)); if (!length(n) || is.na(n)) n <- 0L
  called <- which(!as.character(call) %in% cmp_uncalled & is.finite(sig))
  top <- called[order(-sig[called])][seq_len(min(max(0L, n), length(called)))]
  wanted <- toupper(unlist(strsplit(extra %||% "", "[,;[:space:]]+"))); wanted <- wanted[nzchar(wanted)]
  union(top, which(toupper(d$ORF) %in% wanted | toupper(d$Name) %in% wanted))
}
# One plot of a result group, for any of the three methods. s: ytab_plot_customization_values(); extra: sites (the
# gene's sites.tsv.gz rows, or for the z-score its gene_values.long.csv rows), gene, label_n and label_genes (volcano and
# MA), direction ("both", "depleted" or "enriched": the calls shown in colour and labelled), null_z (the z-score's
# control-control z of the group's background), unit (the z-score's value unit), library_labels (sample-sheet names by
# sample id or sites.tsv column), interactive (TRUE: labels are left to tvc_plotly as annotations; FALSE: drawn with
# ggrepel for downloads).
tvc_plot <- function(choice, genes, method, s, extra = list()) {
  theme <- tvc_theme(s); point <- ggplot2::geom_point(alpha = s$point_opacity, size = s$point_size)
  colours <- ggplot2::scale_colour_manual(values = tvc_call_colours, drop = FALSE)
  zscore <- identical(method, "zscore"); unit <- extra$unit %||% ""
  d <- tvc_tested(genes); direction <- extra$direction %||% "both"
  if (direction %in% c("depleted", "enriched")) d$call[d$call %in% setdiff(c("depleted", "enriched"), direction)] <- "other direction"
  d$call <- factor(d$call, levels = tvc_call_levels(method))
  library_label <- function(x) { lab <- extra$library_labels; out <- if (length(lab)) unname(lab[as.character(x)]) else NA_character_
    ifelse(is.na(out), as.character(x), out) }
  name <- function(d) ifelse(is.na(d$Name) | d$Name == d$ORF, "", d$Name)
  hover <- function(d) if (zscore) sprintf("%s %s<br>log2FC %.2f; z %.2f", d$ORF, name(d), d$log2FC, d$z) else
    sprintf("%s %s<br>p %.3g; q %.3g", d$ORF, name(d), d$p_value, d$q_bh)
  # Labelled genes of a plot whose point positions are columns x and y of `d`; attached to the plot as "tvc_labels".
  labelled <- function(p, d, x, y) cmp_label_plot(p, d, tvc_label_rows(d, method, extra$label_n, extra$label_genes), x, y, s, extra$interactive)
  empty <- function(message) ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = message) + ggplot2::theme_void()
  if (identical(choice, "volcano")) {
    x <- if (identical(method, "gumbel")) "log2_scaled_ratio" else "log2FC"
    d$y <- tvc_significance(d, method); d <- d[is.finite(as.numeric(d[[x]])) & is.finite(d$y), , drop = FALSE]; d$hover <- hover(d)
    p <- ggplot2::ggplot(d, ggplot2::aes(.data[[x]], y, colour = call, text = hover)) + point + colours +
      ggplot2::labs(x = if (identical(x, "log2FC")) "log2FC (treated / control)" else "log2 scaled read ratio (treated / control)",
                    y = if (zscore) "|z| (local-abundance z-score)" else "-log10 p", colour = "Call") + theme
    threshold <- if (zscore) unique(d$z_threshold[is.finite(d$z_threshold)]) else numeric()
    if (length(threshold)) p <- p + ggplot2::geom_hline(yintercept = threshold[[1L]], linetype = "dashed", colour = "grey30")
    return(labelled(p, d, x, "y"))
  }
  if (identical(choice, "ma")) {
    d <- d[is.finite(d$log2FC) & is.finite(d$mean_abundance) & d$mean_abundance > 0, , drop = FALSE]
    d$A <- log2(d$mean_abundance); d$hover <- hover(d)
    return(labelled(ggplot2::ggplot(d, ggplot2::aes(A, log2FC, colour = call, text = hover)) + ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
      point + colours +
      ggplot2::labs(x = if (zscore) sprintf("log2 abundance (%s)", if (nzchar(unit)) unit else "per million") else "log2 mean normalized count per site (abundance)",
                    y = "log2FC (treated / control)", colour = "Call") + theme, d, "A", "log2FC"))
  }
  # The MA plot's two means drawn against each other, log10(x + 1) on both axes; the diagonal is no change.
  if (identical(choice, "scatter")) {
    d <- d[is.finite(d$mean_control) & is.finite(d$mean_treatment), , drop = FALSE]
    d$x <- log10(d$mean_control + 1); d$y <- log10(d$mean_treatment + 1); d$hover <- hover(d)
    return(labelled(ggplot2::ggplot(d, ggplot2::aes(x, y, colour = call, text = hover)) + ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40") +
      point + colours + ggplot2::labs(x = "Control (log10 value + 1)", y = "Treated (log10 value + 1)", colour = "Call") + theme, d, "x", "y"))
  }
  if (identical(choice, "pvalues")) {
    d <- d[is.finite(d$p_value), , drop = FALSE]
    return(ggplot2::ggplot(d, ggplot2::aes(p_value)) + ggplot2::geom_histogram(breaks = seq(0, 1, .05), fill = "#7fa7c9") +
      ggplot2::geom_hline(yintercept = nrow(d) / 20, linetype = "dashed", colour = "grey30") +
      ggplot2::labs(x = "p-value", y = "Genes") + theme)
  }
  if (identical(choice, "nullz")) {
    z <- extra$null_z; if (is.null(z) || !nrow(z)) return(empty("No control-control z-scores in this result"))
    z <- z[is.finite(z$z), , drop = FALSE]; threshold <- unique(genes$z_threshold[is.finite(genes$z_threshold)])
    p <- ggplot2::ggplot(z, ggplot2::aes(z)) + ggplot2::geom_histogram(bins = 80, fill = "#7fa7c9") +
      ggplot2::labs(x = "Control-control z", y = "Genes") + theme
    if (length(threshold)) p <- p + ggplot2::geom_vline(xintercept = c(-1, 1) * threshold[[1L]], linetype = "dashed", colour = "grey30")
    return(p)
  }
  if (identical(choice, "runs")) {
    d <- d[is.finite(d$r) & is.finite(d$null_mean_run), , drop = FALSE]
    d$hover <- sprintf("%s<br>observed run %s sites; expected %.1f", hover(d), d$r, d$null_mean_run)
    return(ggplot2::ggplot(d, ggplot2::aes(null_mean_run, r, colour = call, text = hover)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40") +
      ggplot2::geom_jitter(width = .05, height = .15, alpha = s$point_opacity, size = s$point_size) + colours +
      ggplot2::labs(x = "Expected longest depleted run (sites; control model)", y = "Observed longest depleted run (sites)", colour = "Call") + theme)
  }
  if (identical(choice, "sites")) {
    sites <- extra$sites; gene <- extra$gene %||% ""
    shown <- genes$Name[match(gene, genes$ORF)]; if (length(shown) && !is.na(shown) && nzchar(shown)) gene <- shown
    if (is.null(sites) || !nrow(sites)) return(empty(if (zscore) "No values for this gene" else "No sites for this gene"))
    # z-score: a paired dot plot on the scale the test uses. Each library is a point at log2(value + pseudocount), named on
    # the x axis (controls, a gap, then treated libraries, so replicate names never collide); a line
    # marks each condition at log2(mean + pseudocount), since the z-score compares log2 of the condition means; thin
    # lines join each control-treated pair of the comparison design (extra$pairs: parent_sample, treated_sample).
    if (zscore) {
      pc <- as.numeric(extra$pseudocount %||% .5); orf <- extra$gene %||% ""
      sites$role <- factor(sites$role, levels = c("control", "treated")); sites$library <- library_label(sites$sample)
      sites <- sites[order(sites$role, sites$library), , drop = FALSE]
      sites$x <- seq_len(nrow(sites)) + (sites$role == "treated")
      sites$y <- log2(sites$value + pc)
      sites$hover <- sprintf("%s<br>%s (%s)<br>%.2f %s; log2(value + %s) = %.2f", gene, sites$library, sites$role, sites$value, unit, pc, sites$y)
      means <- stats::aggregate(value ~ role, sites, mean); means$y <- log2(means$value + pc)
      means$from <- tapply(sites$x, sites$role, min)[as.character(means$role)] - .35; means$to <- tapply(sites$x, sites$role, max)[as.character(means$role)] + .35
      means$hover <- sprintf("%s mean %.2f %s; log2(mean + %s) = %.2f", means$role, means$value, unit, pc, means$y)
      pairs <- extra$pairs; ia <- match(pairs$parent_sample, sites$sample); ib <- match(pairs$treated_sample, sites$sample); ok <- !is.na(ia) & !is.na(ib)
      pairs <- data.frame(x = sites$x[ia[ok]], y = sites$y[ia[ok]], xend = sites$x[ib[ok]], yend = sites$y[ib[ok]])
      row <- genes[genes$ORF == orf, , drop = FALSE]
      title <- if (nrow(row)) sprintf("%s: log2FC %.2f, z %.2f, %s", gene, row$log2FC[[1L]], row$z[[1L]], row$call[[1L]]) else gene
      colours_role <- c(control = "grey35", treated = "#d6604d")
      p <- ggplot2::ggplot(sites, ggplot2::aes(x, y))
      if (nrow(pairs)) p <- p + ggplot2::geom_segment(data = pairs, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, colour = "grey70", linewidth = .5)
      return(p + ggplot2::geom_segment(data = means, ggplot2::aes(x = from, xend = to, y = y, yend = y, colour = role, text = hover), inherit.aes = FALSE, linewidth = 1.2) +
        ggplot2::geom_point(ggplot2::aes(colour = role, text = hover), size = s$point_size + 2) +
        ggplot2::scale_colour_manual(values = colours_role, labels = c(control = "Control", treated = "Treated"), drop = FALSE) +
        ggplot2::scale_x_continuous(breaks = sites$x, labels = sites$library, expand = ggplot2::expansion(add = .6)) +
        ggplot2::labs(x = NULL, y = sprintf("log2(gene value + %s)", format(pc)), colour = NULL, title = title) +
        theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)))
    }
    libs <- grep("^(control|treatment)_[0-9]+$", names(sites), value = TRUE)
    long <- do.call(rbind, lapply(libs, function(l) data.frame(position = sites$position, library = l, count = as.numeric(sites[[l]]),
      condition = ifelse(startsWith(l, "control"), "control", "treated"), depleted = if ("depleted" %in% names(sites)) sites$depleted %in% c("True", TRUE) else FALSE)))
    long$library <- library_label(long$library)
    long$library[long$library %in% libs] <- sub("^treatment_", "treated ", sub("^control_", "control ", long$library[long$library %in% libs]))
    controls <- sort(unique(long$library[long$condition == "control"])); treated <- sort(unique(long$library[long$condition == "treated"]))
    shades <- c(stats::setNames(grDevices::colorRampPalette(c("grey20", "grey62"))(max(2L, length(controls)))[seq_along(controls)], controls),
                stats::setNames(grDevices::colorRampPalette(c("#b2182b", "#f4a582"))(max(2L, length(treated)))[seq_along(treated)], treated))
    long$hover <- sprintf("%s<br>position %s<br>%s %.1f", gene, format(long$position, big.mark = ","), long$library, long$count)
    p <- ggplot2::ggplot(long, ggplot2::aes(position, count, colour = library, text = hover)) +
      ggplot2::geom_point(size = s$point_size + .5, alpha = max(s$point_opacity, .8)) + ggplot2::scale_colour_manual(values = shades) +
      ggplot2::labs(x = sprintf("%s position (bp)", gene), y = if (identical(method, "gumbel")) "Reads per site" else "Normalized count per site",
                    colour = "Library") + theme
    if (any(long$depleted)) p <- p + ggplot2::geom_rug(data = long[long$depleted & long$condition == "control", , drop = FALSE], sides = "b", colour = YTAB_DIRECTION_COLOURS[["depleted"]], inherit.aes = FALSE, ggplot2::aes(x = position))
    return(p)
  }
  empty("Choose a plot")
}

# zscore_request: a z-score result (run or analysis id) to open, set by the Run subtab's "Open result" buttons.
tvc_server <- function(input, output, session, repo_root, active, active_project_path, norm_runs, start,
                       design_support = reactive(data.frame()), zscore_request = reactive("")) {
  tvc_runs <- reactivePoll(3000, session,
    checkFunc = function() { x <- ytab_analysis_runs(tryCatch(file.path(active()$project_root, "treated_vs_control"), error = function(e) "")); paste(x$path, x$modified, collapse = "|") },
    valueFunc = function() ytab_analysis_runs(tryCatch(file.path(active()$project_root, "treated_vs_control"), error = function(e) "")))
  # z-score results, newest first; the check only lists their result tables (discovery reads them all).
  zscore_tables <- function() { root <- tryCatch(file.path(active()$project_root, "treated_vs_parent"), error = function(e) "")
    x <- if (dir.exists(root)) list.files(root, pattern = "^treated_vs_parent_results\\.csv$", recursive = TRUE, full.names = TRUE) else character()
    paste(x, file.mtime(x), collapse = "|") }
  zscore_runs <- reactivePoll(3000, session, checkFunc = zscore_tables, valueFunc = function() {
    r <- tryCatch(discover_fitness_results(active()), error = function(e) list())
    r[order(-vapply(r, function(x) as.numeric(file.mtime(x$table)), 0))] })
  # Results open on a method the project has results for: the z-score when it has one, else the newest resampling or
  # Gumbel run.
  observeEvent(active()$project_root, {
    zscore <- length(isolate(zscore_runs())) > 0L
    x <- isolate(tvc_runs())
    if (!zscore && nrow(x)) {
      method <- (ytab_read_manifest(x$path[[1L]]) %||% list())$method %||% "resampling"
      updateRadioButtons(session, "tvc_results_method", selected = method)
    }
  })
  design <- reactive({
    path <- file.path(active()$project_root, "config", "comparison_design.csv")
    if (!file.exists(path)) return(data.frame())
    d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    d[tolower(trimws(as.character(d$include))) %in% c("1", "true", "yes", "y"), , drop = FALSE]
  })

  # Z-score input: total reads (CPM) or any method of a normalization run ("<run>|<method>").
  output$fitness_normalization_selector <- renderUI({
    # Default: BGC of the latest run that has it (total reads only when no run has BGC); a valid earlier choice is kept.
    choices <- c(list("Without a normalization run" = c("Total reads (CPM)" = "cpm")), tvc_norm_choices(norm_runs()))
    selected <- if (length(choices) > 1L) tvc_norm_default(choices[-1L], "betageom", isolate(input$fitness_normalization)) else "cpm"
    if (identical(isolate(input$fitness_normalization), "cpm")) selected <- "cpm"
    selectInput("fitness_normalization", "Normalization", choices, selected)
  })
  # Results matches against this choice, so it must exist before the Run subtab is opened.
  outputOptions(output, "fitness_normalization_selector", suspendWhenHidden = FALSE)

  # One normalization choice, grouped by run like the z-score's; defaults follow the pipeline (BGC for resampling, TTR
  # for the Gumbel).
  output$tvc_inputs <- renderUI({
    runs <- norm_runs(); method <- input$tvc_method %||% "resampling"; allowed <- tvc_allowed_methods[[method]] %||% character()
    if (!nrow(runs)) return(tags$p(class = "text-warning", "No normalization run yet. Run Normalization first."))
    choices <- tvc_norm_choices(runs, allowed)
    if (!length(choices)) return(tags$p(class = "text-warning", sprintf("No normalization run has a normalization the %s accepts (%s). Add one in Normalization.",
      if (identical(method, "gumbel")) "Gumbel" else "resampling", paste(as.character(norm_method_label(allowed)), collapse = ", "))))
    selectInput("tvc_norm_choice", "Normalization", choices, tvc_norm_default(choices, if (identical(method, "gumbel")) "ttr" else "betageom", isolate(input$tvc_norm_choice)))
  })
  tvc_norm <- reactive({ parts <- strsplit(input$tvc_norm_choice %||% "", "|", fixed = TRUE)[[1L]]
    if (length(parts) == 2L) list(run = parts[[1L]], method = parts[[2L]]) else NULL })
  output$tvc_run_title <- renderUI(names(tvc_method_choices)[tvc_method_choices == (input$tvc_method %||% "resampling")])
  output$tvc_comparison_selector <- renderUI({
    d <- design()
    if (!nrow(d)) return(tags$p(class = "text-warning", "No comparison design yet. Create it in Design."))
    selectizeInput("tvc_comparisons", "Comparisons", stats::setNames(d$comparison_id, sprintf("%s (%s, %s)", d$comparison_id, d$background, d$treated_condition)),
                   selected = d$comparison_id, multiple = TRUE)
  })
  observeEvent(input$tvc_run, {
    method <- input$tvc_method; comparisons <- input$tvc_comparisons
    if (!length(comparisons)) { showNotification("Select at least one comparison.", type = "error"); return() }
    norm <- tvc_norm()
    if (is.null(norm) || !norm$method %in% tvc_allowed_methods[[method]]) { showNotification("Choose a normalization this method accepts.", type = "error"); return() }
    out <- ytab_new_run_dir(file.path(active()$project_root, "treated_vs_control"), method)
    args <- c("--project-config", active_project_path(), "--method", method, "--normalization-run", norm$run,
              "--normalization", norm$method, "--comparisons", paste(comparisons, collapse = ","),
              "--fdr", as.character(input$tvc_fdr %||% .05), "--output", out)
    start("ytab_treated_vs_control.py", args, "treated_vs_control", sprintf("%s started.", names(tvc_method_choices)[tvc_method_choices == method]))
  }, ignoreInit = TRUE)

  # ---- Results -----------------------------------------------------------------------------------------------------
  # One view for the three methods. Resampling and Gumbel results are run folders with a genes.tsv; a z-score result is
  # a treated_vs_parent folder, read through tvc_zscore_genes().
  results_method <- reactive(input$tvc_results_method %||% "zscore")
  method_runs <- reactive({
    if (identical(results_method(), "zscore")) { r <- zscore_runs()
      return(data.frame(path = vapply(r, `[[`, "", "output_dir"), label = vapply(r, fitness_result_display_name, ""), stringsAsFactors = FALSE)) }
    x <- tvc_runs()
    x <- x[vapply(x$path, function(p) identical(ytab_read_manifest(p)$method, results_method()), logical(1)), , drop = FALSE]
    data.frame(path = x$path, label = names(ytab_run_choices(x)) %||% character(), stringsAsFactors = FALSE)
  })
  # A requested z-score result: selected in the selector if it is on the page, or when the selector next renders (after the
  # method switches to the z-score). It is cleared once the browser reports it selected, not when the selector renders,
  # so a re-render in the same flush cannot fall back to the previous choice.
  pending <- reactiveVal("")
  observeEvent(zscore_request(), { id <- zscore_request() %||% ""; if (!nzchar(id)) return()
    r <- Find(function(x) identical(x$run_id, id) || identical(x$analysis_id, id), zscore_runs())
    if (!is.null(r)) { pending(r$output_dir); updateSelectInput(session, "tvc_result", selected = r$output_dir) } })
  observeEvent(input$tvc_result, if (identical(input$tvc_result, pending())) pending(""))
  output$tvc_result_selector <- renderUI({
    x <- method_runs()
    if (!nrow(x)) return(tags$p(class = "text-muted", "No result for this method yet; run it from Treated vs Control: Run."))
    want <- isolate(pending()); kept <- isolate(input$tvc_result)
    choice <- if (want %in% x$path) want else if (!is.null(kept) && kept %in% x$path) kept else x$path[[1L]]
    selectInput("tvc_result", "Result", stats::setNames(x$path, x$label), choice)
  })
  selected <- reactive({ x <- method_runs(); req(nrow(x) > 0); path <- input$tvc_result %||% x$path[[1L]]; req(path %in% x$path); path })
  data <- reactive({
    run <- selected(); method <- results_method()
    if (identical(method, "zscore")) {
      r <- Find(function(x) identical(x$output_dir, run), zscore_runs()); req(r)
      manifest <- list(method = "zscore", method_label = names(tvc_method_choices)[tvc_method_choices == "zscore"], current = isTRUE(r$current),
                       normalization = r$normalization, value_unit = r$value_unit, pseudocount = r$settings$pseudocount %||% .5)
      genes <- tvc_zscore_genes(run)
    } else { manifest <- ytab_read_manifest(run); genes <- ytab_read_analysis_tsv(file.path(run, "genes.tsv")) }
    if (nrow(genes) && "ORF" %in% names(genes)) genes$Name <- ytab_display_gene_names(genes$ORF, genes$Name, repo_root)
    list(method = method, manifest = manifest, genes = genes)
  })
  output$tvc_group_selector <- renderUI({
    g <- data()$genes; req(nrow(g) > 0); groups <- unique(g[, c("group", "background", "condition")])
    if (nrow(groups) < 2L) return(NULL)
    selectInput("tvc_group", "Background and condition", stats::setNames(groups$group, paste(groups$background, groups$condition, sep = ", ")))
  })
  group_genes <- reactive({ g <- data()$genes; req(nrow(g) > 0); group <- input$tvc_group %||% g$group[[1L]]
    if (!group %in% g$group) group <- g$group[[1L]]; g[g$group == group, , drop = FALSE] })
  output$tvc_summary <- renderUI({
    gate <- ytab_design_gate(design_support(), "Treated vs control", "Treated versus control is not available for this dataset")
    if (!is.null(gate)) return(gate)
    if (!nrow(method_runs())) return(tagList(ytab_design_merge_warning(design_support()),
      ytab_empty_state("No result yet", "Choose a normalization and comparisons in Treated vs Control: Run.")))
    d <- data(); m <- d$manifest; req(m)
    if (identical(d$method, "zscore")) {
      g <- d$genes; groups <- unique(g[, c("group", "background", "condition")])
      cards <- lapply(seq_len(nrow(groups)), function(i) { x <- g[g$group == groups$group[[i]], , drop = FALSE]
        ytab_metric_card(paste(groups$background[[i]], groups$condition[[i]], sep = ", "), sprintf("%s ↓ / %s ↑", sum(x$call == "depleted"), sum(x$call == "enriched")), "depleted / enriched") })
      norm <- if (toupper(m$normalization %||% "CPM") == "CPM") "total reads (CPM)" else as.character(norm_method_label(m$normalization))
      return(tagList(tags$div(class = "ytab-stat-grid", cards),
        tags$p(class = "text-muted", sprintf("%s; normalization %s; a gene is called when |z| passes the 99th percentile of control-control |z|.", m$method_label, norm)),
        if (!isTRUE(m$current)) ytab_warning_banner("Earlier pipeline version", "This result was made before the current output contract; rerun it for every plot.")))
    }
    cards <- unlist(lapply(names(m$groups), function(k) { g <- m$groups[[k]]
      if (!identical(g$status, "success")) return(list(ytab_metric_card(paste(g$background, g$condition, sep = ", "), "failed", g$message %||% "")))
      calls <- g$calls
      value <- if (identical(m$method, "gumbel")) format(calls[["depleted run"]] %||% 0, big.mark = ",") else sprintf("%s ↓ / %s ↑", calls$depleted %||% 0, calls$enriched %||% 0)
      list(ytab_metric_card(paste(g$background, g$condition, sep = ", "), value, if (identical(m$method, "gumbel")) "depleted runs" else "depleted / enriched")) }), recursive = FALSE)
    warnings <- Filter(Negate(is.null), lapply(m$groups, `[[`, "warning"))
    tagList(tags$div(class = "ytab-stat-grid", cards),
      tags$p(class = "text-muted", sprintf("%s; normalization %s; FDR %s.", m$method_label, as.character(norm_method_label(m$normalization)), m$fdr)),
      if (length(warnings)) ytab_warning_banner("Unequal replicate noise", paste(unique(unlist(warnings)), collapse = " ")))
  })
  output$tvc_plot_selector <- renderUI({ choices <- tvc_plot_choices(data()$method); kept <- isolate(input$tvc_plot)
    selectInput("tvc_plot", "Plot", choices, if (!is.null(kept) && kept %in% choices) kept else "volcano") })
  output$tvc_gene_selector <- renderUI({
    req(identical(input$tvc_plot, "sites")); g <- tvc_tested(group_genes()); req(nrow(g) > 0)
    g <- g[order(-tvc_significance(g, data()$method)), , drop = FALSE]
    # Thousands of genes: the choices are searched on the server, filled in once the empty selector is on the page.
    choices <- stats::setNames(g$ORF, ifelse(is.na(g$Name) | g$Name == g$ORF, g$ORF, paste0(g$ORF, " (", g$Name, ")")))
    session$onFlushed(function() updateSelectizeInput(session, "tvc_gene", choices = choices, selected = g$ORF[[1L]], server = TRUE), once = TRUE)
    selectizeInput("tvc_gene", "Gene", choices = NULL, options = list(placeholder = "Type a gene name or ID"))
  })
  # Gene view rows: the group's sites.tsv.gz (resampling, Gumbel) or, for the z-score, the gene values of the group's
  # control and treated libraries, marked by role.
  group_sites <- reactive({
    req(identical(input$tvc_plot, "sites")); g <- group_genes()
    if (identical(data()$method, "zscore")) {
      path <- file.path(selected(), "tables", "gene_values.long.csv"); if (!file.exists(path)) return(data.frame())
      v <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
      controls <- strsplit(g$controls[[1L]], ";", fixed = TRUE)[[1L]]; treated <- strsplit(g$treated[[1L]], ";", fixed = TRUE)[[1L]]
      v <- v[v$sample %in% c(controls, treated), , drop = FALSE]; v$role <- ifelse(v$sample %in% controls, "control", "treated")
      return(v)
    }
    path <- file.path(selected(), g$group[[1L]], "sites.tsv.gz")
    if (!file.exists(path)) return(data.frame())
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  })
  gene_sites <- reactive({ s <- group_sites(); gene <- input$tvc_gene %||% ""; if (!nrow(s) || !nzchar(gene)) return(s[0, , drop = FALSE])
    if (identical(data()$method, "zscore")) return(s[s$feature_id == gene, , drop = FALSE])
    s[vapply(strsplit(as.character(s$genes), "[;,]"), function(x) gene %in% x, logical(1)), , drop = FALSE] })
  # The z-score's control-treated pairs (comparison design) of the group's background, for the gene-values plot.
  zscore_pairs <- reactive({
    req(identical(input$tvc_plot, "sites"), identical(data()$method, "zscore"))
    path <- file.path(selected(), "tables", "treated_vs_parent.by_pool.log2fc_z.csv"); empty <- data.frame(parent_sample = character(), treated_sample = character())
    if (!file.exists(path)) return(empty)
    x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    unique(x[x$background == group_genes()$background[[1L]], c("parent_sample", "treated_sample")])
  })
  # The z-score's control-control z of the group's background: the null its call threshold is taken from.
  null_z <- reactive({
    req(identical(input$tvc_plot, "nullz")); path <- file.path(selected(), "tables", "control_control_z_scores.csv")
    if (!file.exists(path)) return(data.frame())
    z <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE); z[z$background == group_genes()$background[[1L]], , drop = FALSE]
  })
  output$tvc_plot_caption <- renderUI({ zscore <- identical(data()$method, "zscore"); gumbel <- identical(data()$method, "gumbel")
    switch(input$tvc_plot %||% "volcano",
    ma = if (zscore) "Each point is a gene: its abundance in the noise model against log2 fold change; colours are the calls." else
      "Each point is a gene: mean normalized count per site (abundance) against log2 fold change; colours are calls at the chosen FDR.",
    scatter = sprintf("Each gene's mean in the treated libraries against its mean in the controls, as log10(value + 1) of %s; genes off the diagonal changed, colours are the calls.",
      if (zscore) sprintf("the gene value (%s)", data()$manifest$value_unit %||% "per million") else "the normalized count per site"),
    pvalues = "p-values of the tested genes. With no effect they are uniform (dashed line: expected count per bin); a peak near 0 holds the affected genes.",
    nullz = "z-scores between the two control libraries, where no gene should differ: the null the call threshold comes from (dashed lines: ±threshold, the 99th percentile of |z|).",
    runs = "Each gene's longest run of depleted sites against the longest run expected by chance from the controls' count model; genes far above the line are depleted runs.",
    sites = if (zscore) sprintf("Each point is one library's gene value (normalized site counts summed over the gene, in %s) on the z-score's log2 scale; the bars are each condition's mean, which the z-score compares, and grey lines join the control-treated pairs of the comparison design.", data()$manifest$value_unit %||% "reads per million") else if (gumbel)
      "Reads at each insertion site of the gene, per library; red ticks mark sites the test counts as depleted." else
      "Normalized count at each insertion site of the gene, per library (before recentering).",
    if (zscore) "Each point is a gene: log2 fold change against |z|; genes above the dashed line (the control-control threshold) are called." else if (gumbel)
      "Each point is a gene tested for a depleted run of sites; the x axis is the gene's treated/control read ratio after scaling." else
      "Each point is a gene; colours are calls at the chosen FDR.") })
  display <- reactive(ytab_plot_customization_values(input, "tvc"))
  output$tvc_plot_frame <- renderUI({ req(nrow(group_genes()) > 0); s <- display()
    ytab_plot_frame(plotly::plotlyOutput("tvc_plotly", height = s$height), s$width) })
  # Sample-sheet names of the group's libraries, by sample id and (resampling, Gumbel) by sites.tsv column: the manifest
  # lists the controls and treatments in column order.
  library_labels <- reactive({
    d <- data(); g <- group_genes(); root <- active()$project_root
    if (identical(d$method, "zscore")) { ids <- unique(unlist(strsplit(c(g$controls[[1L]], g$treated[[1L]]), ";", fixed = TRUE)))
      return(stats::setNames(ytab_sample_labels(ids, root), ids)) }
    m <- d$manifest$groups[[g$group[[1L]]]]; controls <- unlist(m$controls); treated <- unlist(m$treatments)
    stats::setNames(ytab_sample_labels(c(controls, treated), root),
                    c(sprintf("control_%d", seq_along(controls)), sprintf("treatment_%d", seq_along(treated))))
  })
  # interactive: labels as draggable annotations with leader lines (tvc_plotly); downloads draw them with ggrepel.
  build_plot <- function(interactive) { choice <- input$tvc_plot %||% "volcano"; d <- data()
    tvc_plot(choice, group_genes(), d$method, display(), list(
      sites = if (identical(choice, "sites")) gene_sites(), gene = input$tvc_gene, null_z = if (identical(choice, "nullz")) null_z(),
      label_n = input$tvc_label_n %||% 10, label_genes = input$tvc_label_genes, unit = d$manifest$value_unit %||% "", pseudocount = d$manifest$pseudocount,
      pairs = if (identical(choice, "sites") && identical(d$method, "zscore")) zscore_pairs(),
      direction = if (identical(d$method, "gumbel")) "both" else input$tvc_direction %||% "both", library_labels = library_labels(),
      interactive = interactive)) }
  current_plot <- reactive(build_plot(TRUE))
  output$tvc_plotly <- plotly::renderPlotly(tvc_plotly(current_plot()))
  output$tvc_table <- DT::renderDT({
    g <- tvc_tested(group_genes()); g <- g[order(-tvc_significance(g, data()$method)), , drop = FALSE]
    keep <- intersect(c("ORF", "Name", "n_sites", "n_eligible", "n_depleted", "run_span_bp", "mean_control", "mean_treatment",
                        "log2FC", "log2_scaled_ratio", "z", "p_value", "q_bh", "call"), names(g))
    g <- g[, keep, drop = FALSE]
    for (column in intersect(c("mean_control", "mean_treatment", "log2FC", "log2_scaled_ratio", "z"), names(g))) g[[column]] <- signif(g[[column]], 3)
    for (column in intersect(c("p_value", "q_bh"), names(g))) g[[column]] <- formatC(g[[column]], format = "g", digits = 3)
    ytab_datatable(g, page_length = 15L)
  })
  name <- function(ext) paste0(project_display_label(active()), "_treated_vs_control_", basename(selected()), ext)
  for (format in c("png", "pdf")) local({ ext <- format
    output[[paste0("tvc_download_", ext)]] <- downloadHandler(filename = function() name(paste0(".", ext)), content = function(file) { s <- display()
      ytab_save_ggplot(file, build_plot(FALSE), width = switch(s$width, wide = 13, full = 16, 10), height = max(4, as.numeric(sub("px", "", s$height)) / 75)) }) })
  output$tvc_download_genes <- downloadHandler(filename = function() name("_genes.tsv"), content = function(file) {
    if (identical(data()$method, "zscore")) utils::write.table(data()$genes, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    else file.copy(file.path(selected(), "genes.tsv"), file) })
  tvc_agreement_server(input, output, session, repo_root, active, start, zscore_runs, tvc_runs)
}

# ---- Compare methods ---------------------------------------------------------------------------------------------------
# How the three methods' calls agree on the same comparisons: `ytab tvc-agreement` (ytab.fitness.agreement) writes the
# per-gene table, the pairwise numbers and the call intersections into <project>/treated_vs_control/agreement/<run>;
# the app only reads them.
tvc_agreement_methods <- c(zscore = "z-score", resampling = "Resampling", gumbel = "Gumbel", published = "Published")
tvc_agreement_colours <- c(zscore = "#0072B2", resampling = "#E69F00", gumbel = "#009E73", published = "#CC79A7")

tvc_agreement_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Compare methods",
    uiOutput("tvc_ag_inputs"),
    actionButton("tvc_ag_run", "Compare", class = "btn-primary"),
    tags$hr(),
    uiOutput("tvc_ag_result_selector"),
    uiOutput("tvc_ag_group_selector"),
    selectInput("tvc_ag_plot", "Plot", c("Recovery matrix" = "matrix", "Call overlap (Venn)" = "venn", "Call overlap (UpSet)" = "sets",
                                         "Effect scatter" = "effect", "Significance scatter" = "score")),
    conditionalPanel("input.tvc_ag_plot == 'venn'", uiOutput("tvc_ag_venn_selector")),
    conditionalPanel("input.tvc_ag_plot == 'matrix' || input.tvc_ag_plot == 'sets' || input.tvc_ag_plot == 'venn'",
      selectInput("tvc_ag_direction", "Calls", c("Depleted" = "depleted", "Enriched (the Gumbel is depletion only)" = "enriched"))),
    conditionalPanel("input.tvc_ag_plot == 'effect' || input.tvc_ag_plot == 'score'", uiOutput("tvc_ag_pair_selector")),
    ytab_plot_downloads_ui("tvcag", c(genes = "Download per-gene table (TSV)", pairs = "Download agreement (TSV)")),
    ytab_plot_customization_controls("tvcag", include_points = TRUE, include_labels = FALSE, open = FALSE)),
  main = tagList(uiOutput("tvc_ag_summary"),
    ytab_plot_card("Method agreement", uiOutput("tvc_ag_plot_frame"), description = uiOutput("tvc_ag_caption")),
    ytab_table_card("Agreement", DT::DTOutput("tvc_ag_table")),
    ytab_methods_drawer(
      tags$p("Each pair of methods is compared on the genes both tested, per background and treated condition. ",
             tags$b("Calls"), " are compared by recovery: of one method's calls in a direction, the share the other method also calls (with published calls as the first method, this is the recall of the published hits), shown with how many genes each calls. ",
             tags$b("Effects"), " and ", tags$b("signed significance"), " (z for the z-score; sign(effect) × −log10 p for resampling; −(−log10 p) for the depletion-only Gumbel) are compared by Spearman rank correlation, since the methods' scales and tails differ. The Gumbel tests depletion only and is left out of enriched comparisons."),
      tags$p("Command: ", tags$code("ytab tvc-agreement --zscore RUN --resampling RUN --gumbel RUN --output DIR"), "."))))

# Recovery matrix: tile (row, column) = of the row method's calls in `direction`, the share the column method also calls
# (the published row is the recall of the published hits); the diagonal gives each method's number of calls.
tvc_agreement_matrix <- function(pairs, direction, s) {
  n <- function(k) paste0(direction, "_", k)
  both <- rbind(data.frame(row = pairs$method_a, col = pairs$method_b, calls = pairs[[n("n_a")]], other = pairs[[n("n_b")]], shared = pairs[[n("n_both")]],
                           value = pairs[[n("a_recovered_by_b")]], rho = pairs$spearman_score),
                data.frame(row = pairs$method_b, col = pairs$method_a, calls = pairs[[n("n_b")]], other = pairs[[n("n_a")]], shared = pairs[[n("n_both")]],
                           value = pairs[[n("b_recovered_by_a")]], rho = pairs$spearman_score))
  m <- unique(c(pairs$method_a, pairs$method_b)); lab <- tvc_agreement_methods[m]
  both$r <- factor(tvc_agreement_methods[both$row], rev(lab)); both$c <- factor(tvc_agreement_methods[both$col], lab)
  both$label <- ifelse(is.na(both$value), "depletion\nonly", sprintf("%.0f%%\n%s of %s", 100 * both$value, both$shared, both$calls))
  both$tip <- sprintf("Of %s's %s %s calls, %s (%.0f%%) are also called by %s<br>Spearman of signed significance %.2f", tvc_agreement_methods[both$row], both$calls, direction,
                      both$shared, 100 * both$value, tvc_agreement_methods[both$col], both$rho)
  own <- do.call(rbind, lapply(m, function(k) { x <- both$calls[both$row == k]; x <- suppressWarnings(max(x, na.rm = TRUE))
    data.frame(r = factor(tvc_agreement_methods[[k]], levels(both$r)), c = factor(tvc_agreement_methods[[k]], levels(both$c)),
               label = if (is.finite(x)) sprintf("%s calls", format(x, big.mark = ",")) else "depletion\nonly", tip = tvc_agreement_methods[[k]]) }))
  ggplot2::ggplot(both, ggplot2::aes(c, r, fill = value, text = tip)) + ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_tile(data = own, ggplot2::aes(c, r, text = tip), fill = "grey92", colour = "white", linewidth = 1, inherit.aes = FALSE) +
    ggplot2::geom_text(data = rbind(both[, c("r", "c", "label")], own[, c("r", "c", "label")]), ggplot2::aes(c, r, label = label), size = s$text_size / 3.6,
                       fontface = "bold", inherit.aes = FALSE) +
    ggplot2::scale_fill_gradient(low = "white", high = YTAB_FOCAL, limits = c(0, 1), na.value = "grey92", labels = function(x) paste0(100 * x, "%"),
                                 name = sprintf("Row's %s calls\nalso called by column", direction)) +
    ggplot2::coord_equal() + ggplot2::labs(x = "Also called by", y = "Calls of") + tvc_theme(s) + ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# Venn of the calls in one direction for two or three methods: region counts are the genes called by exactly that set of
# the chosen methods (from sets.tsv, genes tested by every method of the comparison).
tvc_agreement_venn <- function(sets, methods, s) {
  validate(need(length(methods) %in% 2:3, "Choose two or three methods for the Venn diagram."))
  member <- lapply(strsplit(sets$methods, "&", fixed = TRUE), function(x) intersect(x, methods))
  key <- vapply(member, function(x) paste(sort(match(x, methods)), collapse = ""), "")
  counts <- tapply(sets$genes, key, sum); get <- function(k) if (k %in% names(counts)) counts[[k]] else 0L
  circle <- function(x, y, r = 1.15, k) { t <- seq(0, 2 * pi, length.out = 120); data.frame(x = x + r * cos(t), y = y + r * sin(t), method = k) }
  centres <- if (length(methods) == 2L) list(c(-.6, 0), c(.6, 0)) else list(c(-.65, .45), c(.65, .45), c(0, -.6))
  circles <- do.call(rbind, lapply(seq_along(methods), function(i) circle(centres[[i]][1], centres[[i]][2], k = methods[[i]])))
  regions <- if (length(methods) == 2L) data.frame(k = c("1", "2", "12"), x = c(-1.15, 1.15, 0), y = 0) else
    data.frame(k = c("1", "2", "3", "12", "13", "23", "123"), x = c(-1.2, 1.2, 0, 0, -.75, .75, 0), y = c(.85, .85, -1.25, .95, -.25, -.25, .1))
  regions$label <- format(vapply(regions$k, get, 0), big.mark = ",")
  names_at <- do.call(rbind, lapply(seq_along(methods), function(i) data.frame(x = centres[[i]][1] * 1.9, y = centres[[i]][2] * 1.9 + if (centres[[i]][2] < 0) -.55 else .5,
                                                                             label = tvc_agreement_methods[[methods[[i]]]])))
  ggplot2::ggplot() + ggplot2::geom_polygon(data = circles, ggplot2::aes(x, y, group = method, fill = method), alpha = .25, colour = "black", linewidth = .6) +
    ggplot2::geom_text(data = regions, ggplot2::aes(x, y, label = label), size = s$text_size / 2.6, fontface = "bold") +
    ggplot2::geom_text(data = names_at, ggplot2::aes(x, y, label = label), size = s$text_size / 3, fontface = "bold") +
    ggplot2::scale_fill_manual(values = tvc_agreement_colours, guide = "none") + ggplot2::coord_equal() + ggplot2::theme_void()
}

# Effect or significance of one method against another, genes tested by both; colour: which of the two called the gene.
tvc_agreement_scatter <- function(genes, a, b, what, s) {
  x <- paste0(a, "_", what); y <- paste0(b, "_", what)
  d <- genes[genes[[paste0(a, "_tested")]] %in% c(TRUE, "True") & genes[[paste0(b, "_tested")]] %in% c(TRUE, "True") & is.finite(genes[[x]]) & is.finite(genes[[y]]), , drop = FALSE]
  # with a depletion-only method in the pair, only depleted calls count (enrichment cannot be called by both)
  called <- if (any(c(a, b) == "gumbel")) "depleted" else c("depleted", "enriched")
  ca <- d[[paste0(a, "_call")]] %in% called; cb <- d[[paste0(b, "_call")]] %in% called
  na <- tvc_agreement_methods[[a]]; nb <- tvc_agreement_methods[[b]]
  d$agree <- factor(ifelse(ca & cb, "called by both", ifelse(ca, paste(na, "only"), ifelse(cb, paste(nb, "only"), "neither"))),
                    c("called by both", paste(na, "only"), paste(nb, "only"), "neither"))
  rho <- suppressWarnings(stats::cor(d[[x]], d[[y]], method = "spearman"))
  d$tip <- sprintf("%s %s<br>%s %.2f (%s)<br>%s %.2f (%s)", d$ORF, ifelse(is.na(d$Name) | d$Name == d$ORF, "", d$Name), na, d[[x]], d[[paste0(a, "_call")]], nb, d[[y]], d[[paste0(b, "_call")]])
  d <- d[order(d$agree != "neither"), , drop = FALSE]
  label <- function(m) if (what == "effect") (if (m == "gumbel") "log2 scaled read ratio" else "log2FC") else (if (m == "zscore") "z" else if (m == "gumbel") "−(−log10 p)" else "signed −log10 p")
  ggplot2::ggplot(d, ggplot2::aes(.data[[x]], .data[[y]], colour = agree, text = tip)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80") + ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
    ggplot2::geom_point(alpha = s$point_opacity, size = s$point_size) +
    ggplot2::scale_colour_manual(values = stats::setNames(c("#222222", tvc_agreement_colours[[a]], tvc_agreement_colours[[b]], "grey82"), levels(d$agree)), drop = FALSE) +
    ggplot2::labs(x = sprintf("%s: %s", na, label(a)), y = sprintf("%s: %s", nb, label(b)), colour = NULL,
                  title = sprintf("Spearman ρ = %.3f (%s genes)", rho, format(nrow(d), big.mark = ","))) + tvc_theme(s)
}

# UpSet of the calls in one direction: bars are genes in exactly the marked set of methods; dots mark the methods.
tvc_agreement_upset <- function(sets, methods, s) {
  validate(need(nrow(sets) > 0, "No calls to compare in this direction."))
  d <- sets[order(-sets$genes), , drop = FALSE]; d$set <- factor(seq_len(nrow(d)))
  members <- strsplit(d$methods, "&", fixed = TRUE)
  d$group <- ifelse(d$n_methods == length(methods), "all", ifelse(d$n_methods == 1, vapply(members, `[`, "", 1L), "some"))
  fills <- c(all = YTAB_FOCAL, some = YTAB_COMPARATOR, tvc_agreement_colours[methods])
  ymax <- max(d$genes); step <- ymax * .1; rows <- -step * seq_along(methods)
  dots <- do.call(rbind, lapply(seq_along(methods), function(i) data.frame(set = d$set, y = rows[[i]], on = vapply(members, function(m) methods[[i]] %in% m, TRUE))))
  links <- do.call(rbind, lapply(seq_len(nrow(d)), function(k) { on <- rows[methods %in% members[[k]]]; if (length(on) > 1) data.frame(set = d$set[[k]], ymin = min(on), ymax = max(on)) }))
  d$tip <- sprintf("Called by: %s<br>%s genes", gsub("&", ", ", d$methods), format(d$genes, big.mark = ","))
  ticks <- pretty(c(0, ymax)); ticks <- ticks[ticks <= ymax * 1.05]
  p <- ggplot2::ggplot(d, ggplot2::aes(set, genes)) + ggplot2::geom_col(ggplot2::aes(fill = group, text = tip), width = .7) +
    ggplot2::geom_text(ggplot2::aes(y = genes + ymax * .045, label = format(genes, big.mark = ",")), size = s$text_size / 3.8, fontface = "bold") +
    ggplot2::geom_point(data = dots, ggplot2::aes(set, y, colour = on), size = 3.2, inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey85"), guide = "none") +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::scale_y_continuous(breaks = c(ticks, rows), labels = c(format(ticks, big.mark = ","), tvc_agreement_methods[methods]), expand = ggplot2::expansion(mult = c(.02, .1))) +
    ggplot2::scale_x_discrete(labels = NULL) + ggplot2::labs(x = "Methods that call the genes (exactly this set)", y = "Genes") +
    tvc_theme(s) + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  if (!is.null(links) && nrow(links)) p <- p + ggplot2::geom_segment(data = links, ggplot2::aes(x = set, xend = set, y = ymin, yend = ymax), linewidth = .8, inherit.aes = FALSE)
  p
}

tvc_agreement_server <- function(input, output, session, repo_root, active, start, zscore_runs, tvc_runs) {
  root <- reactive(file.path(active()$project_root, "treated_vs_control", "agreement"))
  runs <- reactivePoll(3000, session, checkFunc = function() { x <- ytab_analysis_runs(tryCatch(root(), error = function(e) "")); paste(x$path, x$modified, collapse = "|") },
                       valueFunc = function() ytab_analysis_runs(tryCatch(root(), error = function(e) "")))
  output$tvc_ag_inputs <- renderUI({
    z <- zscore_runs(); x <- tvc_runs(); kind <- vapply(x$path, function(p) ytab_read_manifest(p)$method %||% "", "")
    none <- c("Do not include" = "")
    zc <- c(stats::setNames(vapply(z, `[[`, "", "output_dir"), vapply(z, fitness_result_display_name, "")), none)
    pick <- function(m) { y <- x[kind == m, , drop = FALSE]; c(ytab_run_choices(y), none) }
    tagList(selectInput("tvc_ag_zscore", "z-score result", zc), selectInput("tvc_ag_resampling", "Resampling result", pick("resampling")),
            selectInput("tvc_ag_gumbel", "Gumbel result", pick("gumbel")),
            helpText("The latest result of each method is chosen; pick results made on the same comparisons."))
  })
  observeEvent(input$tvc_ag_run, {
    chosen <- c(zscore = input$tvc_ag_zscore %||% "", resampling = input$tvc_ag_resampling %||% "", gumbel = input$tvc_ag_gumbel %||% "")
    chosen <- chosen[nzchar(chosen)]
    if (length(chosen) < 2L) { showNotification("Choose results of at least two methods.", type = "error"); return() }
    out <- ytab_new_run_dir(root(), paste(names(chosen), collapse = "_"))
    start("ytab_tvc_agreement.py", c(as.vector(rbind(paste0("--", names(chosen)), unname(chosen))), "--output", out), "tvc_agreement", "Method comparison started.")
  }, ignoreInit = TRUE)
  output$tvc_ag_result_selector <- renderUI({ x <- runs()
    if (!nrow(x)) return(tags$p(class = "text-muted", "No comparison yet: choose results above and click Compare."))
    selectInput("tvc_ag_result", "Comparison", ytab_run_choices(x)) })
  selected <- reactive({ x <- runs(); req(nrow(x) > 0); p <- input$tvc_ag_result %||% x$path[[1L]]; req(p %in% x$path); p })
  data <- reactive({ run <- selected(); genes <- ytab_read_analysis_tsv(file.path(run, "genes.tsv"))
    if (nrow(genes)) genes$Name <- ytab_display_gene_names(genes$ORF, genes$Name, repo_root)
    list(manifest = ytab_read_manifest(run), genes = genes, pairs = ytab_read_analysis_tsv(file.path(run, "pairs.tsv")), sets = ytab_read_analysis_tsv(file.path(run, "sets.tsv"))) })
  output$tvc_ag_group_selector <- renderUI({ p <- data()$pairs; req(nrow(p) > 0); g <- unique(p[, c("group", "background", "condition")])
    if (nrow(g) < 2L) return(NULL)
    selectizeInput("tvc_ag_group", "Background and condition", stats::setNames(g$group, paste(g$background, g$condition, sep = ", "))) })
  group <- reactive({ p <- data()$pairs; req(nrow(p) > 0); g <- input$tvc_ag_group %||% p$group[[1L]]; if (!g %in% p$group) p$group[[1L]] else g })
  output$tvc_ag_pair_selector <- renderUI({ m <- unlist(data()$manifest$methods); req(length(m) >= 2)
    pairs <- utils::combn(m, 2, simplify = FALSE)
    selectInput("tvc_ag_pair", "Methods", stats::setNames(vapply(pairs, paste, "", collapse = "|"), vapply(pairs, function(x) paste(tvc_agreement_methods[x], collapse = " vs "), ""))) })
  output$tvc_ag_venn_selector <- renderUI({ m <- unlist(data()$manifest$methods); req(length(m) >= 2)
    selectizeInput("tvc_ag_venn", "Methods (two or three)", stats::setNames(m, tvc_agreement_methods[m]), selected = utils::head(m, 3), multiple = TRUE, options = list(maxItems = 3)) })
  output$tvc_ag_summary <- renderUI({
    if (!nrow(runs())) return(ytab_empty_state("No method comparison yet", "Choose a result of each method and click Compare; the Gumbel is optional."))
    p <- data()$pairs; p <- p[p$group == group(), , drop = FALSE]; req(nrow(p) > 0)
    cards <- lapply(seq_len(nrow(p)), function(i) ytab_metric_card(paste(tvc_agreement_methods[[p$method_a[[i]]]], "vs", tvc_agreement_methods[[p$method_b[[i]]]]),
      sprintf("%s shared", p$depleted_n_both[[i]]), sprintf("depleted: %s and %s calls; ρ %.2f (significance)", p$depleted_n_a[[i]], p$depleted_n_b[[i]], p$spearman_score[[i]])))
    tags$div(class = "ytab-stat-grid", cards)
  })
  display <- reactive(ytab_plot_customization_values(input, "tvcag"))
  plot_obj <- reactive({
    d <- data(); s <- display(); g <- group(); choice <- input$tvc_ag_plot %||% "matrix"; dir <- input$tvc_ag_direction %||% "depleted"
    if (choice == "matrix") return(tvc_agreement_matrix(d$pairs[d$pairs$group == g, , drop = FALSE], dir, s))
    if (choice == "venn") { m <- input$tvc_ag_venn %||% utils::head(unlist(d$manifest$methods), 3); if (dir == "enriched") m <- setdiff(m, unlist(d$manifest$depletion_only))
      return(tvc_agreement_venn(d$sets[d$sets$group == g & d$sets$direction == dir, , drop = FALSE], m, s)) }
    if (choice == "sets") { m <- unlist(d$manifest$methods); if (dir == "enriched") m <- setdiff(m, unlist(d$manifest$depletion_only))
      return(tvc_agreement_upset(d$sets[d$sets$group == g & d$sets$direction == dir, , drop = FALSE], m, s)) }
    pair <- strsplit(input$tvc_ag_pair %||% paste(unlist(d$manifest$methods)[1:2], collapse = "|"), "|", fixed = TRUE)[[1L]]
    tvc_agreement_scatter(d$genes[d$genes$group == g, , drop = FALSE], pair[[1L]], pair[[2L]], choice, s)
  })
  output$tvc_ag_plot_frame <- renderUI({ req(nrow(runs()) > 0); s <- display(); ytab_plot_frame(plotly::plotlyOutput("tvc_ag_plotly", height = s$height), s$width) })
  output$tvc_ag_plotly <- plotly::renderPlotly({ built <- ytab_ggplotly(plot_obj())   # a Venn has no axes
    if (identical(input$tvc_ag_plot, "venn")) built <- plotly::layout(built, xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)); built })
  output$tvc_ag_caption <- renderUI(switch(input$tvc_ag_plot %||% "matrix",
    matrix = "Of the row method's calls in the chosen direction, the share the column method also calls (with the count); the published row is the recall of the published hits. The diagonal gives each method's number of calls.",
    venn = "Genes called in the chosen direction by each combination of the chosen methods; only genes tested by every method of the comparison. The Gumbel is depletion only.",
    effect = "Each gene's effect in one method against the other, over genes both tested; colours show which method called it (either direction; depleted only when the Gumbel is one of the pair).",
    score = "Each gene's signed significance in one method against the other (z; sign(effect) × −log10 p; the Gumbel's depletion as −(−log10 p)).",
    sets = "Genes called in the chosen direction by exactly the marked set of methods, largest sets first; only genes tested by every method shown."))
  output$tvc_ag_table <- DT::renderDT({ p <- data()$pairs; p <- p[p$group == group(), , drop = FALSE]; req(nrow(p) > 0)
    out <- data.frame(Methods = paste(tvc_agreement_methods[p$method_a], "vs", tvc_agreement_methods[p$method_b]), `Genes tested by both` = p$genes_tested_by_both,
      `Spearman (effect)` = round(p$spearman_effect, 3), `Spearman (significance)` = round(p$spearman_score, 3),
      `Depleted: A / B / shared` = sprintf("%s / %s / %s", p$depleted_n_a, p$depleted_n_b, p$depleted_n_both),
      `Depleted: A recovered by B` = sprintf("%.0f%%", 100 * p$depleted_a_recovered_by_b), `Depleted: B recovered by A` = sprintf("%.0f%%", 100 * p$depleted_b_recovered_by_a),
      `Enriched: A / B / shared` = ifelse(is.na(p$enriched_n_a), "depletion only", sprintf("%s / %s / %s", p$enriched_n_a, p$enriched_n_b, p$enriched_n_both)), check.names = FALSE)
    ytab_datatable(out, page_length = 10L) })
  name <- function(ext) paste0(project_display_label(active()), "_method_agreement_", basename(selected()), ext)
  for (format in c("png", "pdf")) local({ ext <- format
    output[[paste0("tvcag_download_", ext)]] <- downloadHandler(filename = function() name(paste0(".", ext)), content = function(file) { s <- display()
      ytab_save_ggplot(file, plot_obj(), width = switch(s$width, wide = 13, full = 16, 10), height = max(4, as.numeric(sub("px", "", s$height)) / 75)) }) })
  for (t in c("genes", "pairs")) local({ tab <- t
    output[[paste0("tvcag_download_", tab)]] <- downloadHandler(filename = function() name(paste0("_", tab, ".tsv")), content = function(file) file.copy(file.path(selected(), paste0(tab, ".tsv")), file)) })
}
