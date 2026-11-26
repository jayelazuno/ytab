# Treated vs Control: two ways to call genes depleted or enriched in treated libraries relative to their controls.
# The local-abundance z-score (default) uses the Fitness run and results views (ytab fitness); resampling runs on a
# normalization run (ytab treated-vs-control) and has its own views here.
#
# The control-conditioned Gumbel is no longer offered. It scores runs of depleted sites rather than the gene's change
# in abundance, which made it the odd member of a chooser that reads as "three ways to do one thing", and the measured
# case for it was thin: calls on 2 of 14 public datasets, depletion only, and on Gale not one published hit that the
# z-score and resampling miss between them. It remains a CLI method (`ytab control-gumbel`) and a reported result; see
# docs/methods. Its name and colour stay below so comparisons already made with it still display.

tvc_method_choices <- c("Local-abundance z-score" = "zscore",
                        "Resampling" = "resampling")   # validation status: Methods drawer

# Labels for every method a stored run may name, including ones no longer offered.
tvc_method_labels <- c(tvc_method_choices, "Control-conditioned Gumbel (depletion only)" = "gumbel")

# Normalizations each method accepts: resampling takes the scalar methods and the beta-geometric corrections. The
# Gumbel's entry is kept so a stored Gumbel run still resolves its normalization.
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

# One drawer per method, beside that method's own controls. A single drawer describing both meant the page
# always carried the method you were not using, and its references and command line with it.
# One drawer per method, written as a list rather than a paragraph: these are four separate facts -- what the
# statistic is, what the noise model is, what goes in, what comes out -- and a reader wants whichever one they
# came for, not a sentence to parse to find it.
tvc_zscore_methods_text <- function() ytab_methods_drawer(
  tags$ul(class = "ytab-method-points",
    tags$li(tags$b("Statistic. "), "Replicate libraries of each condition are averaged; a gene's log2 fold change between the treated and control means is divided by the noise expected at that gene's abundance."),
    tags$li(tags$b("Noise model. "), "Fitted to the controls of each background: running SD of control\u2013control log2 ratios against abundance, SD = m1 + m2\u00b7A^m3."),
    tags$li(tags$b("Call threshold. "), "The 99th percentile of control\u2013control |z|."),
    tags$li(tags$b("Input. "), "BGC from the latest normalization run by default; gene totals of its normalized site counts, per million of the mean library total, so the 0.5 pseudocount means the same as with CPM. Any other normalization can be chosen."),
    tags$li(tags$b("Output. "), "One call per gene \u2014 depleted, enriched or none \u2014 from the replicate means. Per-pool z-scores are kept as a diagnostic."),
    tags$li(tags$b("Where it is weak. "), "The threshold is how much two controls disagree, not an error rate, and the treated libraries are noisier than the controls it was fitted on. On the benchmark it made 64 false calls per null experiment with BGC, and its false discovery rate with planted effects was 15\u201355%. It is kept because it ", tags$b("ranks"), " better than anything else tested (precision at K +0.086 over resampling): use it to prioritise genes, and resampling for a list with an error rate.")),
  ytab_reference("Follows Gale et al. 2020 (G3 10:3859\u20133870)."),
  tags$p("Command: ", tags$code("ytab fitness"), "."))

tvc_resampling_methods_text <- function() ytab_methods_drawer(
  tags$ul(class = "ytab-method-points",
    tags$li(tags$b("Statistic. "), "A permutation test on a gene's normalized site counts in control versus treated libraries, permuted within each pool when there are several. Two-sided, Benjamini\u2013Hochberg FDR."),
    tags$li(tags$b("Null. "), "Enumerated exactly by grid convolution; Monte Carlo is the fallback."),
    tags$li(tags$b("Input. "), "BGC by default. With BGC the treated libraries are recentred (trimmed sum) so that typical genes show no change."),
    tags$li(tags$b("Validation. "), "Passed its pre-registered calibration gates with BGC: 0.05\u20130.075 false calls per null experiment, false discovery rate 6.3% against a 5% target. Across the other normalizations it ranged 6\u201312%."),
    tags$li(tags$b("Where it is weak. "), "Power comes from the number of sites, so a gene with few sites cannot reach significance however far it moved.")),
  ytab_reference("Adapted from TRANSIT's resampling (DeJesus et al. 2015, PLoS Comput Biol 11:e1004401)."),
  tags$p("Command: ", tags$code("ytab treated-vs-control --method resampling"), "."))

# The same stages the z-score's Run tab has: choose the inputs, see what will be written, preview it, run it,
# watch it. Resampling had only a Run button, which on the costlier method is the wrong way round -- its default
# is exact enumeration over every gene, so knowing what a run will do before it starts is worth more here.
tvc_run_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel(uiOutput("tvc_run_title", inline = TRUE),
    uiOutput("tvc_inputs"),
    uiOutput("tvc_normalization_note"),
    uiOutput("tvc_comparison_selector"),
    uiOutput("tvc_selection_summary"),
    tags$div(class = "ytab-actions",
      actionButton("tvc_select_all", "Select all", class = "btn-secondary"),
      actionButton("tvc_clear_selection", "Clear selection", class = "btn-secondary")),
    numericInput("tvc_fdr", "FDR", .05, min = .001, max = .25, step = .01),
    radioButtons("tvc_execution_mode", "Execution mode",
                 c("Preview analysis" = "preview", "Run analysis" = "run"), "preview", inline = TRUE),
    uiOutput("tvc_run_summary"),
    uiOutput("tvc_run_action"),
    job_progress_ui("tvc_job")),
  main = tvc_resampling_methods_text())

# One results view for the three methods: the same controls, display options, plots and gene table, whichever method
# made the result. Only the plots that test something method-specific differ (the z-score's control-control null, the
# p-values of resampling and the Gumbel, the Gumbel's runs).
tvc_results_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Results",
    uiOutput("tvc_result_selector"),
    uiOutput("tvc_group_selector"),
    uiOutput("tvc_plot_selector"),
    uiOutput("tvc_gene_selector"),
    conditionalPanel("input.tvc_plot == 'volcano' || input.tvc_plot == 'ma' || input.tvc_plot == 'scatter' || input.tvc_plot == 'comet'",
      conditionalPanel("input.tvc_results_method != 'gumbel'",
        selectInput("tvc_direction", "Calls shown", c("Depleted and enriched" = "both", "Depleted" = "depleted", "Enriched" = "enriched"))),
      numericInput("tvc_label_n", "Label top hits", 10, min = 0, max = 50, step = 1),
      conditionalPanel("input.tvc_plot == 'comet'",
        helpText("On the comet this is per gene set, so a large family cannot take every label.")),
      textAreaInput("tvc_label_genes", "Also label", "", rows = 2, placeholder = "Gene names or IDs, separated by commas or new lines")),
  conditionalPanel("input.tvc_plot == 'comet'",
    uiOutput("tvc_comet_source"), uiOutput("tvc_comet_sets_selector"),
    textAreaInput("tvc_comet_genes", "Highlight my genes", "", rows = 2,
                  placeholder = "Gene names or IDs, separated by commas or new lines")),
    ytab_plot_downloads_ui("tvc", c(genes = "Download genes (TSV)"), capture = TRUE),
    ytab_plot_customization_controls("tvc", include_points = TRUE, include_labels = FALSE, open = FALSE)),
  main = tagList(uiOutput("tvc_summary"),
    ytab_plot_card("Treated vs control", uiOutput("tvc_plot_frame"), description = uiOutput("tvc_plot_caption")),
    ytab_table_card("Genes", DT::DTOutput("tvc_table")), uiOutput("tvc_results_methods_text")))

# The method chooser sits above the method's own card, so it needs the surrounding chrome a card gives: left
# alone it reads as a stray control floating on the page background.
tvc_run_tab_ui <- function() tagList(
  tags$div(class = "ytab-method-switch",
    radioButtons("tvc_method", "Method", tvc_method_choices, "zscore", inline = TRUE)),
  conditionalPanel("input.tvc_method == 'zscore'", fitness_run_ui()),
  conditionalPanel("input.tvc_method != 'zscore'", tvc_run_ui()))

tvc_results_tab_ui <- function() tagList(
  tags$div(class = "ytab-method-switch",
    radioButtons("tvc_results_method", "Method", c(tvc_method_choices, "Compare methods" = "agreement"), "zscore",
                 inline = TRUE)),
  conditionalPanel("input.tvc_results_method != 'agreement'", tvc_results_ui()),
  conditionalPanel("input.tvc_results_method == 'agreement'", tvc_agreement_ui()))

tvc_plot_choices <- function(method) switch(method %||% "resampling",
  zscore = c("Volcano" = "volcano", "MA plot" = "ma", "Condition vs control" = "scatter", "Comet (gene sets)" = "comet", "Control-control z (call threshold)" = "nullz", "Gene values" = "sites"),
  gumbel = c("Volcano" = "volcano", "Observed vs expected run" = "runs", "p-value histogram" = "pvalues", "Gene sites" = "sites"),
  c("Volcano" = "volcano", "MA plot" = "ma", "Condition vs control" = "scatter", "Comet (gene sets)" = "comet", "p-value histogram" = "pvalues", "Gene sites" = "sites"))

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
  # The comet: the same two means as the scatter, coloured by gene set rather than by call. Built by
  # comet_ggplot so that this tab, Multi-condition and any later one draw the identical plot.
  if (identical(choice, "comet")) {
    sets <- extra$comet
    if (is.null(sets) || !nrow(sets))
      return(empty("No gene sets selected.\nChoose one or more sets, or name genes to highlight."))
    p <- comet_ggplot(d, sets, s, theme, hover, extra$interactive)
    return(p %||% empty("No genes with both means in this comparison"))
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
  # The drawer follows the method whose result is on screen, so the page never carries the one you are not reading.
  output$tvc_results_methods_text <- renderUI(
    if (identical(input$tvc_results_method %||% "zscore", "resampling")) tvc_resampling_methods_text()
    else tvc_zscore_methods_text())
  output$tvc_run_title <- renderUI(names(tvc_method_labels)[tvc_method_labels == (input$tvc_method %||% "resampling")])

  # What this method will accept, said before the Run button rather than in a red toast after it.
  output$tvc_normalization_note <- renderUI({
    method <- input$tvc_method %||% "resampling"
    allowed <- tvc_allowed_methods[[method]] %||% character()
    if (!length(allowed)) return(NULL)
    # The short names only: Site Normalization is where these are chosen and explained.
    helpText(sprintf("%s accepts %s; BGC is the default. A normalization that reshapes counts rather than scaling them is not offered, because the test needs the counts it was given.",
                     if (identical(method, "gumbel")) "The Gumbel" else "Resampling",
                     paste(c(betageom = "BGC", aBGC = "aBGC", ttr = "TTR", totreads = "CPM",
                             nzmean = "non-zero mean", midlc = "MidLC", quantile = "quantile",
                             zinfnb = "ZINB")[allowed], collapse = ", ")))
  })
  output$tvc_selection_summary <- renderUI({
    d <- design(); chosen <- input$tvc_comparisons
    if (!nrow(d)) return(NULL)
    tags$p(class = "text-muted", sprintf("%d of %d comparisons selected.", length(chosen), nrow(d)))
  })
  observeEvent(input$tvc_select_all, {
    d <- design(); req(nrow(d) > 0)
    updateSelectizeInput(session, "tvc_comparisons", selected = d$comparison_id)
  }, ignoreInit = TRUE)
  observeEvent(input$tvc_clear_selection,
               updateSelectizeInput(session, "tvc_comparisons", selected = character()), ignoreInit = TRUE)

  # The run that is about to happen, including where it will be written. A run directory is generated, so this is
  # the only place the path can be seen before it exists.
  output$tvc_run_summary <- renderUI({
    method <- input$tvc_method %||% "resampling"; norm <- tvc_norm(); chosen <- input$tvc_comparisons
    if (!length(chosen) || is.null(norm)) return(NULL)
    tags$dl(class = "ytab-run-summary",
      tags$dt("Method"), tags$dd(names(tvc_method_labels)[tvc_method_labels == method]),
      tags$dt("Comparisons"), tags$dd(sprintf("%d selected", length(chosen))),
      tags$dt("Normalization"), tags$dd(as.character(norm_method_label(norm$method))),
      tags$dt("FDR"), tags$dd(format(input$tvc_fdr %||% .05)),
      tags$dt("Output"), tags$dd(tags$code(file.path("treated_vs_control", paste0(method, "_<timestamp>")))))
  })
  output$tvc_run_action <- renderUI({
    preview <- identical(input$tvc_execution_mode %||% "preview", "preview")
    actionButton("tvc_run", if (preview) "Preview analysis" else "Run analysis", class = "btn-primary")
  })
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
    args <- c("--project-config", active_project_path(), "--method", method, ytab_flag("--normalization-run", norm$run, required = TRUE),
              ytab_flag("--normalization", norm$method, required = TRUE), "--comparisons", paste(comparisons, collapse = ","),
              "--fdr", as.character(input$tvc_fdr %||% .05), "--output", out)
    # Preview is the pipeline's own --dry-run: it validates the inputs and reports what it would write.
    if (identical(input$tvc_execution_mode %||% "preview", "preview")) args <- c(args, "--dry-run")
    start("ytab_treated_vs_control.py", args, "treated_vs_control",
          sprintf("%s %s.", names(tvc_method_labels)[tvc_method_labels == method],
                  if (identical(input$tvc_execution_mode %||% "preview", "preview")) "preview started" else "started"))
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
      manifest <- list(method = "zscore", method_label = names(tvc_method_labels)[tvc_method_labels == "zscore"], current = isTRUE(r$current),
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
  # Gene sets shipped with the project's species, and the sets the user chose to overlay. Assembled here so that
  # choosing "Comet" draws a plot straight away: the axes are columns the results already carry, and the overlay is
  # a lookup over the same GMT and orthology files `ytab comet` reads, which is why the two agree gene for gene
  # (verified on Gale 2020: 5,234 of 5,234). No statistic is computed; this is colouring, like colouring by call.
  comet_species <- reactive(list(root = active()$repo_root,
                                 dir = (active()$reference %||% list())$species_dir %||% ""))
  comet_files <- reactive(comet_set_files(comet_species()$root, comet_species()$dir))
  comet_sets <- reactive({
    files <- comet_files(); if (!length(files)) return(list())
    chosen <- input$tvc_comet_file %||% unname(files)[[1L]]
    comet_read_gmt(if (chosen %in% files) chosen else unname(files)[[1L]])
  })
  # Only offered when the species ships more than one collection; one collection needs no choice.
  output$tvc_comet_source <- renderUI({ files <- comet_files()
    if (length(files) < 2L) return(NULL)
    selectInput("tvc_comet_file", "Gene set file", choices = files, selected = input$tvc_comet_file) })
  # A searchable multi-select: the number of sets is data-length, so it must not be a fixed list of checkboxes.
  output$tvc_comet_sets_selector <- renderUI({ sets <- comet_sets()
    if (!length(sets)) return(helpText(sprintf("No gene sets for this species. Add a GMT file under %s/genesets/.",
                                               comet_species()$dir)))
    kept <- intersect(input$tvc_comet_sets %||% utils::head(names(sets), 4L), names(sets))
    selectizeInput("tvc_comet_sets", "Overlay sets", choices = names(sets), selected = kept, multiple = TRUE,
                   options = list(plugins = list("remove_button"), placeholder = "Search gene sets")) })

  comet <- reactive({
    req(identical(input$tvc_plot, "comet"))
    g <- group_genes(); req(nrow(g) > 0)
    sets <- comet_sets()
    comet_overlay(g, sets, intersect(input$tvc_comet_sets %||% utils::head(names(sets), 4L), names(sets)),
                  comet_orthologs(comet_species()$root, comet_species()$dir, g$ORF),
                  input$tvc_comet_genes, input$tvc_label_n %||% 10L,
                  direction = input$tvc_direction %||% "both")
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
    comet = "Each gene's treated mean against its control mean, coloured by gene set: the dense body on the diagonal is the genes that did not respond, and genes spraying away from it did. Sets and labels come from the run's comet table (ytab comet); labelled names are the ortholog's where one exists.",
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
      comet = if (identical(choice, "comet")) comet(),
      label_n = input$tvc_label_n %||% 10, label_genes = input$tvc_label_genes, unit = d$manifest$value_unit %||% "", pseudocount = d$manifest$pseudocount,
      pairs = if (identical(choice, "sites") && identical(d$method, "zscore")) zscore_pairs(),
      direction = if (identical(d$method, "gumbel")) "both" else input$tvc_direction %||% "both", library_labels = library_labels(),
      interactive = interactive)) }
  current_plot <- reactive(build_plot(TRUE))
  output$tvc_plotly <- plotly::renderPlotly(tvc_plotly(current_plot()))
  # Shaped once, for the screen and the download alike. The download used to write every gene and every column
  # of the run, ignoring the chosen comparison -- a different table from the one being looked at.
  tvc_genes_display <- reactive({
    g <- tvc_tested(group_genes()); g <- g[order(-tvc_significance(g, data()$method)), , drop = FALSE]
    keep <- intersect(c("ORF", "Name", "n_sites", "n_eligible", "n_depleted", "run_span_bp", "mean_control", "mean_treatment",
                        "log2FC", "log2_scaled_ratio", "z", "p_value", "q_bh", "call"), names(g))
    g <- g[, keep, drop = FALSE]
    for (column in intersect(c("mean_control", "mean_treatment", "log2FC", "log2_scaled_ratio", "z"), names(g))) g[[column]] <- signif(g[[column]], 3)
    for (column in intersect(c("p_value", "q_bh"), names(g))) g[[column]] <- formatC(g[[column]], format = "g", digits = 3)
    g
  })
  output$tvc_table <- DT::renderDT(ytab_datatable(tvc_genes_display(), page_length = 15L))
  name <- function(ext) paste0(project_display_label(active()), "_treated_vs_control_", basename(selected()), ext)
  for (format in c("png", "pdf")) local({ ext <- format
    output[[paste0("tvc_download_", ext)]] <- downloadHandler(filename = function() name(paste0(".", ext)), content = function(file) { s <- display()
      ytab_save_ggplot(file, build_plot(FALSE), width = switch(s$width, wide = 13, full = 16, 10), height = max(4, as.numeric(sub("px", "", s$height)) / 75)) }) })
    ytab_plot_capture_server(input, output, session, "tvc", "tvc_plotly", function() sub("[.]$", "", name("")))
  output$tvc_download_genes <- downloadHandler(filename = function() name("_genes.tsv"),
    content = function(file) utils::write.table(tvc_genes_display(), file, sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
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
      uiOutput("tvc_ag_direction_selector")),
    conditionalPanel("input.tvc_ag_plot == 'effect' || input.tvc_ag_plot == 'score'", uiOutput("tvc_ag_pair_selector")),
    ytab_plot_downloads_ui("tvcag", c(genes = "Download per-gene table (TSV)", pairs = "Download agreement (TSV)"), capture = TRUE),
    ytab_plot_customization_controls("tvcag", include_points = TRUE, include_labels = FALSE, open = FALSE)),
  main = tagList(uiOutput("tvc_ag_summary"),
    ytab_plot_card("Method agreement", uiOutput("tvc_ag_plot_frame"), description = uiOutput("tvc_ag_caption")),
    ytab_table_card("Agreement", DT::DTOutput("tvc_ag_table")),
    ytab_methods_drawer(
      tags$p("Each pair of methods is compared on the genes both tested, per background and treated condition. ",
             tags$b("Calls"), " are compared by recovery: of one method's calls in a direction, the share the other method also calls (with published calls as the first method, this is the recall of the published hits), shown with how many genes each calls. ",
             tags$b("Effects"), " and ", tags$b("signed significance"), " (z for the z-score; sign(effect) × −log10 p for resampling; −(−log10 p) for the depletion-only Gumbel) are compared by Spearman rank correlation, since the methods' scales and tails differ. A method that tests depletion only is left out of the enriched comparisons, and the directions on offer follow from the methods in the comparison."),
      tags$p("Recovery between two methods is often well below 100%, and reading it needs the rank correlation beside it."),
      tags$p(tags$b("The same quantity, a different threshold. "),
             "The z-score and resampling both score a gene's overall change in abundance, so they rank genes nearly identically (Spearman of effect near 1) and still share a minority of their calls: recovery depends on where each one's cutoff falls in that shared ranking. ",
             "A fixed z cutoff answers to effect size alone, so it keeps genes whose few insertion sites leave a permutation test short of significance; a permutation test turns many sites into power, so it keeps genes whose effect is modest. ",
             "The two surpluses are not equally trustworthy: the z-score's calls are not error-controlled and it makes many calls on simulated null comparisons, where resampling makes almost none, so the genes only the z-score calls are the ones to treat with suspicion."),
      tags$p(tags$b("A different quantity, in earlier comparisons. "),
             "A comparison that includes a stored control-conditioned Gumbel result is not reading a threshold difference at all: that test scores the longest run of consecutive depleted sites rather than the gene's overall shift, so it finds a different shape of depletion \u2014 a gene with a contiguous dead stretch but only a mild average change is its kind, and a gene whose sites are all modestly down with no dead run is resampling's. ",
             "A high Spearman of effect next to a low Spearman of signed significance is the sign of it: the pair agrees on how much each gene changed and disagrees on what counts as evidence. ",
             "Low recovery there is complementary coverage rather than error, but it does not make the union a hit list: joining two error-controlled lists keeps neither one's error rate."),
      tags$p("Command: ", tags$code("ytab tvc-agreement --zscore RUN --resampling RUN --output DIR"), "."))))

# Chosen directions, in a fixed order, as panel titles.
tvc_agreement_directions <- c(depleted = "Depleted calls", enriched = "Enriched calls")
tvc_agreement_direction_order <- function(directions) intersect(names(tvc_agreement_directions), directions)
# A light-to-dark ramp around the focal colour. White-to-one-colour leaves 30-60% — where recovery usually sits —
# hard to tell apart, so the stops are spread over that range instead.
TVC_RECOVERY_RAMP <- c("#f4fafa", "#d3e9ea", "#a5d4d6", "#6ab8bc", "#2f979d", "#087f85", "#055c61")

# Recovery matrix: tile (row, column) = of the row method's calls in that direction, the share the column method also
# calls (the published row is the recall of the published hits); the diagonal gives each method's own number of calls.
# Both directions are drawn as side-by-side panels, and a depletion-only method's enriched cells are left empty.
tvc_agreement_matrix <- function(pairs, directions, s) {
  validate(need(nrow(pairs) > 0, "No method pairs to compare in this comparison."))
  directions <- tvc_agreement_direction_order(directions)
  m <- unique(c(pairs$method_a, pairs$method_b)); lab <- unname(tvc_agreement_methods[m])
  side <- function(dir, a, b, suffix) { col <- function(k) pairs[[paste0(dir, "_", k)]] %||% rep(NA_real_, nrow(pairs))
    data.frame(direction = dir, row = pairs[[a]], col = pairs[[b]], calls = col(paste0("n_", suffix[[1L]])), shared = col("n_both"),
               value = col(paste0(suffix[[1L]], "_recovered_by_", suffix[[2L]])), rho = pairs$spearman_effect, rho_s = pairs$spearman_score) }
  off <- do.call(rbind, lapply(directions, function(dir)
    rbind(side(dir, "method_a", "method_b", c("a", "b")), side(dir, "method_b", "method_a", c("b", "a")))))
  own <- do.call(rbind, lapply(directions, function(dir) do.call(rbind, lapply(m, function(k) {
    x <- suppressWarnings(max(off$calls[off$direction == dir & off$row == k], na.rm = TRUE))
    data.frame(direction = dir, row = k, col = k, calls = if (is.finite(x)) x else NA_real_) }))))
  place <- function(d) { d$r <- factor(unname(tvc_agreement_methods[d$row]), rev(lab)); d$c <- factor(unname(tvc_agreement_methods[d$col]), lab)
    d$panel <- factor(unname(tvc_agreement_directions[d$direction]), unname(tvc_agreement_directions[directions])); d }
  off <- place(off); own <- place(own)
  count <- function(x) formatC(x, format = "d", big.mark = ",")
  # The percentage carries the cell; the counts behind it are in the tooltip, so three lines of text do not compete.
  off$label <- ifelse(is.na(off$value), "–", sprintf("%.0f%%", 100 * off$value))
  off$tip <- ifelse(is.na(off$value),
    sprintf("%s does not call enrichment", unname(tvc_agreement_methods[off$row])),
    sprintf("Of %s's %s %s calls, %s (%.0f%%) are also called by %s<br>Spearman: effect %.3f, signed significance %.2f",
            unname(tvc_agreement_methods[off$row]), count(off$calls), off$direction, count(off$shared), 100 * off$value,
            unname(tvc_agreement_methods[off$col]), off$rho, off$rho_s))
  own$label <- ifelse(is.na(own$calls), "–", paste0(count(own$calls), "\ncalls"))
  own$tip <- ifelse(is.na(own$calls), sprintf("%s does not call enrichment", unname(tvc_agreement_methods[own$row])),
                    sprintf("%s calls %s genes in this direction", unname(tvc_agreement_methods[own$row]), count(own$calls)))
  # Dark tiles need light text and light tiles dark text, or the percentage disappears at one end of the ramp.
  off$ink <- ifelse(!is.na(off$value) & off$value > .72, "white", "#10292b")
  ggplot2::ggplot(off, ggplot2::aes(c, r, fill = value, text = tip)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.4) +
    ggplot2::geom_tile(data = own, ggplot2::aes(c, r, text = tip), fill = "grey95", colour = "white", linewidth = 1.4, inherit.aes = FALSE) +
    ggplot2::geom_text(data = own, ggplot2::aes(c, r, label = label), colour = "grey30", size = s$text_size / 3.8, fontface = "bold", lineheight = .95, inherit.aes = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = label, colour = ink), size = s$text_size / 2.9, fontface = "bold", show.legend = FALSE) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradientn(colours = TVC_RECOVERY_RAMP, limits = c(0, 1), na.value = "grey95", breaks = seq(0, 1, .25),
                                  labels = function(x) paste0(100 * x, "%"), name = "Also called\nby the column") +
    ggplot2::coord_equal() + ggplot2::labs(x = "Also called by", y = "Calls of") +
    (if (length(directions) > 1L) ggplot2::facet_wrap(~panel) else ggplot2::facet_null()) +
    # Bold black axis text and the blank ticks a tile plot wants come from ytab_style_ggplot.
    tvc_theme(s) + ggplot2::theme(panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)), panel.spacing = grid::unit(1.4, "lines"),
      legend.key.height = grid::unit(1.1, "lines"))
}

# Venn of the calls in one direction for two or three methods: region counts are the genes called by exactly that set of
# the chosen methods (from sets.tsv, genes tested by every method of the comparison).
tvc_agreement_venn <- function(sets, methods, directions, s, depletion_only = character()) {
  validate(need(length(methods) %in% 2:3, "Choose two or three methods for the Venn diagram."))
  directions <- tvc_agreement_direction_order(directions)
  # A depletion-only method keeps its circle in the enriched panel, because the geometry is shared between panels.
  # Drawn faint and dashed with its regions dashed out, so it reads as "cannot call this" and not "called nothing".
  inactive <- function(dir) if (identical(dir, "enriched")) intersect(methods, depletion_only) else character()
  circle <- function(x, y, r = 1.15, k) { t <- seq(0, 2 * pi, length.out = 120); data.frame(x = x + r * cos(t), y = y + r * sin(t), method = k) }
  centres <- if (length(methods) == 2L) list(c(-.6, 0), c(.6, 0)) else list(c(-.65, .45), c(.65, .45), c(0, -.6))
  layout <- if (length(methods) == 2L) data.frame(k = c("1", "2", "12"), x = c(-1.15, 1.15, 0), y = 0) else
    data.frame(k = c("1", "2", "3", "12", "13", "23", "123"), x = c(-1.2, 1.2, 0, 0, -.75, .75, 0), y = c(.85, .85, -1.25, .95, -.25, -.25, .1))
  # The geometry is the same in every panel; only the region counts change with the direction.
  panel <- function(d, dir) { d$panel <- factor(unname(tvc_agreement_directions[dir]), unname(tvc_agreement_directions[directions])); d }
  per <- function(dir) { x <- sets[sets$direction == dir, , drop = FALSE]; dead <- match(inactive(dir), methods)
    member <- lapply(strsplit(x$methods, "&", fixed = TRUE), function(k) intersect(k, methods))
    key <- vapply(member, function(k) paste(sort(match(k, methods)), collapse = ""), "")
    counts <- tapply(x$genes, key, sum); get <- function(k) if (k %in% names(counts)) counts[[k]] else 0L
    r <- layout; involved <- vapply(r$k, function(k) any(as.integer(strsplit(k, "")[[1L]]) %in% dead), TRUE)
    named <- vapply(r$k, function(k) paste(tvc_agreement_methods[methods[as.integer(strsplit(k, "")[[1L]])]], collapse = " and "), "")
    r$label <- ifelse(involved, "–", formatC(vapply(r$k, get, 0), format = "d", big.mark = ","))
    r$ink <- ifelse(involved, "grey55", "black")
    r$tip <- ifelse(involved, sprintf("%s does not call enrichment", paste(tvc_agreement_methods[inactive(dir)], collapse = " and ")),
                    sprintf("%s: %s genes called %s", named, r$label, dir))
    panel(r, dir) }
  regions <- do.call(rbind, lapply(directions, per))
  circles <- do.call(rbind, lapply(directions, function(dir) { d <- panel(do.call(rbind, lapply(seq_along(methods), function(i)
    circle(centres[[i]][1], centres[[i]][2], k = methods[[i]]))), dir); d$off <- d$method %in% inactive(dir); d }))
  names_at <- do.call(rbind, lapply(directions, function(dir) panel(do.call(rbind, lapply(seq_along(methods), function(i)
    data.frame(x = centres[[i]][1] * 1.9, y = centres[[i]][2] * 1.9 + if (centres[[i]][2] < 0) -.55 else .5,
               label = tvc_agreement_methods[[methods[[i]]]],
               ink = if (methods[[i]] %in% inactive(dir)) "grey55" else "black"))), dir)))
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = circles, ggplot2::aes(x, y, group = interaction(method, panel), fill = method, alpha = off, linetype = off),
                          colour = "black", linewidth = .6) +
    ggplot2::geom_text(data = regions, ggplot2::aes(x, y, label = label, colour = ink, text = tip), size = s$text_size / 2.6, fontface = "bold") +
    ggplot2::geom_text(data = names_at, ggplot2::aes(x, y, label = label, colour = ink), size = s$text_size / 3, fontface = "bold") +
    ggplot2::scale_fill_manual(values = tvc_agreement_colours, guide = "none") + ggplot2::scale_colour_identity() +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = .25, `TRUE` = .05), guide = "none") +
    ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") + ggplot2::coord_equal() +
    (if (length(directions) > 1L) ggplot2::facet_wrap(~panel) else ggplot2::facet_null()) +
    ggplot2::theme_void() + ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)))
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
tvc_agreement_upset <- function(sets, methods, directions, s, depletion_only = character()) {
  directions <- tvc_agreement_direction_order(directions)
  d <- sets[sets$direction %in% directions, , drop = FALSE]
  validate(need(nrow(d) > 0, "No calls to compare in this direction."))
  d$panel <- factor(unname(tvc_agreement_directions[d$direction]), unname(tvc_agreement_directions[directions]))
  d <- d[order(d$panel, -d$genes), , drop = FALSE]
  d$set <- factor(seq_len(nrow(d)), levels = seq_len(nrow(d)))
  members <- strsplit(d$methods, "&", fixed = TRUE)
  # "Called by all" means every method that can call in that direction, so a depletion-only method is not counted
  # against the enriched sets.
  possible <- vapply(d$direction, function(dir) length(if (identical(dir, "enriched")) setdiff(methods, depletion_only) else methods), 0L)
  d$group <- ifelse(d$n_methods == possible, "all", ifelse(d$n_methods == 1, vapply(members, `[`, "", 1L), "some"))
  fills <- c(all = YTAB_FOCAL, some = YTAB_COMPARATOR, tvc_agreement_colours[methods])
  # One y scale across the panels, so the method rows line up and the bars stay comparable between directions.
  ymax <- max(d$genes); step <- ymax * .1; rows <- -step * seq_along(methods)
  dots <- do.call(rbind, lapply(seq_along(methods), function(i)
    data.frame(set = d$set, panel = d$panel, y = rows[[i]], on = vapply(members, function(m) methods[[i]] %in% m, TRUE))))
  links <- do.call(rbind, lapply(seq_len(nrow(d)), function(k) { on <- rows[methods %in% members[[k]]]
    if (length(on) > 1) data.frame(set = d$set[[k]], panel = d$panel[[k]], ymin = min(on), ymax = max(on)) }))
  d$tip <- sprintf("Called %s by: %s<br>%s genes", d$direction, gsub("&", ", ", d$methods), formatC(d$genes, format = "d", big.mark = ","))
  ticks <- pretty(c(0, ymax)); ticks <- ticks[ticks <= ymax * 1.05]
  p <- ggplot2::ggplot(d, ggplot2::aes(set, genes)) + ggplot2::geom_col(ggplot2::aes(fill = group, text = tip), width = .7) +
    ggplot2::geom_text(ggplot2::aes(y = genes + ymax * .045, label = formatC(genes, format = "d", big.mark = ",")), size = s$text_size / 3.8, fontface = "bold") +
    ggplot2::geom_point(data = dots, ggplot2::aes(set, y, colour = on), size = 3.2, inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey85"), guide = "none") +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::scale_y_continuous(breaks = c(ticks, rows), labels = c(formatC(ticks, format = "d", big.mark = ","), tvc_agreement_methods[methods]), expand = ggplot2::expansion(mult = c(.06, .1))) +
    ggplot2::scale_x_discrete(labels = NULL) + ggplot2::labs(x = "Methods that call the genes (exactly this set)", y = "Genes") +
    (if (length(directions) > 1L) ggplot2::facet_wrap(~panel, scales = "free_x") else ggplot2::facet_null()) +
    tvc_theme(s) + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)), panel.spacing = grid::unit(1.2, "lines"))
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
    # A stored Gumbel result can still be compared, but only where one exists: the method is no longer offered, so the
    # selector appears for projects that already have such a run and is absent everywhere else.
    gumbel <- pick("gumbel")
    tagList(selectInput("tvc_ag_zscore", "z-score result", zc), selectInput("tvc_ag_resampling", "Resampling result", pick("resampling")),
            if (length(gumbel) > 1L) selectInput("tvc_ag_gumbel", "Gumbel result (earlier runs only)", gumbel),
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
  # The directions on offer come from the comparison itself: a method that tests depletion only has no enriched calls,
  # so unless two of the methods can call enrichment there is nothing to choose between.
  all_methods <- reactive(unlist(data()$manifest$methods) %||% character())
  depletion_only <- reactive(unlist(data()$manifest$depletion_only) %||% character())
  enriched_methods <- reactive(setdiff(all_methods(), depletion_only()))
  output$tvc_ag_direction_selector <- renderUI({ req(nrow(runs()) > 0)
    if (length(enriched_methods()) < 2L) return(helpText(sprintf("Depleted calls only: %s tests depletion.",
      if (length(depletion_only())) paste(tvc_agreement_methods[depletion_only()], collapse = " and ") else "this comparison")))
    selectInput("tvc_ag_direction", "Calls", c("Both, side by side" = "both", "Depleted only" = "depleted", "Enriched only" = "enriched")) })
  direction <- reactive(if (length(enriched_methods()) < 2L) "depleted" else input$tvc_ag_direction %||% "both")
  directions <- reactive(if (identical(direction(), "both")) c("depleted", "enriched") else direction())
  output$tvc_ag_summary <- renderUI({
    if (!nrow(runs())) return(ytab_empty_state("No method comparison yet", "Choose a result of each method and click Compare; the Gumbel is optional."))
    p <- data()$pairs; p <- p[p$group == group(), , drop = FALSE]; req(nrow(p) > 0)
    pct <- function(x) ifelse(is.na(x), "–", sprintf("%.0f%%", 100 * x))
    count <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
    # Recovery and the effect-rank correlation sit on the same card: a pair can rank the genes almost identically and
    # still share a minority of its calls, and only the two numbers together say which is happening.
    cards <- lapply(seq_len(nrow(p)), function(i) { q <- p[i, , drop = FALSE]
      a <- tvc_agreement_methods[[q$method_a]]; b <- tvc_agreement_methods[[q$method_b]]
      ytab_metric_card(paste(a, "vs", b), sprintf("%s / %s", pct(q$depleted_a_recovered_by_b), pct(q$depleted_b_recovered_by_a)),
        sprintf("of %s's %s and %s's %s depleted calls, %s shared · effect ranks ρ %.3f · signed significance ρ %.2f · %s genes tested by both",
                a, count(q$depleted_n_a), b, count(q$depleted_n_b), count(q$depleted_n_both), q$spearman_effect, q$spearman_score,
                count(q$genes_tested_by_both))) })
    # Low recovery has two causes worth telling apart, and the two rank correlations tell them apart: a pair that
    # agrees on the effect and on significance differs only in where its cutoff falls, while a pair that agrees on the
    # effect but not on significance is scoring different evidence altogether.
    worst <- suppressWarnings(pmin(p$depleted_a_recovered_by_b, p$depleted_b_recovered_by_a, na.rm = TRUE))
    low <- !is.na(p$spearman_effect) & p$spearman_effect >= .9 & !is.na(worst) & worst < .8
    named <- function(i) paste(sprintf("%s and %s", tvc_agreement_methods[p$method_a[i]], tvc_agreement_methods[p$method_b[i]]), collapse = "; ")
    threshold <- low & !is.na(p$spearman_score) & p$spearman_score >= .8
    evidence <- low & (is.na(p$spearman_score) | p$spearman_score < .8)
    tags$div(tags$div(class = "ytab-stat-grid", cards),
      if (any(threshold)) tags$p(class = "text-muted", style = "margin-top:.75em",
        tags$b(named(threshold)), " rank these genes almost identically and agree on significance, yet share fewer of their calls: they differ on where the cutoff falls, not on which genes respond."),
      if (any(evidence)) tags$p(class = "text-muted", style = "margin-top:.4em",
        tags$b(named(evidence)), " agree on how much each gene changed (effect ρ ≥ 0.9) but not on what counts as evidence (signed significance ρ < 0.8): they score different things, so their calls overlap only in part, by design."),
      if (any(low)) tags$p(class = "text-muted", style = "margin-top:.4em",
        "The methods note below says what drives each, and which surplus calls carry error control."))
  })
  display <- reactive(ytab_plot_customization_values(input, "tvcag"))
  plot_obj <- reactive({
    d <- data(); s <- display(); g <- group(); choice <- input$tvc_ag_plot %||% "matrix"; dirs <- directions()
    if (choice == "matrix") return(tvc_agreement_matrix(d$pairs[d$pairs$group == g, , drop = FALSE], dirs, s))
    sets <- d$sets[d$sets$group == g, , drop = FALSE]
    # On the enriched panel alone a depletion-only method would be an empty circle, so it is dropped; with both
    # directions the geometry is shared between the panels, so it stays and simply calls nothing enriched.
    if (choice == "venn") { m <- input$tvc_ag_venn %||% utils::head(all_methods(), 3)
      if (identical(dirs, "enriched")) m <- setdiff(m, depletion_only())
      return(tvc_agreement_venn(sets, m, dirs, s, depletion_only())) }
    if (choice == "sets") { m <- all_methods(); if (identical(dirs, "enriched")) m <- setdiff(m, depletion_only())
      return(tvc_agreement_upset(sets, m, dirs, s, depletion_only())) }
    pair <- strsplit(input$tvc_ag_pair %||% paste(unlist(d$manifest$methods)[1:2], collapse = "|"), "|", fixed = TRUE)[[1L]]
    tvc_agreement_scatter(d$genes[d$genes$group == g, , drop = FALSE], pair[[1L]], pair[[2L]], choice, s)
  })
  output$tvc_ag_plot_frame <- renderUI({ req(nrow(runs()) > 0); s <- display(); ytab_plot_frame(plotly::plotlyOutput("tvc_ag_plotly", height = s$height), s$width) })
  output$tvc_ag_plotly <- plotly::renderPlotly({ built <- ytab_ggplotly(plot_obj())
    # A Venn has no axes, and faceting it gives a second pair (xaxis2/yaxis2), so every axis present is hidden.
    if (identical(input$tvc_ag_plot, "venn")) { built <- plotly::plotly_build(built)
      for (ax in grep("^[xy]axis[0-9]*$", names(built$x$layout), value = TRUE))
        built$x$layout[[ax]] <- utils::modifyList(built$x$layout[[ax]], list(visible = FALSE)) }
    built })
  output$tvc_ag_caption <- renderUI({ both <- length(directions()) > 1L
    panels <- if (both) " Depleted and enriched calls are drawn as side-by-side panels." else ""
    switch(input$tvc_ag_plot %||% "matrix",
    matrix = paste0("Of the row method's calls, the share the column method also calls; the published row is the recall of the published hits. The diagonal gives each method's own number of calls, and the counts behind each percentage are in the tooltip, with the pair's rank correlations.",
                    panels, " A method that tests depletion only has no enriched cells, shown as a dash."),
    venn = paste0("Genes called by each combination of the chosen methods; only genes tested by every method of the comparison.",
                  panels, if (both) " A method that tests depletion only keeps its circle in the enriched panel, where it calls nothing." else ""),
    effect = "Each gene's effect in one method against the other, over genes both tested; colours show which method called it (either direction; depleted only when the Gumbel is one of the pair).",
    score = "Each gene's signed significance in one method against the other (z; sign(effect) × −log10 p; the Gumbel's depletion as −(−log10 p)).",
    sets = paste0("Genes called by exactly the marked set of methods, largest sets first; only genes tested by every method shown.",
                  panels, if (both) " The panels share one count axis, so the bars are comparable between directions." else "")) })
  # The screen shows the chosen group's pairs under readable headings; the download copied every group raw.
  tvc_ag_display <- reactive({ p <- data()$pairs; p <- p[p$group == group(), , drop = FALSE]; req(nrow(p) > 0)
    out <- data.frame(Methods = paste(tvc_agreement_methods[p$method_a], "vs", tvc_agreement_methods[p$method_b]), `Genes tested by both` = p$genes_tested_by_both,
      `Spearman (effect)` = round(p$spearman_effect, 3), `Spearman (significance)` = round(p$spearman_score, 3),
      `Depleted: A / B / shared` = sprintf("%s / %s / %s", p$depleted_n_a, p$depleted_n_b, p$depleted_n_both),
      `Depleted: A recovered by B` = sprintf("%.0f%%", 100 * p$depleted_a_recovered_by_b), `Depleted: B recovered by A` = sprintf("%.0f%%", 100 * p$depleted_b_recovered_by_a),
      `Enriched: A / B / shared` = ifelse(is.na(p$enriched_n_a), "depletion only", sprintf("%s / %s / %s", p$enriched_n_a, p$enriched_n_b, p$enriched_n_both)), check.names = FALSE)
    out })
  output$tvc_ag_table <- DT::renderDT(ytab_datatable(tvc_ag_display(), page_length = 10L))
  name <- function(ext) paste0(project_display_label(active()), "_method_agreement_", basename(selected()), ext)
  for (format in c("png", "pdf")) local({ ext <- format
    output[[paste0("tvcag_download_", ext)]] <- downloadHandler(filename = function() name(paste0(".", ext)), content = function(file) { s <- display()
      ytab_save_ggplot(file, plot_obj(), width = switch(s$width, wide = 13, full = 16, 10), height = max(4, as.numeric(sub("px", "", s$height)) / 75)) }) })
  # "pairs" is the table on screen, so it is written as shown; "genes" has no table of its own here and stays
  # the run's own export.
  ytab_plot_capture_server(input, output, session, "tvcag", "tvc_ag_plotly",
                           function() sub("[.]$", "", name("")))
  output$tvcag_download_pairs <- downloadHandler(filename = function() name("_pairs.tsv"),
    content = function(file) utils::write.table(tvc_ag_display(), file, sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
  output$tvcag_download_genes <- downloadHandler(filename = function() name("_genes.tsv"),
    content = function(file) file.copy(file.path(selected(), "genes.tsv"), file))
}
