# Essentiality module: label-free window HMM (default), plus the hybrid model when curated labels exist. The HMM is
# transposon-agnostic: it observes genomic windows, not a fixed site grid, and has been run on Hermes and Ac/Ds libraries.
# (`ytab essentiality`). Combine Parents, Combined Summary and the Classifier are kept as the hybrid's inputs.

essentiality_state_labels <- c(ES = "Essential (ES)", GD = "Growth defect (GD)", NE = "Non-essential (NE)", GA = "Growth advantage (GA)")
essentiality_state_colours <- c(ES = "#b2182b", GD = "#ef8a62", NE = "#67a9cf", GA = "#1a9850")
# HMM vs hybrid agreement: agreeing calls take the state colours; disagreements their own
essentiality_agreement_labels <- c(`consensus essential` = "Consensus essential", `HMM only` = "HMM only", `hybrid only` = "Hybrid only", `not essential` = "Not essential")
essentiality_agreement_colours <- c(`consensus essential` = "#b2182b", `HMM only` = "#984ea3", `hybrid only` = "#ff7f00", `not essential` = "#67a9cf")

# One parent library (with a hit file) per pool of the strain with most parents. Parent detection uses "parent"
# metadata only, so mock controls split for a treatment are not chosen.
essentiality_default_libraries <- function(inventory) {
  x <- inventory[inventory$Included & inventory$Detected_parent & inventory$Eligible == "Yes", , drop = FALSE]
  if (!nrow(x)) return(character())
  strain <- names(sort(table(x$Background), decreasing = TRUE))[[1L]]
  x <- x[x$Background == strain, , drop = FALSE]
  x$Sample[!duplicated(x$Pool)]
}

essentiality_run_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Essentiality",
    # The model comes first, as it does on every other analysis tab. What used to lead here was the HMM's
    # fit-or-reuse choice -- an execution detail -- while the hybrid, the second model, was a checkbox below a
    # collapsed drawer, so the tab read as though the HMM were the only model on offer.
    uiOutput("ess_model_selector"),
    # Choosing a model already says the HMM will be fitted, so there is no "fit a new one" to pick: fitting is
    # what the Run button does. The only departure worth a control is reusing a fit that already exists, which
    # skips the slow step -- and it is offered only when this dataset has one.
    uiOutput("ess_hmm_reuse_toggle"),
    conditionalPanel("!input.ess_reuse_fit", uiOutput("ess_library_selector"),
      helpText("Choose one untreated library per independent pool, all from one strain. Several pools are fitted jointly.")),
    conditionalPanel("input.ess_reuse_fit", uiOutput("ess_hmm_reuse_selector")),
    conditionalPanel("input.ess_model === 'hybrid'", uiOutput("ess_hybrid_inputs")),
    tags$details(class = "ytab-more-options", tags$summary("Settings"),
      numericInput("ess_window", "Window (bp)", 100, min = 50, step = 50),
      numericInput("ess_min_reads", "Reads for a site to count as an insertion", 1, min = 1, step = 1)),
    actionButton("ess_run", "Run essentiality", class = "btn-primary")),
  main = tagList(
    uiOutput("ess_run_status"),
    ytab_methods_drawer(
      tags$table(class = "table table-sm",
        tags$thead(tags$tr(tags$th("Model"), tags$th("Needs labels?"), tags$th("Use when"), tags$th("Main caveat"))),
        tags$tbody(
          tags$tr(tags$td(tags$b("HMM (default)")), tags$td("No"), tags$td("Always; the only option for organisms without essentiality labels"),
                  tags$td("Lower precision: accessibility cold spots and severe growth defects can look essential")),
          tags$tr(tags$td(tags$b("Hybrid")), tags$td("Yes"), tags$td("Curated labels exist for the species or a close relative"),
                  tags$td("Inherits the biases of its labels")))),
      tags$p(tags$b("HMM. "), "The genome is cut into windows (100 bp). Each window's insertion count, across the libraries fitted jointly, is modelled as coming from one of four hidden states: essential (ES), growth defect (GD), non-essential (NE) and growth advantage (GA), with window covariates for insertion-site preference (TTAA-like motif count, AT content, fraction outside coding sequence) and non-unique 50-mers masked. A gene is called ES when its longest run of ES windows has posterior ≥ 0.5 (P_ES_run). No essentiality labels are used."),
      tags$p(tags$b("Hybrid. "), "A random forest on SummaryTable gene features (hits, reads, freedom and neighbourhood indices, longest free interval) plus the HMM posteriors, trained on the curated labels (silver tier by default) with chromosome-grouped cross-validation; training genes receive out-of-fold scores."),
      tags$p(tags$b("Hybrid gene features. "), "Built automatically from the same parent libraries as the HMM, normalized with the chosen method (BGC by default; MidLC thins them to a common saturation), then combined and summarized (SummaryTable) (", tags$code("ytab hybrid-inputs"), "). SummaryTable needs whole numbers: normalized counts are rounded, and a site with insertions keeps at least 1, so site-based features match the data and read-based features carry the method's scaling. Results for the same method and parents are reused."),
      tags$p(tags$b(HTML("Validation (<em>C. glabrata</em>, Gale et al. 2020). ")), "Silver-label AUC: HMM 0.922, hybrid 0.947 (gold 0.974). The hybrid's gain over a random forest retrained without HMM features was small and not significant (+0.003, 95% CI −0.0004 to 0.006). The choice of normalization for the gene features changes the hybrid's AUC by at most 0.002. Genes called essential by both HMM and hybrid are the most reliable; HMM-only calls are candidates for species-specific essentiality or accessibility artefacts."),
      tags$p(tags$b("Inputs for other species. "), HTML("Use the haploid reference the reads were aligned to (for <em>C. albicans</em> the HapA FASTA); the run stops when more than 50% of windows are unmappable or more than 10% of insertion sites are off the reference. Use haploid strains: in a diploid, a heterozygous insertion in an essential gene is viable.")),
      ytab_reference("The HMM is adapted from the HMM of TRANSIT (DeJesus & Ioerger 2013, BMC Bioinformatics 14:303)."),
      tags$p("Command: ", tags$code("ytab essentiality"), " (which runs ", tags$code("ytab hermes-hmm"), " and the model table)."))))

essentiality_model_results_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Results",
    uiOutput("ess_result_selector"),
    uiOutput("ess_view_selector"),
    uiOutput("ess_coordinates_action"),
    uiOutput("ess_plot_selector"),
    conditionalPanel("input.ess_plot === 'track'", uiOutput("ess_chromosome_selector")),
    ytab_reference_projects_ui("ess", display_modes = FALSE),
    ytab_plot_downloads_ui("ess", c(genes = "Download gene table (TSV)"), capture = TRUE)),
  main = tagList(
    uiOutput("ess_summary"),
    ytab_plot_card("Essentiality", uiOutput("ess_plot_frame"), description = uiOutput("ess_plot_caption")),
    ytab_table_card("Genes", DT::DTOutput("ess_genes_table"))))

essentiality_module_ui <- function() navset_tab(id = "essentiality_tabs",
  nav_panel("Run", value = "run", essentiality_run_ui()),
  nav_panel("Results", value = "results", essentiality_model_results_ui()))

ess_read_windows <- function(run_dir) {
  path <- file.path(run_dir, "plots", "windows.tsv.gz")
  if (!file.exists(path)) return(data.frame())
  utils::read.delim(gzfile(path), stringsAsFactors = FALSE)
}

# Interactive chromosome track. Drawn natively with WebGL lines and a one-row state strip (a heatmap renders as one
# image): converting the stacked-area ggplot took ~6 s and drawing ~20,000 SVG area points took ~12 s per chromosome.
ess_track_plotly <- function(windows, chromosome, chrom_label = chromosome) {
  w <- windows[windows$chromosome == chromosome, , drop = FALSE]
  position <- format(w$start, big.mark = ",")
  states <- names(essentiality_state_labels)
  # the pipeline's most probable state per window ("unmappable" for masked windows)
  label <- as.character(w$state); best <- match(label, states, nomatch = 0L)
  hover <- sprintf("%s:%s<br>%d insertions<br>most probable: %s<br>P(ES) %.2f  P(GD) %.2f  P(NE) %.2f  P(GA) %.2f", chrom_label, position, w$insertions, label,
                   w$P_ES, w$P_GD, w$P_NE, w$P_GA)
  colours <- c("#d9d9d9", unname(essentiality_state_colours[states])); edges <- seq(0, 1, length.out = 6)
  scale <- do.call(rbind, lapply(seq_along(colours), function(k) rbind(c(edges[k], colours[k]), c(edges[k + 1], colours[k]))))
  insertions <- plotly::plot_ly(x = w$start, y = w$insertions, type = "scattergl", mode = "lines", line = list(shape = "hv", color = "#595959", width = 1),
                                text = hover, hoverinfo = "text", showlegend = FALSE) |> plotly::layout(yaxis = list(title = "Insertions"))
  strip <- plotly::plot_ly(x = w$start, y = "state", z = matrix(best, nrow = 1), type = "heatmap", zmin = -.5, zmax = 4.5, colorscale = scale,
                           showscale = FALSE, text = matrix(hover, nrow = 1), hoverinfo = "text") |> plotly::layout(yaxis = list(title = "", showticklabels = FALSE))
  lines <- plotly::plot_ly()
  for (state in c("ES", "GD")) lines <- plotly::add_trace(lines, x = w$start, y = w[[paste0("P_", state)]], type = "scattergl", mode = "lines",
    line = list(color = essentiality_state_colours[[state]], width = 1.2), name = sprintf("P(%s)", state), text = hover, hoverinfo = "text")
  lines <- plotly::layout(lines, yaxis = list(title = "Posterior", range = c(0, 1)))
  legend <- paste(sprintf("<span style='color:%s'>&#9632;</span> %s", c(colours[-1], colours[1]), c(essentiality_state_labels[states], "unmappable")), collapse = "  ")
  plotly::subplot(insertions, strip, lines, nrows = 3, shareX = TRUE, heights = c(.32, .1, .58), titleY = TRUE) |>
    plotly::layout(xaxis = list(title = sprintf("Position on %s (bp)", chrom_label)), legend = list(orientation = "h", y = -0.12),
                   annotations = list(list(text = paste("Most probable state:", legend), xref = "paper", yref = "paper", x = 0, y = 1.06,
                                           showarrow = FALSE, xanchor = "left", font = list(size = 11)))) |>
    ytab_plotly_style(native = TRUE)
}

# The same track read through the hybrid. The top panel is the libraries' insertions, identical in every view
# because it is a property of the data and not of a model. Below it the hybrid's calls and scores are drawn as
# segments spanning each gene, because the hybrid scores genes where the HMM scores windows -- the gene's span
# comes from the pipeline's gene table (start/end).
ess_gene_track_plotly <- function(windows, genes, chromosome, chrom_label = chromosome) {
  w <- windows[windows$chromosome == chromosome, , drop = FALSE]
  g <- genes[!is.na(genes$chromosome) & genes$chromosome == chromosome, , drop = FALSE]
  if (!all(c("start", "end") %in% names(g))) g$start <- NA_real_
  g <- g[!is.na(g$start) & !is.na(g$end) & !is.na(g$hybrid_score), , drop = FALSE]
  validate(need(nrow(g) > 0, paste(
    "This run's gene table has no coordinates, so the hybrid cannot be drawn against position.",
    "Use \u201cAdd gene coordinates\u201d above to write them from the project's annotation.")))
  g$call <- ifelse(as.logical(g$hybrid_essential), "Essential", "Not essential")
  g$tip <- sprintf("%s %s<br>%s:%s-%s<br>hybrid %.2f (%s)", g$ORF, ifelse(is.na(g$Name), "", g$Name),
                   chrom_label, format(g$start, big.mark = ","), format(g$end, big.mark = ","), g$hybrid_score, g$call)
  fills <- c(Essential = unname(essentiality_state_colours[["ES"]]),
             `Not essential` = unname(essentiality_state_colours[["NE"]]))
  # One trace per call rather than one per gene: thousands of traces will not render.
  segments <- function(d, y) { y <- rep_len(y, nrow(d))   # a constant row needs recycling to the gene count
    list(x = as.vector(rbind(d$start, d$end, NA)),
         y = as.vector(rbind(y, y, NA)),
         text = as.vector(rbind(d$tip, d$tip, NA))) }
  insertions <- plotly::plot_ly(x = w$start, y = w$insertions, type = "scattergl", mode = "lines",
                                line = list(shape = "hv", color = "#595959", width = 1), hoverinfo = "none", showlegend = FALSE) |>
    plotly::layout(yaxis = list(title = "Insertions"))
  strip <- plotly::plot_ly()
  for (call in names(fills)) {
    d <- g[g$call == call, , drop = FALSE]; if (!nrow(d)) next
    xy <- segments(d, 0)
    strip <- plotly::add_trace(strip, x = xy$x, y = xy$y, type = "scattergl", mode = "lines",
                               line = list(color = fills[[call]], width = 8), name = call, text = xy$text, hoverinfo = "text")
  }
  strip <- plotly::layout(strip, yaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE, range = c(-1, 1)))
  scores <- plotly::plot_ly()
  for (call in names(fills)) {
    d <- g[g$call == call, , drop = FALSE]; if (!nrow(d)) next
    xy <- segments(d, d$hybrid_score)
    scores <- plotly::add_trace(scores, x = xy$x, y = xy$y, type = "scattergl", mode = "lines",
                                line = list(color = fills[[call]], width = 2), showlegend = FALSE, text = xy$text, hoverinfo = "text")
  }
  scores <- plotly::layout(scores, yaxis = list(title = "Hybrid score", range = c(0, 1)),
                           shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = .5, y1 = .5,
                                              line = list(dash = "dash", width = 1, color = "#444"))))
  plotly::subplot(insertions, strip, scores, nrows = 3, shareX = TRUE, heights = c(.32, .1, .58), titleY = TRUE) |>
    plotly::layout(xaxis = list(title = sprintf("Position on %s (bp)", chrom_label)), legend = list(orientation = "h", y = -0.12)) |>
    ytab_plotly_style(native = TRUE)
}

# What each view paints along a chromosome. The HMM states a window; the hybrid calls a gene; the comparison
# labels a gene by whether the models agree. Returned as segments so one trace can carry a whole category.
ess_landscape_bands <- function(view, windows, genes, chromosome) {
  if (identical(view, "hmm")) {
    w <- windows[windows$chromosome == chromosome, , drop = FALSE]
    if (!nrow(w)) return(NULL)
    width <- suppressWarnings(diff(sort(w$start)[1:2])); if (!is.finite(width) || width <= 0) width <- 100
    # Neighbouring windows in the same state are one band, not one band each. Drawn window by window the
    # genome overview is ~127,000 segments and takes seconds to build; merged it is a few thousand and looks
    # the same, because adjacent segments of one colour are indistinguishable from the run they form.
    w <- w[order(w$start), , drop = FALSE]
    label <- as.character(w$state)
    run <- cumsum(c(TRUE, label[-1] != label[-length(label)] | diff(w$start) != width))
    first <- !duplicated(run); last <- !duplicated(run, fromLast = TRUE)
    return(list(data = data.frame(start = w$start[first], end = w$start[last] + width,
                                  category = label[first], stringsAsFactors = FALSE),
                colours = c(essentiality_state_colours, unmappable = "#d9d9d9"),
                labels = c(essentiality_state_labels, unmappable = "unmappable"), unit = "window"))
  }
  g <- genes[!is.na(genes$chromosome) & genes$chromosome == chromosome, , drop = FALSE]
  if (!all(c("start", "end") %in% names(g))) return(NULL)
  g <- g[!is.na(g$start) & !is.na(g$end), , drop = FALSE]
  if (identical(view, "hybrid")) {
    g <- g[!is.na(g$hybrid_score), , drop = FALSE]; if (!nrow(g)) return(NULL)
    return(list(data = data.frame(start = g$start, end = g$end, name = g$ORF,
                                  category = ifelse(as.logical(g$hybrid_essential), "Essential", "Not essential"), stringsAsFactors = FALSE),
                colours = c(Essential = unname(essentiality_state_colours[["ES"]]),
                            `Not essential` = unname(essentiality_state_colours[["NE"]])),
                labels = c(Essential = "Essential", `Not essential` = "Not essential"), unit = "gene"))
  }
  g <- g[!is.na(g$agreement) & nzchar(as.character(g$agreement)), , drop = FALSE]; if (!nrow(g)) return(NULL)
  list(data = data.frame(start = g$start, end = g$end, name = g$ORF, category = as.character(g$agreement), stringsAsFactors = FALSE),
       colours = essentiality_agreement_colours, labels = essentiality_agreement_labels, unit = "gene")
}

# One row per chromosome: the genome's landscape at a glance, not a full reading of any one chromosome.
ess_genome_overview_plotly <- function(windows, genes, chromosomes, labels, view = "hmm") {
  bands <- lapply(chromosomes, function(chrom) ess_landscape_bands(view, windows, genes, chrom))
  names(bands) <- chromosomes
  keep <- !vapply(bands, is.null, logical(1))
  validate(need(any(keep), "Nothing to draw for this model yet."))
  bands <- bands[keep]; rows <- labels[keep]
  palette <- bands[[1]]$colours; legend_labels <- bands[[1]]$labels
  order <- rev(seq_along(bands))            # first chromosome at the top
  fig <- plotly::plot_ly()
  seen <- character(0)
  for (category in names(palette)) {
    xs <- numeric(0); ys <- numeric(0); tips <- character(0)
    for (i in seq_along(bands)) {
      d <- bands[[i]]$data; d <- d[d$category == category, , drop = FALSE]; if (!nrow(d)) next
      tip <- sprintf("%s<br>%s-%s<br>%s", rows[[i]], format(d$start, big.mark = ","), format(d$end, big.mark = ","),
                     legend_labels[[category]] %||% category)
      xs <- c(xs, as.vector(rbind(d$start, d$end, NA)))
      ys <- c(ys, as.vector(rbind(rep(order[i], nrow(d)), rep(order[i], nrow(d)), NA)))
      tips <- c(tips, as.vector(rbind(tip, tip, NA)))
    }
    if (!length(xs)) next
    seen <- c(seen, category)
    fig <- plotly::add_trace(fig, x = xs, y = ys, type = "scattergl", mode = "lines",
                             line = list(color = palette[[category]], width = 14),
                             name = legend_labels[[category]] %||% category, text = tips, hoverinfo = "text")
  }
  validate(need(length(seen) > 0, "Nothing to draw for this model yet."))
  plotly::layout(fig,
    xaxis = list(title = "Position (bp)", zeroline = FALSE),
    yaxis = list(title = "", tickmode = "array", tickvals = order, ticktext = rows,
                 range = c(0.3, length(bands) + 0.7), zeroline = FALSE, showgrid = FALSE),
    legend = list(orientation = "h", y = -0.18)) |> ytab_plotly_style(native = TRUE)
}

# Several chromosomes together: each keeps its insertion profile, because comparing those landscapes is the
# reason to put them side by side; the band below each says what the model made of it.
ess_multi_chromosome_plotly <- function(windows, genes, chromosomes, labels, view = "hmm") {
  panels <- list(); heights <- numeric(0)
  for (i in seq_along(chromosomes)) {
    chrom <- chromosomes[[i]]
    w <- windows[windows$chromosome == chrom, , drop = FALSE]
    band <- ess_landscape_bands(view, windows, genes, chrom)
    profile <- plotly::plot_ly(x = w$start, y = w$insertions, type = "scattergl", mode = "lines",
                               line = list(shape = "hv", color = "#595959", width = 1), hoverinfo = "none", showlegend = FALSE) |>
      plotly::layout(yaxis = list(title = labels[[i]]))
    strip <- plotly::plot_ly()
    if (!is.null(band)) for (category in names(band$colours)) {
      d <- band$data[band$data$category == category, , drop = FALSE]; if (!nrow(d)) next
      tip <- sprintf("%s<br>%s-%s<br>%s", labels[[i]], format(d$start, big.mark = ","), format(d$end, big.mark = ","),
                     band$labels[[category]] %||% category)
      strip <- plotly::add_trace(strip, x = as.vector(rbind(d$start, d$end, NA)),
                                 y = as.vector(rbind(rep(0, nrow(d)), rep(0, nrow(d)), NA)),
                                 type = "scattergl", mode = "lines", line = list(color = band$colours[[category]], width = 10),
                                 name = band$labels[[category]] %||% category,
                                 text = as.vector(rbind(tip, tip, NA)), hoverinfo = "text",
                                 showlegend = i == 1L, legendgroup = category)
    }
    strip <- plotly::layout(strip, yaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE, range = c(-1, 1)))
    panels <- c(panels, list(profile, strip)); heights <- c(heights, .78, .22)
  }
  plotly::subplot(panels, nrows = length(panels), shareX = TRUE, titleY = TRUE,
                  heights = heights/sum(heights), margin = .012) |>
    plotly::layout(xaxis = list(title = "Position (bp)"), legend = list(orientation = "h", y = -0.1)) |>
    ytab_plotly_style(native = TRUE)
}

# The downloadable form of the chromosome track. It takes the same inputs the on-screen track does and
# follows the same three modes, so a download is a picture of what is being looked at rather than of a
# default that was chosen for it.
# Model colours: the label-free HMM in teal, the hybrid in brick red, the current random forest in grey.
essentiality_model_colours <- c(`HMM (label-free)` = YTAB_FOCAL, Hybrid = YTAB_ACCENT, `Random forest` = YTAB_COMPARATOR)

# `view` is which model the results are read through: "hmm", "hybrid", or "compare" for the two side by side.
# The comparison used to be unavoidable wherever a hybrid existed, which presented the two models as co-equal.
ess_ggplot <- function(choice, genes, windows = NULL, chromosome = NULL, manifest = list(), run_dir = NULL,
                       reference = data.frame(), project_label = "This dataset", view = "hmm") {
  base <- ggplot2::theme_minimal(base_size = 12)
  if (identical(choice, "target")) {
    tm <- ytab_read_analysis_tsv(file.path(run_dir %||% "", "plots", "target_model.tsv"))
    validate(need(nrow(tm) > 0, "This run has no insertion-target summary. Run essentiality again (or `ytab essentiality --refresh-plots RUN`)."))
    used <- isTRUE(as.logical(tm$used_in_fit[[1L]]))
    tm$tip <- sprintf("offset %+d<br>%.2f bits<br>A %.2f  C %.2f  G %.2f  T %.2f", tm$offset, tm$information_bits, tm$A, tm$C, tm$G, tm$T)
    return(ggplot2::ggplot(tm, ggplot2::aes(offset, information_bits, text = tip)) + ggplot2::geom_col(fill = YTAB_FOCAL, width = .8) +
      # A base is drawn only where the position is informative relative to the rest of the profile. A fixed
      # 0.05-bit cut labelled 10 of 28 positions on Hermes -- two real (T at +1, A at +6, ~0.58 bits) and eight
      # at 0.05-0.11, which read as a motif that is not there. A share of the profile's own maximum keeps this
      # right for other transposons, where the informative positions differ in number (TTAA has four).
      # Clear of the bar: vjust alone left the base sitting on the bar's top edge. The gap is a share of the
      # tallest bar, so it holds whatever the information range of the transposon turns out to be.
      ggplot2::geom_text(data = tm[tm$information_bits >= max(.1, .25 * max(tm$information_bits, na.rm = TRUE)), , drop = FALSE],
                         ggplot2::aes(label = preferred), vjust = 0, nudge_y = .04 * max(tm$information_bits, na.rm = TRUE),
                         fontface = "bold", size = 4.5) +
      # Every offset is labelled. Ticking every second one left the odd positions bare -- including +1, which
      # for Hermes is one of the two positions that carry the motif.
      ggplot2::scale_x_continuous(breaks = sort(unique(tm$offset))) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .25))) +
      ggplot2::labs(x = "Offset from the insertion site (0 = first base of the target-site duplication)", y = "Information (bits)",
                    title = if (used) "Learned from these libraries and used in the HMM fit" else "Learned from these libraries (diagnostic; the fit used the Hermes motif)") + base)
  }
  if (identical(choice, "roc")) {
    roc <- ytab_read_analysis_tsv(file.path(run_dir %||% "", "plots", "roc.tsv"))
    validate(need(nrow(roc) > 0, "ROC curves need curated labels: run essentiality with the hybrid model (labels), or refresh an earlier run's plots."))
    # roc.tsv carries a curve per model, so a single model's view shows its own rather than borrowing the pair.
    if (identical(view, "hmm")) roc <- roc[grepl("HMM", roc$model, fixed = TRUE), , drop = FALSE]
    if (identical(view, "hybrid")) roc <- roc[grepl("Hybrid", roc$model, fixed = TRUE), , drop = FALSE]
    validate(need(nrow(roc) > 0, "This run has no ROC curve for the chosen model."))
    roc$tier <- factor(ifelse(roc$tier == "all", "All labelled genes", paste(tools::toTitleCase(roc$tier), "labels")))
    roc$tier <- factor(roc$tier, unique(c("All labelled genes", levels(roc$tier))))
    auc <- unique(roc[, c("model", "tier", "genes", "essential", "auc")])
    auc$label <- sprintf("%s %.3f", sub(" \\(label-free\\)", "", auc$model), auc$auc)
    # labels centred inside the panel: plotly ignores hjust, so right-aligned labels were clipped at the panel edge
    auc$y <- .05 + .08 * (stats::ave(seq_len(nrow(auc)), auc$tier, FUN = seq_along) - 1)
    roc$tip <- sprintf("%s<br>%s<br>FPR %.3f, TPR %.3f<br>AUC %.3f (%s genes, %s essential)", roc$model, roc$tier, roc$fpr, roc$tpr, roc$auc,
                       format(roc$genes, big.mark = ","), format(roc$essential, big.mark = ","))
    return(ggplot2::ggplot(roc, ggplot2::aes(fpr, tpr, colour = model, group = model, text = tip)) +
      ggplot2::geom_abline(linetype = "dotted", colour = "grey50") + ggplot2::geom_step(linewidth = .9) +
      ggplot2::geom_text(data = auc, ggplot2::aes(x = .7, y = y, label = label, colour = model), fontface = "bold", size = 3.6, inherit.aes = FALSE, show.legend = FALSE) +
      ggplot2::facet_wrap(~tier, nrow = 1) + ggplot2::scale_colour_manual(values = essentiality_model_colours) + ggplot2::coord_equal() +
      ggplot2::scale_x_continuous(breaks = c(0, .5, 1)) + ggplot2::scale_y_continuous(breaks = c(0, .5, 1)) +
      ggplot2::labs(x = "False positive rate", y = "True positive rate", colour = NULL) + base + ggplot2::theme(legend.position = "bottom"))
  }
  if (identical(choice, "track")) {
    w <- windows[windows$chromosome == chromosome, , drop = FALSE]
    posterior <- do.call(rbind, lapply(names(essentiality_state_labels), function(s)
      data.frame(start = w$start, state = s, value = w[[paste0("P_", s)]])))
    posterior$state <- factor(posterior$state, levels = names(essentiality_state_labels))
    counts <- data.frame(start = w$start, insertions = w$insertions, mappable = w$mappable == 1)
    top <- ggplot2::ggplot(counts, ggplot2::aes(start, insertions, text = sprintf("%s:%s<br>%d insertions%s", chromosome, format(start, big.mark = ","), insertions,
                                                                                  ifelse(mappable, "", " (unmappable window)")))) +
      ggplot2::geom_col(ggplot2::aes(fill = mappable), width = diff(w$start[1:2]) %||% 100) +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "grey35", `FALSE` = "grey80"), guide = "none") +
      ggplot2::scale_y_sqrt() + ggplot2::labs(x = NULL, y = "Insertions (all libraries)") + base
    bottom <- ggplot2::ggplot(posterior, ggplot2::aes(start, value, fill = state, group = state, text = sprintf("%s:%s<br>%s %.2f", chromosome, format(start, big.mark = ","), state, value))) +
      ggplot2::geom_area(position = "stack") + ggplot2::scale_fill_manual(values = essentiality_state_colours, labels = essentiality_state_labels) +
      ggplot2::labs(x = sprintf("Position on %s (bp)", chromosome), y = "HMM state posterior", fill = NULL) + base
    return(list(top = top, bottom = bottom))
  }
  if (identical(choice, "scores")) {
    if (identical(view, "compare") && "hybrid_score" %in% names(genes) && any(!is.na(genes$hybrid_score))) {
      d <- genes[!is.na(genes$hybrid_score), , drop = FALSE]
      return(ggplot2::ggplot(d, ggplot2::aes(hmm_score, hybrid_score, colour = agreement,
          text = sprintf("%s %s<br>HMM %.2f (%s)<br>hybrid %.2f<br>%s", ORF, ifelse(is.na(Name), "", Name), hmm_score, hmm_call, hybrid_score, agreement))) +
        ggplot2::geom_point(alpha = .6, size = 1.4) + ggplot2::geom_hline(yintercept = .5, linetype = "dashed") + ggplot2::geom_vline(xintercept = .5, linetype = "dashed") +
        ggplot2::scale_colour_manual(values = essentiality_agreement_colours, labels = essentiality_agreement_labels) +
        ggplot2::labs(x = "HMM score (P_ES_run)", y = "Hybrid score", colour = "Agreement") + base)
    }
    # One model, either model: the same histogram of the score coloured by that model's call. Only the column
    # and its palette change, so the two views are read the same way rather than being different plot types.
    hybrid_view <- identical(view, "hybrid") && "hybrid_score" %in% names(genes) && any(!is.na(genes$hybrid_score))
    score <- if (hybrid_view) "hybrid_score" else "hmm_score"
    x_lab <- if (hybrid_view) "Hybrid score" else "HMM score (P_ES_run)"
    d <- genes[!is.na(genes[[score]]), , drop = FALSE]
    d$.score <- as.numeric(d[[score]])
    if (hybrid_view) {
      d$.call <- ifelse(as.logical(d$hybrid_essential), "Essential", "Not essential")
      fills <- c(Essential = unname(essentiality_state_colours[["ES"]]),
                 `Not essential` = unname(essentiality_state_colours[["NE"]]))
      legend <- "Hybrid call"
    } else {
      d$.call <- d$hmm_call
      fills <- c(essentiality_state_colours, uncertain = "grey70", `N/A` = "grey90")
      legend <- "HMM call"
    }
    if (is.data.frame(reference) && nrow(reference) && !hybrid_view) {
      # Projects score different numbers of genes, so each is its own density: the shape is what carries across.
      long <- rbind(data.frame(.project = project_label, hmm_score = d$.score, stringsAsFactors = FALSE),
                    data.frame(.project = reference$.project, hmm_score = reference$hmm_score, stringsAsFactors = FALSE))
      long <- long[is.finite(long$hmm_score), , drop = FALSE]
      return(ggplot2::ggplot(long, ggplot2::aes(hmm_score, colour = .project)) + ggplot2::geom_density(linewidth = .9) +
        ggplot2::scale_colour_manual(values = ytab_reference_project_fills(unique(long$.project), project_label)) +
        ggplot2::labs(x = x_lab, y = "Density", colour = NULL) + base)
    }
    return(ggplot2::ggplot(d, ggplot2::aes(.score, fill = .call)) + ggplot2::geom_histogram(bins = 50) +
      ggplot2::scale_fill_manual(values = fills) +
      ggplot2::labs(x = x_lab, y = "Genes", fill = legend) + base)
  }
  # counts as written by the pipeline (manifest calls / agreement)
  counts <- function(x, model) if (length(x)) data.frame(call = names(x), Freq = as.integer(unlist(x)), model = model, stringsAsFactors = FALSE)
  # The manifest records the HMM's four states and the two models' agreement, not the hybrid's own calls: the
  # hybrid is a binary classifier, so its counts are the agreement read down the hybrid's side.
  agree <- manifest$agreement
  hybrid_calls <- if (length(agree))
    list(`Essential` = (agree[["consensus essential"]] %||% 0) + (agree[["hybrid only"]] %||% 0),
         `Not essential` = (agree[["HMM only"]] %||% 0) + (agree[["not essential"]] %||% 0))
  calls <- switch(view,
    hybrid  = counts(hybrid_calls, "Hybrid"),
    compare = counts(agree, "HMM vs hybrid"),
    counts(manifest$calls, "HMM"))
  if (is.null(calls) || !nrow(calls)) return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "Call counts are not in this run's manifest; run it again.") + ggplot2::theme_void())
  calls$label <- ifelse(calls$call %in% names(essentiality_state_labels), essentiality_state_labels[calls$call],
                        ifelse(calls$call %in% names(essentiality_agreement_labels), essentiality_agreement_labels[calls$call], calls$call))
  if (is.data.frame(reference) && nrow(reference)) {
    # Projects call different numbers of genes, so calls are shown as a share of each project's own genes.
    hmm <- calls[calls$model == "HMM", , drop = FALSE]
    share <- function(x) if (sum(x) > 0) 100 * x / sum(x) else x
    by_project <- split(reference$hmm_call, reference$.project)
    long <- rbind(
      data.frame(.project = project_label, call = hmm$call, value = share(hmm$Freq), stringsAsFactors = FALSE),
      do.call(rbind, lapply(names(by_project), function(name) {
        counted <- table(by_project[[name]])
        data.frame(.project = name, call = names(counted), value = share(as.numeric(counted)), stringsAsFactors = FALSE)
      })))
    long$.project <- factor(long$.project, levels = c(project_label, setdiff(unique(long$.project), project_label)))
    long$label <- ifelse(long$call %in% names(essentiality_state_labels), essentiality_state_labels[long$call], long$call)
    return(ggplot2::ggplot(long, ggplot2::aes(stats::reorder(call, -value), value, fill = .project,
        text = sprintf("%s<br>%s: %.1f%% of its genes", .project, label, value))) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = .8), width = .75) +
      ggplot2::scale_fill_manual(values = ytab_reference_project_fills(unique(long$.project), project_label)) +
      ggplot2::labs(x = NULL, y = "Genes of the dataset (%)", fill = NULL) + base +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)))
  }
  # The hybrid's own two calls are not in either map: it is a binary classifier, so it reuses the essential and
  # non-essential colours the HMM states already use, keeping one colour per meaning across the views.
  fills <- c(essentiality_state_colours, essentiality_agreement_colours, uncertain = "grey70",
             `Essential` = unname(essentiality_state_colours[["ES"]]),
             `Not essential` = unname(essentiality_state_colours[["NE"]]))
  ggplot2::ggplot(calls, ggplot2::aes(stats::reorder(call, -Freq), Freq, fill = call, text = sprintf("%s: %s genes", label, format(Freq, big.mark = ",")))) +
    ggplot2::geom_col(show.legend = FALSE) + ggplot2::scale_fill_manual(values = fills) +
    ggplot2::scale_x_discrete(labels = function(x) ifelse(x %in% names(essentiality_agreement_labels), essentiality_agreement_labels[x], x)) +
    ggplot2::facet_wrap(~model, scales = "free_x") +
    ggplot2::labs(x = NULL, y = "Genes") + base + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
}

essentiality_module_server <- function(input, output, session, repo_root, active, active_project_path, jobs, python_bin,
                                       projects = reactive(data.frame()), reference_state = NULL) {
  if (is.null(reference_state)) reference_state <- ytab_reference_state()
  ytab_reference_server("ess", input, output, session, repo_root, active, projects, reference_state, display_modes = FALSE)
  reference_essentiality <- reactive(ytab_reference_essentiality(ytab_reference_metas(reference_state(), repo_root)))
  runs_dir <- reactive(file.path(active()$project_root, "essentiality"))
  runs <- reactivePoll(3000, session,
    checkFunc = function() { x <- ytab_analysis_runs(tryCatch(runs_dir(), error = function(e) "")); paste(x$path, x$modified, collapse = "|") },
    valueFunc = function() { x <- ytab_analysis_runs(tryCatch(runs_dir(), error = function(e) "")); x[file.exists(file.path(x$path, "models", "genes.tsv")), , drop = FALSE] })
  inventory <- reactive(essentiality_parent_inventory(active()))
  design_support <- reactive(ytab_design_support(python_bin, repo_root, tryCatch(active_project_path(), error = function(e) "")))
  # Whether this project's libraries are dense enough for an essentiality call to be worth making, read from the
  # newest QC run (`ytab qc` writes it; the app does not recompute density). Without this the tab invites the user
  # to fit an HMM on a library that cannot support one, and the calls come back looking as ordinary as any other.
  readiness <- reactive({
    qc <- ytab_analysis_runs(tryCatch(file.path(active()$project_root, "qc"), error = function(e) ""))
    if (!nrow(qc)) return(data.frame())
    ytab_read_analysis_tsv(file.path(qc$path[[1L]], "essentiality_readiness.tsv"))
  })
  readiness_banner <- reactive({
    r <- readiness()
    if (!nrow(r) || !"sites_per_coding_kb" %in% names(r)) return(NULL)
    r <- r[!is.na(r$sites_per_coding_kb) & r$n_libraries > 0, , drop = FALSE]
    short <- r[!as.logical(r$supports_hybrid), , drop = FALSE]
    if (!nrow(short)) return(NULL)
    ytab_warning_banner("These libraries are sparser than an essentiality call needs",
      paste0(paste(sprintf("Strain %s: %.1f insertion sites per coding kb pooled over %d %s librar%s.",
                           short$background, short$sites_per_coding_kb, short$n_libraries, short$role,
                           ifelse(short$n_libraries == 1, "y", "ies")), collapse = " "),
             sprintf(" The label-free HMM asks for about %g and the hybrid for about %g, reached by pooling independent inductions. ",
                     r$hmm_target[[1L]], r$hybrid_target[[1L]]),
             "A run on libraries this sparse will still produce calls, but they mostly reflect where the transposon happened to land. See the QC tab for the per-library figures."))
  })
  # Normalizations the hybrid's gene features can use: the methods of the project's normalization runs (MidLC
  # common saturation among them). BGC first when available (the app's suggested normalization).
  hybrid_normalizations <- reactive({
    root <- active()$project_root
    methods <- unique(unlist(lapply(ytab_analysis_runs(file.path(root, "normalization"))$path, function(run) ytab_read_manifest(run)$methods %||% character())))
    methods <- setdiff(as.character(methods), "nonorm"); methods <- c(intersect("betageom", methods), setdiff(methods, "betageom"))
    choices <- stats::setNames(methods, as.character(norm_method_label(methods)))
    choices
  })
  label_files <- reactive(list.files(file.path(repo_root, "resources", "species"), pattern = "^labels\\.tsv$", recursive = TRUE, full.names = TRUE))

  # The hybrid trains on curated essentiality labels, which exist for a handful of species. Offering it where
  # there are none produces a run that fails at its last step, so the option states the reason instead.
  output$ess_model_selector <- renderUI({
    have_labels <- length(label_files()) > 0L
    choices <- c("HMM" = "hmm")
    if (have_labels) choices <- c(choices, "HMM + hybrid" = "hybrid")
    tags$div(class = "ytab-method-switch",
      radioButtons("ess_model", "Model", choices, "hmm", inline = TRUE),
      if (!have_labels)
        tags$p(class = "text-muted",
               "The hybrid model needs curated essentiality labels; this installation has none. Add them under resources/species/<species>/essentiality_labels/labels.tsv."))
  })

  output$ess_library_selector <- renderUI({
    x <- inventory(); x <- x[x$Included & !x$Explicit_treated, , drop = FALSE]
    choices <- stats::setNames(x$Sample, sprintf("%s (%s, pool %s)", ytab_sample_labels(x$Sample, active()), x$Background, x$Pool))
    selectizeInput("ess_libraries", "Libraries", choices, selected = essentiality_default_libraries(inventory()), multiple = TRUE)
  })
  output$ess_hmm_reuse_toggle <- renderUI({
    if (!nrow(runs())) return(NULL)
    checkboxInput("ess_reuse_fit", "Reuse an earlier HMM fit (skips the slow fitting step)", FALSE)
  })
  output$ess_hmm_reuse_selector <- renderUI({
    x <- runs(); if (!nrow(x)) return(NULL)
    selectInput("ess_hmm_reuse", "Earlier fit", ytab_run_choices(x))
  })
  output$ess_hybrid_inputs <- renderUI({
    norms <- hybrid_normalizations(); labels <- label_files()
    tagList(
      selectInput("ess_labels", "Curated labels", c(stats::setNames(labels, sub("/essentiality_labels/labels.tsv$", "", sub(".*/resources/species/", "", labels))), "Other file" = "other"),
                  selected = labels[grepl(paste0("/species/", tolower(active()$species %||% ""), "/"), labels, fixed = TRUE)][1] %||% labels[1]),
      conditionalPanel("input.ess_labels === 'other'", textInput("ess_labels_path", "Label file (TSV: gene_id, label, tier)")),
      if (length(norms)) selectInput("ess_normalization", "Normalization for the gene features", norms)
      else tags$p(class = "text-warning", "No normalization in this dataset yet. Run Normalization (or MidLC) first."))
  })

  observeEvent(input$ess_run, {
    if (jobs$job_is_running()) { showNotification("Wait for the current job to finish.", type = "warning"); return() }
    reuse <- isTRUE(input$ess_reuse_fit)
    hmm_args <- if (reuse) {
      run <- input$ess_hmm_reuse %||% ""
      if (!nzchar(run)) { showNotification("Choose an earlier fit.", type = "error"); return() }
      # A manifest without hmm_run would drop out of c() and leave a bare "--hmm-run", which then swallows the
      # next flag as its value and starts a job with the wrong arguments rather than failing.
      fit <- ytab_read_manifest(run)$hmm_run %||% ""
      if (!nzchar(fit)) { showNotification("That result does not record which HMM fit it used.", type = "error"); return() }
      ytab_flag("--hmm-run", fit, required = TRUE)
    } else {
      libraries <- input$ess_libraries %||% character()
      if (!length(libraries)) { showNotification("Choose at least one library.", type = "error"); return() }
      c("--libraries", libraries, "--window", as.character(input$ess_window %||% 100), "--min-reads", as.character(input$ess_min_reads %||% 1))
    }
    hybrid_args <- character()
    if (identical(input$ess_model, "hybrid")) {
      labels <- if (identical(input$ess_labels, "other")) trimws(input$ess_labels_path %||% "") else input$ess_labels %||% ""
      if (!nzchar(input$ess_normalization %||% "") || !nzchar(labels) || !file.exists(labels)) {
        showNotification("The hybrid needs a label file and a normalization (Normalization tab).", type = "error"); return()
      }
      hybrid_args <- c(ytab_flag("--labels", labels, required = TRUE), ytab_flag("--normalization", input$ess_normalization, required = TRUE))
    }
    args <- c(file.path(repo_root, "scripts", "local", "ytab_essentiality.py"), "--project-config", active_project_path(), hmm_args, hybrid_args,
              "--output", ytab_new_run_dir(runs_dir(), "essentiality"))
    tryCatch({
      jobs$start_job(python_bin(), args, active_project_path(), "essentiality_models", wd = repo_root)
      showNotification(if (reuse) "Building the model table." else "Fitting the HMM; this takes several minutes.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error"))
  }, ignoreInit = TRUE)

  output$ess_run_status <- renderUI({
    x <- runs()
    gate <- ytab_design_gate(design_support(), "Essentiality", "Essentiality is not available for this dataset")
    if (!nrow(x) && !is.null(gate)) return(gate)
    if (!nrow(x)) return(tagList(readiness_banner(), ytab_design_merge_warning(design_support()),
      ytab_empty_state("No essentiality result yet", "The HMM needs only the hit files of untreated parent libraries. Choose the libraries and run.")))
    m <- ytab_read_manifest(x$path[[1L]])
    models <- if (isTRUE(m$hybrid)) "HMM and hybrid" else "HMM"
    tagList(readiness_banner(),
      ytab_result_card("Latest result", tags$p(sprintf("%s: %s on %s.", names(ytab_run_choices(x[1L, , drop = FALSE])), models, paste(ytab_sample_labels(unlist(m$libraries), active()), collapse = ", "))),
                       tags$p("See the Results tab.")))
  })

  selected_run <- reactive({ req(nrow(runs()) > 0); path <- input$ess_result %||% runs()$path[[1L]]; req(path %in% runs()$path); path })
  manifest <- reactive(ytab_read_manifest(selected_run()))
  genes <- reactive({ g <- ytab_read_analysis_tsv(file.path(selected_run(), "models", "genes.tsv"))
    if (nrow(g) && "ORF" %in% names(g)) g$Name <- ytab_display_gene_names(g$ORF, g$Name, repo_root)   # display names by the species rule
    g })
  windows <- reactive(ess_read_windows(selected_run()))

  output$ess_result_selector <- renderUI({
    x <- runs(); if (!nrow(x)) return(tags$p(class = "text-muted", "No essentiality result yet."))
    selectInput("ess_result", "Result", ytab_run_choices(x))
  })
  # The ROC is drawn against the curated labels, which only a hybrid run has. Offering it for an HMM-only
  # result gave a plot with nothing in it.
  # Which model the results are read through. The HMM is the default: it is the Run tab's default, it needs no
  # labels, and for most species it is the only model there is. The hybrid and the comparison are offered only
  # where the run produced one.
  output$ess_view_selector <- renderUI({
    if (!isTRUE(manifest()$hybrid)) return(NULL)
    tags$div(class = "ytab-method-switch",
      radioButtons("ess_view", "View", c("HMM" = "hmm", "Hybrid" = "hybrid", "Compare models" = "compare"), "hmm", inline = TRUE))
  })
  ess_view <- reactive(if (isTRUE(manifest()$hybrid)) input$ess_view %||% "hmm" else "hmm")

  # Results produced before the pipeline recorded gene coordinates cannot draw a gene against position. The
  # fix is a pipeline step, so the app runs it: the user presses a button, not a command line.
  ess_has_coordinates <- reactive({
    g <- genes()
    nrow(g) > 0 && all(c("start", "end") %in% names(g)) && any(!is.na(g$start))
  })
  output$ess_coordinates_action <- renderUI({
    if (!nrow(runs()) || ess_has_coordinates()) return(NULL)
    tagList(
      tags$p(class = "text-muted", "This result was made before gene coordinates were recorded, so the chromosome track cannot show genes."),
      actionButton("ess_add_coordinates", "Add gene coordinates", class = "btn-default"))
  })
  observeEvent(input$ess_add_coordinates, {
    if (jobs$job_is_running()) { showNotification("Wait for the current job to finish.", type = "warning"); return() }
    run <- selected_run(); req(nzchar(run %||% ""))
    args <- c(file.path(repo_root, "scripts", "local", "ytab_essentiality.py"),
              "--project-config", active_project_path(), "--refresh-plots", run)
    tryCatch({
      jobs$start_job(python_bin(), args, active_project_path(), "essentiality_models", wd = repo_root)
      showNotification("Writing gene coordinates from the project's annotation.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error"))
  })

  output$ess_plot_selector <- renderUI({
    view <- ess_view()
    # Every view offers the same plots, each drawn for the model being viewed: the HMM's track is per window,
    # the hybrid's per gene, and both sit over the same insertion profile.
    choices <- c("Chromosome track" = "track", "Model scores" = "scores", "Calls" = "calls",
                 "Insertion-target preference" = "target")
    if (isTRUE(manifest()$hybrid)) choices <- c(choices, "ROC against curated labels" = "roc")
    radioButtons("ess_plot", "Plot", choices, selected = if ((input$ess_plot %||% "") %in% choices) input$ess_plot else choices[[1L]])
  })
  # Nothing selected is the useful default: the whole genome at a glance. One chromosome gives the detailed
  # track, several put their insertion landscapes side by side, which is the comparison the single-select
  # could never express.
  output$ess_chromosome_selector <- renderUI({ w <- windows(); req(nrow(w) > 0); chroms <- unique(w$chromosome)
    tagList(
      selectizeInput("ess_chromosome", "Chromosome", stats::setNames(chroms, qc_plot_chromosome_display(chroms, active())),
                     selected = character(0), multiple = TRUE,
                     options = list(plugins = list("remove_button"), placeholder = "All chromosomes (genome overview)")),
      helpText("Leave empty for the whole genome. Choose one for the detailed track, or several to compare them."))
  })
  output$ess_summary <- renderUI({
    if (!nrow(runs())) return(ytab_empty_state("No essentiality result yet", "Run essentiality in the Run tab."))
    m <- manifest(); req(m); calls <- m$hmm$gene_calls
    cards <- lapply(names(essentiality_state_labels), function(s) ytab_metric_card(essentiality_state_labels[[s]], format(calls[[s]] %||% 0, big.mark = ","), "genes (HMM)"))
    view <- ess_view()
    if (identical(view, "compare"))
      cards <- lapply(c("consensus essential", "HMM only", "hybrid only"), function(a)
        ytab_metric_card(essentiality_agreement_labels[[a]], format(m$agreement[[a]] %||% 0, big.mark = ","), "genes (HMM vs hybrid)"))
    if (identical(view, "hybrid")) {
      a <- m$agreement
      cards <- list(
        ytab_metric_card("Essential", format((a[["consensus essential"]] %||% 0) + (a[["hybrid only"]] %||% 0), big.mark = ","), "genes (hybrid)"),
        ytab_metric_card("Not essential", format((a[["HMM only"]] %||% 0) + (a[["not essential"]] %||% 0), big.mark = ","), "genes (hybrid)"))
    }
    notes <- unlist(m$hmm$notes %||% list())
    tagList(tags$div(class = "ytab-stat-grid", cards),
      tags$p(class = "text-muted", sprintf("Libraries: %s. Unmappable windows: %.1f%%. %s", paste(ytab_sample_labels(unlist(m$libraries), active()), collapse = ", "), 100 * m$hmm$masked_fraction,
                                           if (isTRUE(m$hmm$converged)) "The HMM converged." else "The HMM did not converge; see the note.")),
      if (length(notes)) ytab_warning_banner("Notes from the HMM", paste(notes, collapse = " ")))
  })
  output$ess_plot_caption <- renderUI(switch(input$ess_plot %||% "track",
    track = if (!length(intersect(input$ess_chromosome %||% character(0), unique(windows()$chromosome))))
      switch(ess_view(),
        hybrid  = "Every chromosome at a glance, each gene coloured by the hybrid's call. Choose a chromosome for its detailed track, or several to compare them.",
        compare = "Every chromosome at a glance, each gene coloured by whether the two models agree. Choose a chromosome for its detailed track, or several to compare them.",
        "Every chromosome at a glance, each window coloured by its most probable HMM state. Choose a chromosome for its detailed track, or several to compare them.")
    else if (length(intersect(input$ess_chromosome %||% character(0), unique(windows()$chromosome))) > 1L)
      "Insertions per window for each chosen chromosome, with the model's calls banded underneath, so the landscapes can be compared."
    else if (identical(ess_view(), "hybrid"))
      "Insertions per window (top), each gene's hybrid call (strip) and its hybrid score (bottom), drawn across the gene's span. The dashed line is the 0.5 threshold."
      else "Insertions per window (top), the most probable HMM state of each window (strip; grey = unmappable) and the posterior probabilities of essential and growth-defect states (bottom). Zoom with the mouse; double-click to reset. The PNG/PDF download shows all four state posteriors.",
    # The caption follows the view, because the view decides which plot was drawn: describing the HMM-vs-hybrid
    # scatter under a single model's density is how this read before the views existed.
    scores = switch(ess_view(),
      hybrid  = "Distribution of the hybrid score by call.",
      compare = "Each point is a gene. Dashed lines are the 0.5 thresholds used for both models.",
      "Distribution of the HMM score (P_ES_run) by call."),
    calls = switch(ess_view(),
      hybrid  = "Number of genes the hybrid calls essential and not essential.",
      compare = "How the two models' calls agree: genes both call essential, and those only one of them does.",
      "Number of genes per HMM call."),
    target = "Information (bits) at each position around the insertion site, from the libraries' distinct insertion sites against the local base composition, with the preferred base above informative positions. No labels or prior motif are used; Hermes shows T at +1 and A at +6 (nTnnnnAn).",
    roc = if (identical(ess_view(), "compare")) "Each model's score against the curated labels, per label tier and all together, with the area under the curve (AUC). Only genes a model was not fitted to: the HMM uses no labels, and the hybrid's training genes carry out-of-fold scores."
          else "This model's score against the curated labels, per label tier and all together, with the area under the curve (AUC). Only genes the model was not fitted to."))
  # The track grows with what is on it: a glance at the genome is short, stacked chromosomes are not.
  ess_plot_height <- reactive({
    if (!identical(input$ess_plot %||% "track", "track")) return("520px")
    n <- length(intersect(input$ess_chromosome %||% character(0), unique(windows()$chromosome)))
    if (n == 0L) return(sprintf("%dpx", max(320L, 40L * length(unique(windows()$chromosome)) + 120L)))
    if (n == 1L) return("620px")
    sprintf("%dpx", min(1600L, 200L * n + 120L))
  })
  output$ess_plot_frame <- renderUI({ req(nrow(runs()) > 0); plotly::plotlyOutput("ess_plotly", height = ess_plot_height()) })
  output$ess_plotly <- plotly::renderPlotly({
    choice <- input$ess_plot %||% "track"
    if (identical(choice, "track")) {
      w <- windows(); req(nrow(w) > 0)
      chosen <- intersect(input$ess_chromosome %||% character(0), unique(w$chromosome))
      labels <- qc_plot_chromosome_display(chosen, active())
      view <- ess_view()
      # None chosen: the genome at a glance. One: the full track. Several: their landscapes side by side.
      if (!length(chosen)) {
        all_chroms <- unique(w$chromosome)
        return(ess_genome_overview_plotly(w, genes(), all_chroms, qc_plot_chromosome_display(all_chroms, active()), view))
      }
      if (length(chosen) > 1L) return(ess_multi_chromosome_plotly(w, genes(), chosen, labels, view))
      if (identical(view, "hybrid")) return(ess_gene_track_plotly(w, genes(), chosen, labels))
      return(ess_track_plotly(w, chosen, labels))
    }
    out <- ytab_ggplotly(ess_ggplot(choice, genes(), manifest = manifest(), run_dir = selected_run(), view = ess_view()))
    if (identical(choice, "calls")) {
      # ggplotly draws the facet y title as an annotation that runs into the tick labels: use the axis title instead
      out$x$layout$annotations <- Filter(function(a) !identical(trimws(gsub("<[^>]+>", "", as.character(a$text %||% ""))), "Genes"), out$x$layout$annotations)
      out <- plotly::layout(out, showlegend = FALSE, yaxis = list(automargin = TRUE, title = list(text = "<b>Genes</b>", standoff = 18)))
    }
    out
  })
  # One definition of the gene table, used by the screen and by the download: two definitions is how a
  # download ends up carrying columns the user was never shown.
  ess_genes_display <- reactive({
    g <- genes(); req(nrow(g) > 0)
    # The columns follow the view, like the cards and the plots: an Agreement column beside both models' scores
    # is the comparison again, made whether or not it was asked for.
    model_columns <- switch(ess_view(),
      hybrid  = c("hybrid_essential", "hybrid_score"),
      compare = c("hmm_call", "hmm_score", "hybrid_essential", "hybrid_score", "agreement"),
      c("hmm_call", "hmm_score", "hmm_status"))
    keep <- intersect(c("ORF", "Name", "chromosome", model_columns, "mappable_fraction", "reference_label"), names(g))
    g <- g[, keep, drop = FALSE]
    for (column in intersect(c("hmm_score", "hybrid_score", "mappable_fraction"), names(g))) g[[column]] <- round(g[[column]], 3)
    if ("chromosome" %in% names(g)) g$chromosome <- qc_plot_chromosome_display(g$chromosome, active())
    if ("agreement" %in% names(g)) g$agreement <- ifelse(g$agreement %in% names(essentiality_agreement_labels), essentiality_agreement_labels[g$agreement], g$agreement)
    if ("hybrid_essential" %in% names(g))
      g$hybrid_essential <- ifelse(as.logical(g$hybrid_essential), "Essential", "Not essential")
    labels <- c(ORF = "ORF", Name = "Name", chromosome = "Chromosome", hmm_call = "HMM call", hmm_score = "HMM score",
                hybrid_essential = "Hybrid call", hybrid_score = "Hybrid score",
                agreement = "Agreement", hmm_status = "HMM status", mappable_fraction = "Mappable fraction",
                reference_label = "Curated label")
    names(g) <- labels[names(g)]
    g
  })
  output$ess_genes_table <- DT::renderDT(ytab_datatable(ess_genes_display(), page_length = 15L))

  # The track is saved by capturing the figure on screen, so nothing rebuilds it here. What remains is for the
  # plots still written server-side.
  plot_for_download <- function() {
    choice <- input$ess_plot %||% "track"
    validate(need(!identical(choice, "track"), "The chromosome track is saved from the figure on screen."))
    ess_ggplot(choice, genes(), manifest = manifest(), run_dir = selected_run(), view = ess_view(),
               reference = reference_essentiality(), project_label = project_display_label(active()))
  }
  download_name <- function(ext) {
    chosen <- intersect(input$ess_chromosome %||% character(0), unique(windows()$chromosome))
    # The file is named by the chromosomes the project defines, not the reference's accessions.
    where <- if (!identical(input$ess_plot %||% "track", "track")) "" else if (!length(chosen)) "_genome"
             else paste0("_", paste(gsub("[^A-Za-z0-9]+", "", qc_plot_chromosome_display(chosen, active())), collapse = "-"))
    paste0(project_display_label(active()), "_essentiality_", ess_view(), "_", input$ess_plot %||% "track", where, "_", basename(selected_run()), ext)
  }
  for (format in c("png", "pdf")) local({
    ext <- format
    output[[paste0("ess_download_", ext)]] <- downloadHandler(filename = function() download_name(paste0(".", ext)),
      content = function(file) ytab_save_ggplot(file, plot_for_download(), width = 12, height = 6))
  })
  ytab_plot_capture_server(input, output, session, "ess", "ess_plotly",
                           function() sub("[.]$", "", download_name("")))
  output$ess_download_genes <- downloadHandler(
    filename = function() paste0(project_display_label(active()), "_essentiality_genes_", ess_view(), "_", basename(selected_run()), ".tsv"),
    content = function(file) utils::write.table(ess_genes_display(), file, sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
  invisible(list(runs = runs, selected_run = selected_run, manifest = manifest, genes = genes, windows = windows,
                 inventory = inventory, plot_for_download = plot_for_download, genes_display = ess_genes_display))
}
