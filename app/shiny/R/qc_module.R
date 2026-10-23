# Quality Control: the library QC table (`ytab qc`) and every QC plot, one at a time.
# Overview: the flagged table (hover a cell for its reason and normalization hint).
# Reads & Mapping, Insertions & Genes, Saturation & Bias: a plot menu (each table metric across libraries, coloured by its
# flag, plus the theme's curves and summaries), with the shared display options. Library Agreement: correlation between
# libraries within the project or with other projects' libraries.

qc_flag_colours <- c(red = "#F4C7C3", yellow = "#FCE8B2", green = "#B7E1CD", grey = "#EEEEEE")
qc_flag_line <- c(red = "#C5221F", yellow = "#B06000", green = "#188038", grey = "#9AA0A6")


# Each subtab's plot menu, grouped by question: label -> value ("metric:<key>" plots a QC-table metric across
# libraries). Metrics that only rescale another (insertions per kb and catalogue density of insertion sites; mean reads
# per catalogue site of reads at insertion sites; max reads of the top site's share) are not listed separately: they
# show in that plot's hover. Every metric stays in the QC table and its downloads.
qc_plot_menus <- list(
  reads_mapping = list(
    "Reads" = c("Reads sequenced" = "metric:input_reads", "Reads retained after cleanup" = "metric:retained_fraction",
                "Reads with 3' adaptor" = "metric:adapter_fraction"),
    "Mapping" = c("Total and mapped reads" = "read_counts", "Mapped (%)" = "metric:percent_mapped",
                  "Mapped with MAPQ above threshold (%)" = "metric:percent_mapq", "Mean MAPQ" = "metric:mean_mapq")),
  insertions = list(
    "Depth and complexity" = c("Reads at insertion sites" = "metric:total_cts", "Insertion sites" = "metric:unique_sites"),
    "Count distribution" = c("Count histogram" = "dist_histogram", "Counts against a geometric fit (QQ)" = "dist_qq", "Ranked counts" = "dist_ranked",
                             "Mean reads per occupied site" = "metric:NZmean", "Median reads per occupied site" = "metric:NZmedian",
                             "Skewness" = "metric:skewness", "Kurtosis" = "metric:kurtosis", "Tail index (Pickands)" = "metric:pickands_tail_index",
                             "Top site's share of reads" = "metric:max_share", "Top 1% of sites' share of reads" = "metric:top1pct_share",
                             "Strand correlation (W vs C)" = "metric:fr_corr"),
    "Genes" = c("Genes hit (%)" = "metric:percent_features_hit", "Genes hit, libraries of a role combined" = "combined_features",
                "Insertions per gene" = "metric:mean_hits_per_feature", "Genes by insertion count" = "gene_hit_bins",
                "Insertions in genes (%)" = "metric:percent_hits_in_features", "Insertion-neighbourhood correlation" = "metric:neighborhood_spearman",
                "Insertion signal along chromosomes" = "genome_bins")),
  saturation_bias = list(
    "Saturation" = c("Saturation curves" = "midlc", "MidLC estimate" = "metric:midlc_est", "Depth / MidLC" = "metric:depth_ratio"),
    "Jackpots" = c("Reads in jackpot sites" = "metric:jackpot_frac_reads", "Jackpot sites" = "metric:n_jackpot_sites"),
    "Positional bias" = c("Centromere bias" = "centromere", "Metaplot around TSS" = "metaplot_TSS", "Metaplot around TTS" = "metaplot_TTS",
                          "Metaplot around tRNA genes" = "metaplot_tRNA", "Sequence bias (+2/+7)" = "sequence_bias")),
  agreement = list(
    "Library agreement" = c("Correlation between libraries" = "cor_heatmap", "Gene insertion heatmap" = "gene_heatmap",
                            "Two libraries, gene by gene" = "pair_scatter")))

# Rescaled companions shown in a metric plot's hover.
qc_metric_companions <- list(unique_sites = c("insertions_per_kb", "density"), total_cts = c("mean_ct"), max_share = c("max_ct"),
                             jackpot_frac_reads = c("total_cts"), percent_mapped = c("mapped_reads"))

qc_menu_label <- function(id, value) { all <- unlist(unname(qc_plot_menus[[id]])); names(all)[match(value, all)] }

qc_plot_descriptions <- c(
  metric = "One QC metric across the libraries, coloured by its flag (green within the cutoff, yellow worth a look, red outside the cutoff, grey not flagged). Hover a bar for the value and the reason.",
  read_counts = "Alignment records and primary mapped reads per library.",
  dist_histogram = "Histogram of read counts at occupied positions (log bins). Well-behaved libraries have many low counts and few high ones.",
  dist_qq = "Observed count quantiles against a geometric distribution fitted to each library's mean. Strong curvature away from the diagonal means a skewed library; BGC suits those.",
  dist_ranked = "Read counts at positions in decreasing order. A few isolated outliers are expected (TTR is robust to them); a steep head means a few sites dominate.",
  gene_hit_bins = "Number of genes with 0, 1–3, 4–7, 8–9, 10–19 and 20 or more insertion sites, per library.",
  combined_features = "Genes with at least one insertion when all libraries of a role are combined.",
  genome_bins = "Mean insertion signal of genome bins, per chromosome, across libraries; dips mark regions with few insertions.",
  midlc = "Unique insertion sites found as reads are subsampled. A curve still rising at full depth means the library is not saturated.",
  centromere = "Mean reads by distance from the centromere, across chromosome arms.",
  metaplot_TSS = "Insertion signal around transcription start sites.", metaplot_TTS = "Insertion signal around transcription end sites.",
  metaplot_tRNA = "Insertion signal around tRNA genes (Hermes favours tRNA-flanking regions).",
  sequence_bias = "Enrichment of each base pair at positions +2 and +7 of the insertion site over the genome background; the dashed line at 1 is no enrichment. Bars are the mean of each role's libraries, or each library.",
  cor_heatmap = "Spearman correlation of gene insertion density between libraries, ordered by similarity. Hover for the pair type and the 95% interval.",
  gene_heatmap = "Gene insertion density (log scale, centred per gene) for the genes that vary most across the libraries.",
  pair_scatter = "Gene insertion density of two libraries, gene by gene; the Spearman correlation is in the title.")

qc_display_controls <- function(id, heatmap = FALSE) ytab_plot_customization_controls(id, include_points = TRUE, include_bars = TRUE,
  include_heatmap = heatmap, default_height = "medium", default_bar_orientation = "auto", default_show_value_labels = FALSE)

qc_theme_ui <- function(id, extra = NULL) ytab_two_column_layout(
  controls = tagList(ytab_control_panel("Plot",
      uiOutput(paste0("qcp_", id, "_choice")),
      uiOutput(paste0("qcp_", id, "_options")),
      ytab_plot_downloads_ui(paste0("qcp_", id), c(data = "Download plotted data (TSV)")),
      qc_display_controls(paste0("qcp_", id), heatmap = id %in% c("saturation_bias", "agreement"))),
    extra),
  main = tagList(uiOutput(paste0("qcp_", id, "_card"))))

qc_overview_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Library QC table",
    actionButton("qc_table_run", "Build QC table", class = "btn-primary"),
    helpText("Collects every per-library QC metric (reads, mapping, insertion statistics, genes, saturation) from the project's outputs. Takes seconds."),
    tags$hr(),
    uiOutput("qc_table_result_selector"),
    tags$div(class = "ytab-actions",
      downloadButton("qc_table_download_xlsx", "Download table (Excel, coloured)"),
      downloadButton("qc_table_download_tsv", "Download table (TSV)"),
      downloadButton("qc_table_download_flags", "Download flags and comments (TSV)"))),
  main = tagList(
    uiOutput("qc_table_summary"),
    ytab_table_card("QC table", DT::DTOutput("qc_table"),
      description = qc_flag_legend("Hover a cell for its reason, the metric's meaning and, where one applies, a normalization hint. Every metric is plotted across libraries in the subtabs.")),
    ytab_methods_drawer(
      tags$p(tags$b("Rule of thumb. "), "Skewness of the read counts at occupied positions: above 30 often signals a problem (yellow) and below 50 is desired (above 50 red). It can be higher under strong selection."),
      tags$p(tags$b("Compared with the project. "), "Every other flagged metric compares a library with the project's libraries of the same role (parents with parents, treated with treated): robust z = (value − median) / (1.4826 × MAD), in the metric's bad direction only; yellow beyond 2.5 and red beyond 3 (Leys et al. 2013). Needs at least 4 libraries of the role, otherwise the cell is grey. The spread is at least 10% of the median, so near-identical libraries are not flagged for trivial differences."),
      tags$p(tags$b("Hermes. "), "Hermes has no fixed set of possible insertion sites: density is computed over the project's catalogue (positions with an insertion in any library) and insertions per kb measure coverage, so there are no absolute density or depth targets. The strand correlation is shown but not flagged: the two orientations at a position are separate insertions."),
      tags$p(tags$b("Normalization hints. "), "Skewed or jackpot-dominated libraries suit BGC (betageom); TTR is robust to isolated outlier sites."),
      ytab_reference("The count statistics, the count-distribution plots and the skewness rule of thumb are adapted from TRANSIT (DeJesus et al. 2015, PLoS Comput Biol 11:e1004401), whose guidance was written for Himar1."),
      tags$p("Command: ", tags$code("ytab qc"), ". Results are written to the project's ", tags$code("qc"), " folder."))))

qc_library_diagnostics_run_ui <- function() tags$details(class = "ytab-more-options", tags$summary("Run library diagnostics (saturation, jackpots, positional bias)"),
  uiOutput("library_diagnostics_readiness"), uiOutput("qc_current_selection"),
  selectizeInput("qc_samples", "Samples to diagnose", choices = character(), multiple = TRUE, options = list(placeholder = "Select samples…")),
  uiOutput("qc_sample_presets"), uiOutput("qc_selection_message"),
  radioButtons("library_diagnostics_mode", "Execution mode", c("Preview command" = "preview", "Run diagnostics" = "run"), "preview", inline = TRUE),
  uiOutput("qc_execution_summary"), uiOutput("library_diagnostics_cache_ui"), uiOutput("library_diagnostics_action"),
  job_progress_ui("library_diagnostics_job"), uiOutput("library_diagnostics_result"),
  tags$details(class = "ytab-technical-details", tags$summary("Technical details"), verbatimTextOutput("library_diagnostics_technical")))

qc_agreement_run_ui <- function() ytab_control_panel("Compare libraries",
  radioButtons("libcor_scope", NULL, c("Within this project" = "within", "With other projects" = "across"), "within"),
  conditionalPanel("input.libcor_scope === 'across'", uiOutput("libcor_project_selector"), uiOutput("libcor_ortholog_selector")),
  radioButtons("libcor_metric", "Gene summary", c("Insertion density (sites per kb)" = "density", "Reads per kb" = "reads_per_kb"), "density"),
  helpText("Density compares libraries from different pools. Reads per kb also reflects clone abundance, so use it only within a pool."),
  numericInput("libcor_min_reads", "Reads for a site to count as an insertion", 1, min = 1, step = 1),
  actionButton("libcor_run", "Run library correlation", class = "btn-primary"),
  tags$hr(), uiOutput("libcor_result_selector"),
  tags$div(class = "ytab-actions", downloadButton("libcor_download_matrix", "Download matrix (TSV)"), downloadButton("libcor_download_pairs", "Download pairs (TSV)")))

quality_control_ui <- function() navset_tab(id = "qc_tabs",
  nav_panel("Overview", value = "overview", qc_overview_ui()),
  nav_panel("Reads & Mapping", value = "reads_mapping", qc_theme_ui("reads_mapping")),
  nav_panel("Insertions & Genes", value = "insertions", qc_theme_ui("insertions")),
  nav_panel("Saturation & Bias", value = "saturation_bias", qc_theme_ui("saturation_bias", qc_library_diagnostics_run_ui())),
  nav_panel("Library Agreement", value = "agreement", qc_theme_ui("agreement", qc_agreement_run_ui()),
    uiOutput("libcor_summary"), ytab_table_card("Library pairs", DT::DTOutput("libcor_pairs_table")),
    ytab_methods_drawer(
      tags$p("Each library is summarized per gene as insertion sites per kb (or reads per kb). Two libraries are compared by the Spearman correlation of these values across the genes they share; the 95% interval is from a bootstrap over genes. Libraries of different species are compared over 1:1 orthologs."),
      tags$ul(
        tags$li(tags$b("Same pool (shared clones): "), "the libraries share insertion clones, so a high correlation reflects library structure, not biological similarity."),
        tags$li(tags$b("Same strain, different pool: "), "independent libraries; the correlation measures how reproducible the insertion landscape is."),
        tags$li(tags$b("Different strain: "), "strain differences in essentiality or accessibility, plus library noise."),
        tags$li(tags$b("Different species: "), "conservation of the insertion landscape over 1:1 orthologs.")),
      tags$p("Clone evidence uses the read correlation at shared insertion sites: above 0.7 indicates shared clones and below 0.5 independent libraries (calibrated on Gale et al. 2020: same pool 0.89–0.99, independent pools 0.42–0.43). A pair is flagged when this evidence disagrees with the pool labels in the sample sheet."),
      ytab_reference("Generalized from the correlation plot of TRANSIT (DeJesus et al. 2015, PLoS Comput Biol 11:e1004401)."),
      tags$p("Command: ", tags$code("ytab library-correlation"), ". Results are written to the project's ", tags$code("library_correlation"), " folder."))))

# The flag colours in words, each word on its own colour.
qc_flag_legend <- function(...) {
  chip <- function(flag, text) tags$span(style = sprintf("background:%s;padding:1px 6px;border-radius:3px;font-weight:600", qc_flag_colours[[flag]]), text)
  tags$span(chip("green", "Green: within the cutoff"), " ", chip("yellow", "Yellow: worth a look"), " ", chip("red", "Red: outside the cutoff"), " ",
            chip("grey", "Grey: not flagged"), ". ", ...)
}

qc_format_value <- function(x) {
  out <- ifelse(is.na(x), "n/a", ifelse(abs(x) >= 1e5, format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE), formatC(signif(x, 4), format = "fg", flag = "#")))
  sub("\\.$", "", out)
}

qc_html_escape <- function(x) gsub("\"", "&quot;", htmltools::htmlEscape(x), fixed = TRUE)

# Library x metric display table: each cell a coloured block whose tooltip carries the flag's reason and comment.
qc_display_table <- function(table, flags, metrics) {
  out <- data.frame(Library = table$display_name, Role = table$role, check.names = FALSE, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(metrics))) {
    key <- metrics$metric[[i]]; f <- flags[flags$metric == key, , drop = FALSE]; f <- f[match(table$sample, f$sample), , drop = FALSE]
    tip <- paste0(f$label, ": ", f$reason, "\n", f$comment, ifelse(nzchar(f$normalization_hint %||% ""), paste0("\n→ ", f$normalization_hint), ""))
    out[[metrics$label[[i]]]] <- sprintf('<span class="ytab-qc-cell" style="display:block;padding:2px 6px;border-radius:3px;background:%s" title="%s">%s</span>',
      qc_flag_colours[ifelse(is.na(f$flag), "grey", f$flag)], qc_html_escape(tip), qc_format_value(suppressWarnings(as.numeric(f$value))))
  }
  out
}

# Library names shortened by dropping the "-"-separated words every name shares (label mode "compact").
qc_compact_labels <- function(labels) {
  parts <- strsplit(as.character(labels), "[-_ ]+"); common <- Reduce(intersect, parts)
  if (length(labels) < 2L || !length(common)) return(as.character(labels))
  short <- vapply(parts, function(p) paste(setdiff(p, common), collapse = "-"), "")
  ifelse(nzchar(short) & !duplicated(short) & !duplicated(short, fromLast = TRUE), short, as.character(labels))
}

# The shared display options applied to a ggplot: style, text size, label angle and mode, grid.
qc_apply_display <- function(p, s, category_axis = "x") {
  base <- switch(s$style %||% "clean", presentation = ggplot2::theme_bw(base_size = s$text_size * 1.05), ggplot2::theme_minimal(base_size = s$text_size * if (identical(s$style, "compact")) .85 else 1))
  p <- p + base + ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))
  if (identical(s$grid, "hide")) p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
  angle <- if (is.na(s$label_angle)) 30L else s$label_angle
  text <- if (identical(s$label_mode, "hide")) ggplot2::element_blank() else ggplot2::element_text(angle = angle, hjust = if (angle > 0) 1 else .5, vjust = if (angle == 90) .5 else 1)
  if (identical(category_axis, "x")) p + ggplot2::theme(axis.text.x = text) else p + ggplot2::theme(axis.text.y = if (identical(s$label_mode, "hide")) ggplot2::element_blank() else ggplot2::element_text())
}

qc_bar_horizontal <- function(s, labels) identical(s$bar_orientation, "horizontal") || (identical(s$bar_orientation, "auto") && max(nchar(as.character(labels)), 0L) > 12L)

qc_label_libraries <- function(x, s) if (identical(s$label_mode, "compact")) qc_compact_labels(x) else as.character(x)

qc_metric_plot <- function(f, s, companions = NULL) {
  f$value <- suppressWarnings(as.numeric(f$value)); f <- f[!is.na(f$value), , drop = FALSE]; validate(need(nrow(f) > 0, "No values for this metric."))
  f$library <- qc_label_libraries(f$display_name, s); horizontal <- qc_bar_horizontal(s, f$library)
  f$library <- factor(f$library, levels = if (horizontal) rev(unique(f$library)) else unique(f$library))
  extra <- if (is.data.frame(companions) && nrow(companions)) vapply(f$sample, function(x) {
    c2 <- companions[companions$sample == x, , drop = FALSE]; if (!nrow(c2)) return("")
    paste0("\n", paste(sprintf("%s: %s", c2$label, qc_format_value(suppressWarnings(as.numeric(c2$value)))), collapse = "\n"))
  }, "") else ""
  f$tip <- sprintf("%s\n%s: %s%s\n%s", f$display_name, f$label, qc_format_value(f$value), extra, f$reason)
  p <- ggplot2::ggplot(f, ggplot2::aes(x = library, y = value, fill = flag, text = tip)) + ggplot2::geom_col(width = .7) +
    ggplot2::scale_fill_manual(values = qc_flag_line, drop = FALSE) + ggplot2::labs(x = NULL, y = f$label[[1L]], fill = "Flag")
  if (isTRUE(s$show_value_labels)) p <- p + ggplot2::geom_text(ggplot2::aes(label = qc_format_value(value)), size = s$text_size / 4.2,
                                                               hjust = if (horizontal) -.1 else .5, vjust = if (horizontal) .5 else -.4)
  if (horizontal) qc_apply_display(p + ggplot2::coord_flip(), s, "y") else qc_apply_display(p, s, "x")
}

qc_bar_plot <- function(d, x, y, fill, s, ylab) {
  d$.x <- qc_label_libraries(d[[x]], s); horizontal <- qc_bar_horizontal(s, d$.x)
  d$.x <- factor(d$.x, levels = if (horizontal) rev(unique(d$.x)) else unique(d$.x)); d$.tip <- sprintf("%s\n%s: %s", d[[x]], d[[fill]], qc_format_value(d[[y]]))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .x, y = .data[[y]], fill = .data[[fill]], text = .tip)) + ggplot2::geom_col(position = "dodge", width = .75) +
    ggplot2::labs(x = NULL, y = ylab, fill = NULL)
  if (isTRUE(s$show_value_labels)) p <- p + ggplot2::geom_text(ggplot2::aes(label = qc_format_value(.data[[y]])), position = ggplot2::position_dodge(.75), size = s$text_size / 4.2, vjust = -.4)
  if (horizontal) qc_apply_display(p + ggplot2::coord_flip(), s, "y") else qc_apply_display(p, s, "x")
}

qc_line_plot <- function(d, x, y, colour, s, xlab, ylab, logx = FALSE, logy = FALSE, diagonal = FALSE) {
  d$.tip <- sprintf("%s\n%s: %s\n%s: %s", d$.library, xlab, qc_format_value(d[[x]]), ylab, qc_format_value(d[[y]]))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour]], group = .library, text = .tip))
  if (diagonal) p <- p + ggplot2::geom_abline(linetype = 2, colour = "grey50")
  p <- p + ggplot2::geom_line(linewidth = .6, alpha = max(.3, s$point_opacity %||% .8)) + ggplot2::labs(x = xlab, y = ylab, colour = NULL)
  if (logx) p <- p + ggplot2::scale_x_log10()
  if (logy) p <- p + ggplot2::scale_y_log10()
  qc_apply_display(p, utils::modifyList(s, list(label_angle = 0L)), "x")
}

qc_heatmap_fill <- function(palette, midpoint = 0, limits = NULL) switch(palette %||% "blue_red",
  viridis = ggplot2::scale_fill_viridis_c(option = "D", limits = limits, na.value = "grey90"),
  magma = ggplot2::scale_fill_viridis_c(option = "A", limits = limits, na.value = "grey90"),
  gray = ggplot2::scale_fill_gradient(low = "white", high = "black", limits = limits, na.value = "grey90"),
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = midpoint, limits = limits, na.value = "grey90"))

qc_heatmap_theme <- function(p, s) {
  p + ggplot2::theme_minimal(base_size = s$text_size) + ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right",
    axis.text.x = ggplot2::element_text(angle = s$column_label_angle %||% 45, hjust = 1, size = s$row_label_size),
    axis.text.y = if (isTRUE(s$show_row_labels)) ggplot2::element_text(size = s$row_label_size) else ggplot2::element_blank())
}

# Library identity for the plots: display name, role and pool from the sample sheet.
qc_library_info <- function(project) {
  s <- project$samples
  if (!is.data.frame(s) || !nrow(s)) return(data.frame(sample = character(), display = character(), role = character(), pool = character()))
  pick <- function(cols) { for (col in cols) if (col %in% names(s)) { v <- as.character(s[[col]]); if (any(nzchar(v) & !is.na(v))) return(ifelse(is.na(v), "", v)) }; rep("", nrow(s)) }
  include <- if ("include" %in% names(s)) tolower(as.character(s$include)) %in% c("true", "1", "yes") else rep(TRUE, nrow(s))
  out <- data.frame(sample = as.character(s$sample), display = ytab_sample_labels(s$sample, project), role = tolower(pick(c("library_role", "condition", "guessed_condition"))),
                    pool = pick(c("pool", "guessed_pool")), stringsAsFactors = FALSE)[include, , drop = FALSE]
  out$role[!nzchar(out$role)] <- "library"
  out
}

qc_annotate <- function(d, info, sample_col = "sample") {
  i <- match(d[[sample_col]], info$sample)
  d$.library <- ifelse(is.na(i), as.character(d[[sample_col]]), info$display[i]); d$.role <- ifelse(is.na(i), "library", info$role[i]); d$.pool <- ifelse(is.na(i), "", info$pool[i])
  d
}

qc_colour_column <- function(colour_by) switch(colour_by %||% "library", role = ".role", pool = ".pool", ".library")

qc_table_server <- function(input, output, session, repo_root, active, active_project_path, jobs, python_bin) {
  runs_dir <- reactive(file.path(active()$project_root, "qc"))
  runs <- reactivePoll(3000, session,
    checkFunc = function() { x <- ytab_analysis_runs(tryCatch(runs_dir(), error = function(e) "")); paste(x$path, x$modified, collapse = "|") },
    valueFunc = function() ytab_analysis_runs(tryCatch(runs_dir(), error = function(e) "")))
  output$qc_table_result_selector <- renderUI({
    x <- runs()
    if (!nrow(x)) return(tags$p(class = "text-muted", "No QC table for this project yet."))
    selectInput("qc_table_result", "QC table", ytab_run_choices(x))
  })
  observeEvent(input$qc_table_run, {
    if (jobs$job_is_running()) { showNotification("Wait for the current job to finish.", type = "warning"); return() }
    args <- c(file.path(repo_root, "scripts", "local", "ytab_qc.py"), "--project-config", active_project_path(), "--output", ytab_new_run_dir(runs_dir(), "qc"))
    tryCatch({
      jobs$start_job(python_bin(), args, active_project_path(), "library_qc", wd = repo_root)
      showNotification("Building the QC table. It appears here when it finishes.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error"))
  }, ignoreInit = TRUE)

  has_table <- reactive(nrow(runs()) > 0)
  selected_run <- reactive({ req(has_table()); path <- input$qc_table_result %||% runs()$path[[1L]]; req(path %in% runs()$path); path })
  manifest <- reactive(ytab_read_manifest(selected_run()))
  table <- reactive(ytab_read_analysis_tsv(file.path(selected_run(), "qc_table.tsv")))
  metrics <- reactive(ytab_read_analysis_tsv(file.path(selected_run(), "metrics.tsv")))
  flags <- reactive({
    f <- ytab_read_analysis_tsv(file.path(selected_run(), "qc_flags.tsv")); t <- table(); req(nrow(f) > 0, nrow(t) > 0)
    f$display_name <- t$display_name[match(f$sample, t$sample)]; f
  })
  info <- reactive(qc_library_info(active()))
  distribution <- function(kind) {
    validate(need(has_table(), "Build the QC table in Overview first."))
    d <- ytab_read_analysis_tsv(file.path(selected_run(), "plots", paste0(kind, ".tsv.gz")))
    validate(need(nrow(d) > 0, "No count distributions in this QC table."))
    qc_annotate(d, info())
  }

  output$qc_table_summary <- renderUI({
    if (!has_table()) return(ytab_empty_state("No QC table yet", "Build the QC table to collect every library's QC metrics in one place, flagged and commented."))
    m <- manifest(); req(m)
    if (length(m$libraries_without_hits)) ytab_warning_banner("Libraries without a hit file", paste(unlist(m$libraries_without_hits), collapse = ", "))
  })
  output$qc_table <- DT::renderDT({
    t <- table(); f <- flags(); m <- metrics(); req(nrow(t) > 0, nrow(m) > 0)
    DT::datatable(qc_display_table(t, f, m), escape = FALSE, rownames = FALSE, selection = "none",
                  class = "compact cell-border ytab-qc-grid", options = list(scrollX = FALSE, autoWidth = FALSE, pageLength = 25, dom = "tip", ordering = FALSE))
  })
  download_name <- function(suffix) paste0(active()$project_id, "_", basename(selected_run()), suffix)
  output$qc_table_download_xlsx <- downloadHandler(filename = function() download_name(".xlsx"), content = function(file) file.copy(file.path(selected_run(), "qc_table.xlsx"), file))
  output$qc_table_download_tsv <- downloadHandler(filename = function() download_name("_table.tsv"), content = function(file) file.copy(file.path(selected_run(), "qc_table.tsv"), file))
  output$qc_table_download_flags <- downloadHandler(filename = function() download_name("_flags.tsv"), content = function(file) file.copy(file.path(selected_run(), "qc_flags.tsv"), file))

  # --- library correlation results for Library Agreement --------------------------------------------------------------
  libcor_dir <- reactive(file.path(active()$project_root, "library_correlation"))
  libcor_runs <- reactivePoll(3000, session,
    checkFunc = function() { x <- ytab_analysis_runs(tryCatch(libcor_dir(), error = function(e) "")); paste(x$path, x$modified, collapse = "|") },
    valueFunc = function() ytab_analysis_runs(tryCatch(libcor_dir(), error = function(e) "")))
  libcor_results <- reactive({
    x <- libcor_runs(); validate(need(nrow(x) > 0, "No library correlation yet. Choose what to compare and run it."))
    path <- input$libcor_result %||% x$path[[1L]]; req(path %in% x$path)
    list(path = path, matrix = libcor_matrix(path), pairs = ytab_read_analysis_tsv(file.path(path, "pairs.tsv")),
         samples = ytab_read_analysis_tsv(file.path(path, "samples.tsv")), genes = ytab_read_analysis_tsv(file.path(path, "gene_values.tsv.gz")))
  })
  libcor_samples <- reactive(tryCatch(libcor_results()$samples, error = function(e) data.frame()))

  # --- themed subtabs: one plot at a time --------------------------------------------------------------------------
  plot_choices <- function(id) lapply(qc_plot_menus[[id]], as.list)
  sources <- list(
    flags = function() { validate(need(has_table(), "Build the QC table in Overview first.")); flags() },
    distribution = distribution, project = active, info = info, libcor = libcor_results)
  built <- list()
  for (theme_id in names(qc_plot_menus)) local({
    id <- theme_id; prefix <- paste0("qcp_", id)
    output[[paste0(prefix, "_choice")]] <- renderUI({
      tagList(selectInput(paste0(prefix, "_plot"), "Show", plot_choices(id), selected = isolate(input[[paste0(prefix, "_plot")]])),
              if (!identical(id, "agreement") && !has_table()) helpText("Plots of QC-table metrics need the QC table (build it in Overview)."))
    })
    output[[paste0(prefix, "_options")]] <- renderUI(qc_plot_options(id, input[[paste0(prefix, "_plot")]] %||% "", info(), libcor_samples(), dirname(active()$project_root)))
    plot_reactive <- reactive({
      choice <- input[[paste0(prefix, "_plot")]]; req(choice)
      qc_build_plot(choice, input, prefix, ytab_plot_customization_values(input, prefix), sources)
    })
    built[[id]] <<- plot_reactive
    output[[paste0(prefix, "_card")]] <- renderUI({
      choice <- input[[paste0(prefix, "_plot")]] %||% ""; s <- ytab_plot_customization_values(input, prefix)
      key <- if (startsWith(choice, "metric:")) "metric" else choice
      label <- qc_menu_label(id, choice)
      height <- if (choice %in% c("cor_heatmap", "gene_heatmap")) s$heatmap_height else s$height
      ytab_plot_card(if (is.na(label %||% NA)) "Plot" else label, ytab_plot_frame(plotly::plotlyOutput(paste0(prefix, "_plotly"), height = height), s$width),
                     description = if (is.na(qc_plot_descriptions[key])) "" else unname(qc_plot_descriptions[key]))
    })
    output[[paste0(prefix, "_plotly")]] <- plotly::renderPlotly({
      b <- plot_reactive(); ytab_ggplotly(b$plot) |> plotly::layout(legend = list(orientation = "h", y = -0.32, yanchor = "top"), margin = list(b = 90))
    })
    for (format in c("png", "pdf")) local({
      ext <- format
      output[[paste0(prefix, "_download_", ext)]] <- downloadHandler(
        filename = function() paste0(active()$project_id, "_qc_", gsub("[^A-Za-z0-9]+", "_", input[[paste0(prefix, "_plot")]] %||% id), ".", ext),
        content = function(file) { s <- ytab_plot_customization_values(input, prefix)
          ytab_save_ggplot(file, plot_reactive()$plot, width = switch(s$width, wide = 13, full = 16, 10), height = max(4, as.numeric(sub("px", "", s$height)) / 75)) })
    })
    output[[paste0(prefix, "_download_data")]] <- downloadHandler(
      filename = function() paste0(active()$project_id, "_qc_", gsub("[^A-Za-z0-9]+", "_", input[[paste0(prefix, "_plot")]] %||% id), ".tsv"),
      content = function(file) { d <- plot_reactive()$data; d <- d[, !startsWith(names(d), "."), drop = FALSE]
        utils::write.table(d, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "") })
  })

  invisible(list(runs = runs, selected_run = selected_run, table = table, flags = flags, metrics = metrics,
                 info = info, libcor_results = libcor_results, built = built))
}

# Plot-specific options next to the plot menu.
qc_plot_options <- function(id, choice, info, libcor_samples, projects_root = "") {
  roles <- sort(unique(info$role)); prefix <- paste0("qcp_", id)
  colour <- selectInput(paste0(prefix, "_colour"), "Colour by", c("Library" = "library", "Role" = "role", "Pool" = "pool"), "library")
  roles_filter <- if (length(roles) > 1L) checkboxGroupInput(paste0(prefix, "_roles"), "Libraries shown", roles, roles, inline = TRUE)
  if (startsWith(choice, "dist_") || choice %in% c("midlc", "centromere", "metaplot_TSS", "metaplot_TTS", "metaplot_tRNA")) return(tagList(roles_filter, colour))
  if (identical(choice, "sequence_bias"))
    return(tagList(roles_filter, selectInput(paste0(prefix, "_colour"), "Bars", c("Mean of each role" = "role", "Each library" = "library", "Mean of each pool" = "pool"), "role")))
  if (choice %in% c("gene_hit_bins", "read_counts")) return(roles_filter)
  libcor_role_filter <- function() {
    if (!nrow(libcor_samples)) return(NULL)
    present <- sort(unique(libcor_roles(libcor_samples$name, projects_root, libcor_samples)))
    if (length(present) > 1L) checkboxGroupInput(paste0(prefix, "_roles"), "Libraries shown", present, present, inline = TRUE)
  }
  if (identical(choice, "cor_heatmap")) return(libcor_role_filter())
  if (identical(choice, "gene_heatmap")) return(tagList(libcor_role_filter(),
    numericInput(paste0(prefix, "_genes"), "Genes shown (most variable)", 300, min = 20, max = 3000, step = 50)))
  if (identical(choice, "pair_scatter") && nrow(libcor_samples)) {
    names <- stats::setNames(libcor_samples$name, libcor_display_names(libcor_samples$name, projects_root))
    return(tagList(selectInput(paste0(prefix, "_lib_a"), "Library A", names, names[[1L]]),
                   selectInput(paste0(prefix, "_lib_b"), "Library B", names, names[[min(2L, length(names))]])))
  }
  NULL
}

# Library-correlation libraries whose role is ticked in the "Libraries shown" filter (all when there is no filter).
qc_libcor_shown <- function(names, r, input, prefix, projects_root) {
  roles <- input[[paste0(prefix, "_roles")]]
  if (!length(roles)) return(names)
  names[libcor_roles(names, projects_root, r$samples) %in% roles]
}

qc_role_filter <- function(d, input, prefix) { roles <- input[[paste0(prefix, "_roles")]]; if (length(roles)) d[d$.role %in% roles, , drop = FALSE] else d }

qc_distribution_plot <- function(d, kind, s, colour) switch(kind,
  histogram = qc_line_plot(d[d$sites > 0, , drop = FALSE], "bin_start", "sites", colour, s, "Reads at a position (log bins)", "Positions", TRUE, TRUE),
  qq = qc_line_plot(d, "geometric", "observed", colour, s, "Geometric quantile (fit to the mean)", "Observed quantile", TRUE, TRUE, diagonal = TRUE),
  ranked = qc_line_plot(d, "rank", "reads", colour, s, "Rank of position", "Reads at the position", TRUE, TRUE))

# Build one QC plot: list(plot = ggplot, data = plotted data).
qc_build_plot <- function(choice, input, prefix, s, src) {
  colour <- qc_colour_column(input[[paste0(prefix, "_colour")]]); project <- src$project(); info <- src$info()
  need_rows <- function(d, what) { validate(need(is.data.frame(d) && nrow(d) > 0, paste0("No ", what, " for this project yet."))); d }
  if (startsWith(choice, "metric:")) {
    key <- sub("^metric:", "", choice); all <- src$flags(); f <- all[all$metric == key, , drop = FALSE]
    companions <- all[all$metric %in% qc_metric_companions[[key]], , drop = FALSE]
    return(list(plot = qc_metric_plot(f, s, companions), data = rbind(f, companions)))
  }
  if (startsWith(choice, "dist_")) {
    kind <- sub("^dist_", "", choice); d <- qc_role_filter(src$distribution(kind), input, prefix)
    return(list(plot = qc_distribution_plot(d, kind, s, colour), data = d))
  }
  switch(choice,
    read_counts = {
      d <- qc_role_filter(qc_annotate(need_rows(qc_mapping_stats_plot_data(project), "mapping statistics"), info), input, prefix)
      long <- rbind(data.frame(.library = d$.library, Reads = "Alignment records", value = as.numeric(d$total_records)),
                    data.frame(.library = d$.library, Reads = "Primary mapped", value = as.numeric(d$primary_mapped)))
      list(plot = qc_bar_plot(long, ".library", "value", "Reads", s, "Reads") + ggplot2::scale_fill_manual(values = c("#9AA0A6", "#1A73E8")), data = long)
    },
    gene_hit_bins = {
      d <- qc_role_filter(qc_annotate(need_rows(summary_qc_gene_hit_bins_raw(project), "SummaryTable gene counts"), info), input, prefix)
      d$bin <- factor(d[["Hit bin"]], levels = unique(d[["Hit bin"]])); d$.tip <- sprintf("%s\n%s insertion sites: %s genes", d$.library, d$bin, d$Genes)
      p <- ggplot2::ggplot(d, ggplot2::aes(x = bin, y = Genes, fill = .library, text = .tip)) + ggplot2::geom_col(position = "dodge") +
        ggplot2::labs(x = "Insertion sites in the gene", y = "Genes", fill = NULL)
      list(plot = qc_apply_display(p, utils::modifyList(s, list(label_angle = 0L)), "x"), data = d)
    },
    combined_features = {
      d <- need_rows(qc_summary_combined_features_hit_data(project), "SummaryTable results"); d$measure <- "Genes hit (%)"
      list(plot = qc_bar_plot(d, "display_group", "percent_features_hit", "measure", s, "Genes hit (%)") + ggplot2::scale_fill_manual(values = "#1A73E8", guide = "none"), data = d)
    },
    genome_bins = {
      d <- need_rows(qc_summary_binned_plot_data(project), "binned SummaryTable signal"); d$chrom <- factor(d$chrom, levels = unique(d$chrom))
      d$.tip <- sprintf("%s: %s", d$chrom, qc_format_value(d$mean_bin_signal))
      p <- ggplot2::ggplot(d, ggplot2::aes(x = chrom, y = mean_bin_signal)) + ggplot2::geom_boxplot(fill = "#D2E3FC", outlier.shape = NA) +
        ggplot2::geom_point(ggplot2::aes(text = .tip), position = ggplot2::position_jitter(width = .15, seed = 1), size = s$point_size, alpha = s$point_opacity) +
        ggplot2::labs(x = NULL, y = "Mean signal per genome bin")
      list(plot = qc_apply_display(p, s, "x"), data = d)
    },
    midlc = {
      d <- qc_role_filter(qc_annotate(need_rows(qc_library_midlc_plot_data(project), "saturation curves (run library diagnostics)"), info), input, prefix)
      list(plot = qc_line_plot(d, "reads_sampled", "unique_sites", colour, s, "Reads sampled", "Unique insertion sites", logx = TRUE), data = d)
    },
    centromere = {
      d <- qc_role_filter(qc_annotate(need_rows(qc_library_centromere_plot_data(project), "centromere profiles (run library diagnostics)"), info), input, prefix)
      list(plot = qc_line_plot(d, "distance_kb", "mean_reads", colour, s, "Distance from the centromere (kb)", "Mean reads"), data = d)
    },
    metaplot_TSS = , metaplot_TTS = , metaplot_tRNA = {
      feature <- sub("^metaplot_", "", choice); d <- need_rows(qc_library_metaplot_plot_data(project), "metaplots (run library diagnostics)")
      d <- need_rows(d[tolower(d$feature) == tolower(feature), , drop = FALSE], paste(feature, "metaplots"))
      d <- qc_role_filter(qc_annotate(d, info), input, prefix)
      list(plot = qc_line_plot(d, "rel_bp", "mean_reads_per_site", colour, s, sprintf("Position relative to the %s (bp)", feature), "Mean reads per site"), data = d)
    },
    sequence_bias = {
      # Grouped bars as in the former Library Diagnostics plot: one group per dinucleotide at +2/+7, one bar per role
      # (mean over the role's libraries) or per library; dashed line at 1 = no enrichment over the genome background.
      d <- qc_role_filter(qc_annotate(need_rows(qc_library_sequence_bias_plot_data(project), "sequence bias (run library diagnostics)"), info), input, prefix)
      d$enrichment <- as.numeric(d$enrichment); by <- input[[paste0(prefix, "_colour")]] %||% "role"
      d$.group <- if (identical(by, "library")) qc_label_libraries(d$.library, s) else if (identical(by, "pool")) d$.pool else d$.role
      bars <- stats::aggregate(enrichment ~ pair + .group, d, mean)
      bars$n <- stats::aggregate(enrichment ~ pair + .group, d, length)$enrichment
      bars$.tip <- sprintf("%s\n%s: %.2f%s", bars$.group, bars$pair, bars$enrichment, ifelse(bars$n > 1, sprintf(" (mean of %d libraries)", bars$n), ""))
      p <- ggplot2::ggplot(bars, ggplot2::aes(x = pair, y = enrichment, fill = .group, text = .tip)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = .8), width = .75) +
        ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "black") +
        ggplot2::labs(x = "Dinucleotide (+2/+7)", y = "Mean enrichment", fill = NULL)
      if (isTRUE(s$show_value_labels)) p <- p + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", enrichment)), position = ggplot2::position_dodge(width = .8),
                                                                   vjust = -.4, size = s$text_size / 4.2)
      list(plot = qc_apply_display(p, utils::modifyList(s, list(label_angle = 0L)), "x"), data = bars)
    },
    cor_heatmap = {
      r <- src$libcor(); values <- r$matrix; validate(need(!is.null(values), "This correlation result has no matrix."))
      keep <- qc_libcor_shown(rownames(values), r, input, prefix, dirname(project$project_root))
      validate(need(length(keep) >= 2L, "Show at least two libraries.")); values <- values[keep, keep, drop = FALSE]
      names <- rownames(values); long <- expand.grid(a = names, b = names, stringsAsFactors = FALSE); long$rho <- values[cbind(long$a, long$b)]
      pr <- r$pairs; k <- match(paste(long$a, long$b), paste(pr$sample_a, pr$sample_b)); k2 <- match(paste(long$a, long$b), paste(pr$sample_b, pr$sample_a)); k[is.na(k)] <- k2[is.na(k)]
      long$.tip <- ifelse(is.na(k), sprintf("%s\nrho = %.3f", long$a, long$rho),
        sprintf("%s × %s\nrho = %.3f (95%% CI %.3f–%.3f)\n%s", long$a, long$b, long$rho, pr$ci_low[k], pr$ci_high[k], pr$pair_type[k]))
      display <- libcor_display_names(names, dirname(project$project_root)); labels <- make.unique(qc_label_libraries(display, s), sep = " ")
      long$.tip <- paste0(labels[match(long$a, names)], ifelse(long$a == long$b, "", paste0(" × ", labels[match(long$b, names)])), sub("^[^\n]*", "", long$.tip))
      long$a <- factor(labels[match(long$a, names)], levels = labels); long$b <- factor(labels[match(long$b, names)], levels = rev(labels))
      p <- ggplot2::ggplot(long, ggplot2::aes(x = a, y = b, fill = rho, text = .tip)) + ggplot2::geom_tile(colour = "grey35", linewidth = .5) +
        qc_heatmap_fill(s$heatmap_palette, midpoint = .5, limits = c(min(0, min(long$rho, na.rm = TRUE)), 1)) + ggplot2::labs(x = NULL, y = NULL, fill = "Spearman rho") +
        ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0))
      if (isTRUE(s$show_value_labels) || length(names) <= 12L) p <- p + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = max(2.4, s$row_label_size / 3.4))
      list(plot = qc_heatmap_theme(p, s), data = long)
    },
    gene_heatmap = {
      r <- src$libcor(); g <- r$genes; validate(need(nrow(g) > 0, "This correlation result has no gene-level values; run library correlation again."))
      g <- g[g$library %in% qc_libcor_shown(unique(g$library), r, input, prefix, dirname(project$project_root)), , drop = FALSE]
      validate(need(length(unique(g$library)) >= 2L, "Show at least two libraries."))
      species <- names(sort(table(g$species), decreasing = TRUE))[[1L]]; g <- g[g$species == species, , drop = FALSE]
      wide <- tapply(log10(1 + as.numeric(g$density)), list(g$gene, g$library), mean); wide <- wide[stats::complete.cases(wide), , drop = FALSE]
      n <- max(20L, min(as.integer(input[[paste0(prefix, "_genes")]] %||% 300L), nrow(wide)))
      top <- wide[order(apply(wide, 1, stats::var), decreasing = TRUE)[seq_len(n)], , drop = FALSE]; centred <- top - rowMeans(top)
      gene_order <- if (nrow(centred) > 2) stats::hclust(stats::dist(centred))$order else seq_len(nrow(centred))
      lib_order <- if (ncol(centred) > 2) stats::hclust(stats::dist(t(centred)))$order else seq_len(ncol(centred))
      long <- as.data.frame(as.table(centred), stringsAsFactors = FALSE); names(long) <- c("gene", "library", "value")
      labels <- make.unique(qc_label_libraries(libcor_display_names(colnames(centred), dirname(project$project_root)), s), sep = " ")
      long$library <- factor(labels[match(long$library, colnames(centred))], levels = labels[lib_order])
      long$gene <- factor(long$gene, levels = rownames(centred)[gene_order]); long$.tip <- sprintf("%s\n%s\ncentred log10 density: %.2f", long$gene, long$library, long$value)
      p <- ggplot2::ggplot(long, ggplot2::aes(x = library, y = gene, fill = value, text = .tip)) + ggplot2::geom_tile(colour = "white", linewidth = if (n <= 80L) .4 else .1) +
        qc_heatmap_fill(s$heatmap_palette, 0) + ggplot2::labs(x = NULL, y = sprintf("%d most variable genes (%s)", n, species), fill = "Density (centred log10)") +
        ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0))
      list(plot = qc_heatmap_theme(p, utils::modifyList(s, list(show_row_labels = isTRUE(s$show_row_labels) && n <= 80))), data = long)
    },
    pair_scatter = {
      r <- src$libcor(); g <- r$genes; validate(need(nrow(g) > 0, "This correlation result has no gene-level values; run library correlation again."))
      a <- input[[paste0(prefix, "_lib_a")]]; b <- input[[paste0(prefix, "_lib_b")]]; req(a, b)
      ga <- g[g$library == a, , drop = FALSE]; gb <- g[g$library == b, , drop = FALSE]
      validate(need(identical(unique(ga$species), unique(gb$species)), "These libraries are from different species; compare them in the correlation heatmap (over orthologs)."))
      d <- merge(ga[, c("gene", "density")], gb[, c("gene", "density")], by = "gene", suffixes = c("_a", "_b"))
      d$density_a <- as.numeric(d$density_a); d$density_b <- as.numeric(d$density_b); rho <- suppressWarnings(stats::cor(d$density_a, d$density_b, method = "spearman"))
      d$.tip <- sprintf("%s\nA: %.2f\nB: %.2f sites per kb", d$gene, d$density_a, d$density_b)
      p <- ggplot2::ggplot(d, ggplot2::aes(x = density_a + .1, y = density_b + .1, text = .tip)) + ggplot2::geom_abline(linetype = 2, colour = "grey60") +
        ggplot2::geom_point(size = s$point_size, alpha = s$point_opacity, colour = "#1A73E8") + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
        ggplot2::labs(x = paste(libcor_display_names(c(a, b), dirname(project$project_root))[1], "(sites per kb + 0.1)"),
                      y = paste(libcor_display_names(c(a, b), dirname(project$project_root))[2], "(sites per kb + 0.1)"),
                      title = sprintf("Spearman rho = %.3f over %d genes", rho, nrow(d)))
      list(plot = qc_apply_display(p, utils::modifyList(s, list(label_angle = 0L)), "x"), data = d)
    },
    validate(need(FALSE, "Choose a plot.")))
}
