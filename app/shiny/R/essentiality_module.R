# Essentiality module: label-free Hermes window HMM (default), plus the hybrid model when curated labels exist
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
    radioButtons("ess_hmm_mode", "HMM", c("Fit the HMM" = "fit", "Reuse an earlier fit" = "reuse"), "fit"),
    conditionalPanel("input.ess_hmm_mode === 'fit'", uiOutput("ess_library_selector"),
      helpText("Choose one untreated library per independent pool, all from one strain. Several pools are fitted jointly.")),
    conditionalPanel("input.ess_hmm_mode === 'reuse'", uiOutput("ess_hmm_reuse_selector")),
    tags$details(class = "ytab-more-options", tags$summary("Settings"),
      numericInput("ess_window", "Window (bp)", 100, min = 50, step = 50),
      numericInput("ess_min_reads", "Reads for a site to count as an insertion", 1, min = 1, step = 1)),
    checkboxInput("ess_use_labels", "Add the hybrid model (needs curated essentiality labels)", FALSE),
    conditionalPanel("input.ess_use_labels", uiOutput("ess_hybrid_inputs")),
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
      tags$p(tags$b("Validation (C. glabrata, Gale et al. 2020). "), "Silver-label AUC: HMM 0.922, hybrid 0.947 (gold 0.974). The hybrid's gain over a random forest retrained without HMM features was small and not significant (+0.003, 95% CI −0.0004 to 0.006). The choice of normalization for the gene features changes the hybrid's AUC by at most 0.002. Genes called essential by both HMM and hybrid are the most reliable; HMM-only calls are candidates for species-specific essentiality or accessibility artefacts."),
      tags$p(tags$b("Inputs for other species. "), "Use the haploid reference the reads were aligned to (for C. albicans the HapA FASTA); the run stops when more than 50% of windows are unmappable or more than 10% of insertion sites are off the reference. Use haploid strains: in a diploid, a heterozygous insertion in an essential gene is viable."),
      ytab_reference("The HMM is adapted from the HMM of TRANSIT (DeJesus & Ioerger 2013, BMC Bioinformatics 14:303)."),
      tags$p("Command: ", tags$code("ytab essentiality"), " (which runs ", tags$code("ytab hermes-hmm"), " and the model table)."))))

essentiality_model_results_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Results",
    uiOutput("ess_result_selector"),
    radioButtons("ess_plot", "Plot", c("Chromosome track" = "track", "Model scores" = "scores", "Calls" = "calls",
                                       "Insertion-target preference" = "target", "ROC against curated labels" = "roc")),
    conditionalPanel("input.ess_plot === 'track'", uiOutput("ess_chromosome_selector")),
    ytab_plot_downloads_ui("ess", c(genes = "Download gene table (TSV)"))),
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

# Model colours: the label-free HMM in teal, the hybrid in brick red, the current random forest in grey.
essentiality_model_colours <- c(`HMM (label-free)` = YTAB_FOCAL, Hybrid = YTAB_ACCENT, `Random forest` = YTAB_COMPARATOR)

ess_ggplot <- function(choice, genes, windows = NULL, chromosome = NULL, manifest = list(), run_dir = NULL) {
  base <- ggplot2::theme_minimal(base_size = 12)
  if (identical(choice, "target")) {
    tm <- ytab_read_analysis_tsv(file.path(run_dir %||% "", "plots", "target_model.tsv"))
    validate(need(nrow(tm) > 0, "This run has no insertion-target summary. Run essentiality again (or `ytab essentiality --refresh-plots RUN`)."))
    used <- isTRUE(as.logical(tm$used_in_fit[[1L]]))
    tm$tip <- sprintf("offset %+d<br>%.2f bits<br>A %.2f  C %.2f  G %.2f  T %.2f", tm$offset, tm$information_bits, tm$A, tm$C, tm$G, tm$T)
    return(ggplot2::ggplot(tm, ggplot2::aes(offset, information_bits, text = tip)) + ggplot2::geom_col(fill = YTAB_FOCAL, width = .8) +
      ggplot2::geom_text(data = tm[tm$information_bits >= .05, , drop = FALSE], ggplot2::aes(label = preferred), vjust = -.4, fontface = "bold", size = 4.5) +
      ggplot2::scale_x_continuous(breaks = seq(min(tm$offset), max(tm$offset), 2)) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .25))) +
      ggplot2::labs(x = "Offset from the insertion site (0 = first base of the target-site duplication)", y = "Information (bits)",
                    title = if (used) "Learned from these libraries and used in the HMM fit" else "Learned from these libraries (diagnostic; the fit used the Hermes motif)") + base)
  }
  if (identical(choice, "roc")) {
    roc <- ytab_read_analysis_tsv(file.path(run_dir %||% "", "plots", "roc.tsv"))
    validate(need(nrow(roc) > 0, "ROC curves need curated labels: run essentiality with the hybrid model (labels), or refresh an earlier run's plots."))
    roc$tier <- factor(ifelse(roc$tier == "all", "All labelled genes", paste(tools::toTitleCase(roc$tier), "labels")))
    roc$tier <- factor(roc$tier, unique(c("All labelled genes", levels(roc$tier))))
    auc <- unique(roc[, c("model", "tier", "genes", "essential", "auc")])
    auc$label <- sprintf("%s %.3f", sub(" \\(label-free\\)", "", auc$model), auc$auc)
    auc$y <- .05 + .08 * (stats::ave(seq_len(nrow(auc)), auc$tier, FUN = seq_along) - 1)
    roc$tip <- sprintf("%s<br>%s<br>FPR %.3f, TPR %.3f<br>AUC %.3f (%s genes, %s essential)", roc$model, roc$tier, roc$fpr, roc$tpr, roc$auc,
                       format(roc$genes, big.mark = ","), format(roc$essential, big.mark = ","))
    return(ggplot2::ggplot(roc, ggplot2::aes(fpr, tpr, colour = model, group = model, text = tip)) +
      ggplot2::geom_abline(linetype = "dotted", colour = "grey50") + ggplot2::geom_step(linewidth = .9) +
      ggplot2::geom_text(data = auc, ggplot2::aes(x = .98, y = y, label = label, colour = model), hjust = 1, fontface = "bold", size = 3.6, inherit.aes = FALSE, show.legend = FALSE) +
      ggplot2::facet_wrap(~tier, nrow = 1) + ggplot2::scale_colour_manual(values = essentiality_model_colours) + ggplot2::coord_equal() +
      ggplot2::scale_x_continuous(breaks = c(0, .5, 1)) + ggplot2::scale_y_continuous(breaks = c(0, .5, 1)) +
      ggplot2::labs(x = "False positive rate", y = "True positive rate", colour = NULL) + base + ggplot2::theme(legend.position = "top"))
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
    if ("hybrid_score" %in% names(genes) && any(!is.na(genes$hybrid_score))) {
      d <- genes[!is.na(genes$hybrid_score), , drop = FALSE]
      return(ggplot2::ggplot(d, ggplot2::aes(hmm_score, hybrid_score, colour = agreement,
          text = sprintf("%s %s<br>HMM %.2f (%s)<br>hybrid %.2f<br>%s", ORF, ifelse(is.na(Name), "", Name), hmm_score, hmm_call, hybrid_score, agreement))) +
        ggplot2::geom_point(alpha = .6, size = 1.4) + ggplot2::geom_hline(yintercept = .5, linetype = "dashed") + ggplot2::geom_vline(xintercept = .5, linetype = "dashed") +
        ggplot2::scale_colour_manual(values = essentiality_agreement_colours, labels = essentiality_agreement_labels) +
        ggplot2::labs(x = "HMM score (P_ES_run)", y = "Hybrid score", colour = "Agreement") + base)
    }
    d <- genes[!is.na(genes$hmm_score), , drop = FALSE]
    return(ggplot2::ggplot(d, ggplot2::aes(hmm_score, fill = hmm_call)) + ggplot2::geom_histogram(bins = 50) +
      ggplot2::scale_fill_manual(values = c(essentiality_state_colours, uncertain = "grey70", `N/A` = "grey90")) +
      ggplot2::labs(x = "HMM score (P_ES_run)", y = "Genes", fill = "HMM call") + base)
  }
  # counts as written by the pipeline (manifest calls / agreement)
  counts <- function(x, model) if (length(x)) data.frame(call = names(x), Freq = as.integer(unlist(x)), model = model, stringsAsFactors = FALSE)
  calls <- rbind(counts(manifest$calls, "HMM"), counts(manifest$agreement, "HMM vs hybrid"))
  if (is.null(calls) || !nrow(calls)) return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "Call counts are not in this run's manifest; run it again.") + ggplot2::theme_void())
  calls$label <- ifelse(calls$call %in% names(essentiality_state_labels), essentiality_state_labels[calls$call],
                        ifelse(calls$call %in% names(essentiality_agreement_labels), essentiality_agreement_labels[calls$call], calls$call))
  fills <- c(essentiality_state_colours, essentiality_agreement_colours, uncertain = "grey70")
  ggplot2::ggplot(calls, ggplot2::aes(stats::reorder(call, -Freq), Freq, fill = call, text = sprintf("%s: %s genes", label, format(Freq, big.mark = ",")))) +
    ggplot2::geom_col(show.legend = FALSE) + ggplot2::scale_fill_manual(values = fills) +
    ggplot2::scale_x_discrete(labels = function(x) ifelse(x %in% names(essentiality_agreement_labels), essentiality_agreement_labels[x], x)) +
    ggplot2::facet_wrap(~model, scales = "free_x") +
    ggplot2::labs(x = NULL, y = "Genes") + base + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
}

essentiality_module_server <- function(input, output, session, repo_root, active, active_project_path, jobs, python_bin) {
  runs_dir <- reactive(file.path(active()$project_root, "essentiality"))
  runs <- reactivePoll(3000, session,
    checkFunc = function() { x <- ytab_analysis_runs(tryCatch(runs_dir(), error = function(e) "")); paste(x$path, x$modified, collapse = "|") },
    valueFunc = function() { x <- ytab_analysis_runs(tryCatch(runs_dir(), error = function(e) "")); x[file.exists(file.path(x$path, "models", "genes.tsv")), , drop = FALSE] })
  inventory <- reactive(essentiality_parent_inventory(active()))
  # Normalizations the hybrid's gene features can use: the methods of the project's Site Normalization runs (MidLC
  # common saturation among them). BGC first when available (the app's suggested normalization).
  hybrid_normalizations <- reactive({
    root <- active()$project_root
    methods <- unique(unlist(lapply(ytab_analysis_runs(file.path(root, "normalization"))$path, function(run) ytab_read_manifest(run)$methods %||% character())))
    methods <- setdiff(as.character(methods), "nonorm"); methods <- c(intersect("betageom", methods), setdiff(methods, "betageom"))
    choices <- stats::setNames(methods, as.character(norm_method_label(methods)))
    choices
  })
  label_files <- reactive(list.files(file.path(repo_root, "resources", "species"), pattern = "^labels\\.tsv$", recursive = TRUE, full.names = TRUE))

  output$ess_library_selector <- renderUI({
    x <- inventory(); x <- x[x$Included & !x$Explicit_treated, , drop = FALSE]
    choices <- stats::setNames(x$Sample, sprintf("%s (%s, pool %s)", ytab_sample_labels(x$Sample, active()), x$Background, x$Pool))
    selectizeInput("ess_libraries", "Libraries", choices, selected = essentiality_default_libraries(inventory()), multiple = TRUE)
  })
  output$ess_hmm_reuse_selector <- renderUI({
    x <- runs()
    if (!nrow(x)) return(tags$p(class = "text-muted", "No earlier HMM fit in this project."))
    selectInput("ess_hmm_reuse", "Earlier fit", ytab_run_choices(x))
  })
  output$ess_hybrid_inputs <- renderUI({
    norms <- hybrid_normalizations(); labels <- label_files()
    tagList(
      selectInput("ess_labels", "Curated labels", c(stats::setNames(labels, sub("/essentiality_labels/labels.tsv$", "", sub(".*/resources/species/", "", labels))), "Other file" = "other"),
                  selected = labels[grepl(paste0("/species/", tolower(active()$species %||% ""), "/"), labels, fixed = TRUE)][1] %||% labels[1]),
      conditionalPanel("input.ess_labels === 'other'", textInput("ess_labels_path", "Label file (TSV: gene_id, label, tier)")),
      if (length(norms)) selectInput("ess_normalization", "Normalization for the gene features", norms)
      else tags$p(class = "text-warning", "No normalization in this project yet. Run Normalization > Site Normalization (or MidLC) first."))
  })

  observeEvent(input$ess_run, {
    if (jobs$job_is_running()) { showNotification("Wait for the current job to finish.", type = "warning"); return() }
    reuse <- identical(input$ess_hmm_mode, "reuse")
    hmm_args <- if (reuse) {
      run <- input$ess_hmm_reuse %||% ""
      if (!nzchar(run)) { showNotification("Choose an earlier fit.", type = "error"); return() }
      c("--hmm-run", ytab_read_manifest(run)$hmm_run)
    } else {
      libraries <- input$ess_libraries %||% character()
      if (!length(libraries)) { showNotification("Choose at least one library.", type = "error"); return() }
      c("--libraries", libraries, "--window", as.character(input$ess_window %||% 100), "--min-reads", as.character(input$ess_min_reads %||% 1))
    }
    hybrid_args <- character()
    if (isTRUE(input$ess_use_labels)) {
      labels <- if (identical(input$ess_labels, "other")) trimws(input$ess_labels_path %||% "") else input$ess_labels %||% ""
      if (!nzchar(input$ess_normalization %||% "") || !nzchar(labels) || !file.exists(labels)) {
        showNotification("The hybrid needs a label file and a normalization (Normalization tab).", type = "error"); return()
      }
      hybrid_args <- c("--labels", labels, "--normalization", input$ess_normalization)
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
    if (!nrow(x)) return(ytab_empty_state("No essentiality result yet", "The HMM needs only the hit files of untreated parent libraries. Choose the libraries and run."))
    m <- ytab_read_manifest(x$path[[1L]])
    models <- if (isTRUE(m$hybrid)) "HMM and hybrid" else "HMM"
    ytab_result_card("Latest result", tags$p(sprintf("%s: %s on %s.", names(ytab_run_choices(x[1L, , drop = FALSE])), models, paste(ytab_sample_labels(unlist(m$libraries), active()), collapse = ", "))),
                     tags$p("See the Results tab."))
  })

  selected_run <- reactive({ req(nrow(runs()) > 0); path <- input$ess_result %||% runs()$path[[1L]]; req(path %in% runs()$path); path })
  manifest <- reactive(ytab_read_manifest(selected_run()))
  genes <- reactive(ytab_read_analysis_tsv(file.path(selected_run(), "models", "genes.tsv")))
  windows <- reactive(ess_read_windows(selected_run()))

  output$ess_result_selector <- renderUI({
    x <- runs(); if (!nrow(x)) return(tags$p(class = "text-muted", "No essentiality result yet."))
    selectInput("ess_result", "Result", ytab_run_choices(x))
  })
  output$ess_chromosome_selector <- renderUI({ w <- windows(); req(nrow(w) > 0); chroms <- unique(w$chromosome)
    selectInput("ess_chromosome", "Chromosome", stats::setNames(chroms, qc_plot_chromosome_display(chroms, active()))) })
  output$ess_summary <- renderUI({
    if (!nrow(runs())) return(ytab_empty_state("No essentiality result yet", "Run essentiality in the Run tab."))
    m <- manifest(); req(m); calls <- m$hmm$gene_calls
    cards <- lapply(names(essentiality_state_labels), function(s) ytab_metric_card(essentiality_state_labels[[s]], format(calls[[s]] %||% 0, big.mark = ","), "genes (HMM)"))
    if (isTRUE(m$hybrid)) cards <- c(cards, lapply(c("consensus essential", "HMM only", "hybrid only"), function(a)
      ytab_metric_card(essentiality_agreement_labels[[a]], format(m$agreement[[a]] %||% 0, big.mark = ","), "genes (HMM vs hybrid)")))
    notes <- unlist(m$hmm$notes %||% list())
    tagList(tags$div(class = "ytab-stat-grid", cards),
      tags$p(class = "text-muted", sprintf("Libraries: %s. Unmappable windows: %.1f%%. %s", paste(ytab_sample_labels(unlist(m$libraries), active()), collapse = ", "), 100 * m$hmm$masked_fraction,
                                           if (isTRUE(m$hmm$converged)) "The HMM converged." else "The HMM did not converge; see the note.")),
      if (length(notes)) ytab_warning_banner("Notes from the HMM", paste(notes, collapse = " ")))
  })
  output$ess_plot_caption <- renderUI(switch(input$ess_plot %||% "track",
    track = "Insertions per window (top), the most probable HMM state of each window (strip; grey = unmappable) and the posterior probabilities of essential and growth-defect states (bottom). Zoom with the mouse; double-click to reset. The PNG/PDF download shows all four state posteriors.",
    scores = if (isTRUE(manifest()$hybrid)) "Each point is a gene. Dashed lines are the 0.5 thresholds used for both models." else "Distribution of the HMM score (P_ES_run) by call.",
    calls = "Number of genes per HMM call, and agreement between the HMM and the hybrid when the hybrid was run.",
    target = "Information (bits) at each position around the insertion site, from the libraries' distinct insertion sites against the local base composition, with the preferred base above informative positions. No labels or prior motif are used; Hermes shows T at +1 and A at +6 (nTnnnnAn).",
    roc = "Each model's score against the curated labels, per label tier and all together, with the area under the curve (AUC). Only genes a model was not fitted to: the HMM uses no labels, and the hybrid's training genes carry out-of-fold scores."))
  output$ess_plot_frame <- renderUI({ req(nrow(runs()) > 0); plotly::plotlyOutput("ess_plotly", height = if (identical(input$ess_plot, "track")) "620px" else "520px") })
  output$ess_plotly <- plotly::renderPlotly({
    choice <- input$ess_plot %||% "track"
    if (identical(choice, "track")) {
      req(input$ess_chromosome); return(ess_track_plotly(windows(), input$ess_chromosome, qc_plot_chromosome_display(input$ess_chromosome, active())))
    }
    out <- ytab_ggplotly(ess_ggplot(choice, genes(), manifest = manifest(), run_dir = selected_run()))
    if (identical(choice, "calls")) {
      # ggplotly draws the facet y title as an annotation that runs into the tick labels: use the axis title instead
      out$x$layout$annotations <- Filter(function(a) !identical(trimws(gsub("<[^>]+>", "", as.character(a$text %||% ""))), "Genes"), out$x$layout$annotations)
      out <- plotly::layout(out, showlegend = FALSE, yaxis = list(automargin = TRUE, title = list(text = "<b>Genes</b>", standoff = 18)))
    }
    out
  })
  output$ess_genes_table <- DT::renderDT({
    g <- genes(); req(nrow(g) > 0)
    keep <- intersect(c("ORF", "Name", "chromosome", "hmm_call", "hmm_score", "hybrid_score", "agreement", "hmm_status", "mappable_fraction", "reference_label"), names(g))
    g <- g[, keep, drop = FALSE]
    for (column in intersect(c("hmm_score", "hybrid_score", "rf_existing_score", "mappable_fraction"), names(g))) g[[column]] <- round(g[[column]], 3)
    if ("chromosome" %in% names(g)) g$chromosome <- qc_plot_chromosome_display(g$chromosome, active())
    if ("agreement" %in% names(g)) g$agreement <- ifelse(g$agreement %in% names(essentiality_agreement_labels), essentiality_agreement_labels[g$agreement], g$agreement)
    labels <- c(ORF = "ORF", Name = "Name", chromosome = "Chromosome", hmm_call = "HMM call", hmm_score = "HMM score", hybrid_score = "Hybrid score",
                agreement = "Agreement", rf_existing_score = "Current RF score", hmm_status = "HMM status", mappable_fraction = "Mappable fraction",
                reference_label = "Curated label")
    names(g) <- labels[names(g)]
    ytab_datatable(g, page_length = 15L)
  })

  plot_for_download <- function() {
    choice <- input$ess_plot %||% "track"
    if (!identical(choice, "track")) return(ess_ggplot(choice, genes(), manifest = manifest(), run_dir = selected_run()))
    p <- ess_ggplot("track", genes(), windows(), input$ess_chromosome %||% unique(windows()$chromosome)[[1L]])
    p$bottom + ggplot2::geom_col(ggplot2::aes(start, insertions / max(insertions, 1)), data = windows()[windows()$chromosome == (input$ess_chromosome %||% ""), ],
                                 inherit.aes = FALSE, fill = "grey20", alpha = .5) +
      ggplot2::labs(caption = "Dark bars: insertions per window, scaled to the largest window")
  }
  download_name <- function(ext) paste0(active()$project_id, "_essentiality_", input$ess_plot %||% "track", "_", basename(selected_run()), ext)
  for (format in c("png", "pdf")) local({
    ext <- format
    output[[paste0("ess_download_", ext)]] <- downloadHandler(filename = function() download_name(paste0(".", ext)),
      content = function(file) ytab_save_ggplot(file, plot_for_download(), width = 12, height = 6))
  })
  output$ess_download_genes <- downloadHandler(filename = function() paste0(active()$project_id, "_essentiality_genes_", basename(selected_run()), ".tsv"),
    content = function(file) file.copy(file.path(selected_run(), "models", "genes.tsv"), file))
  invisible(list(runs = runs, selected_run = selected_run, manifest = manifest, genes = genes, windows = windows, inventory = inventory))
}
