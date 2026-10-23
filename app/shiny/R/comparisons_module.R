# Condition Comparisons module: Treated vs Control (local-abundance z-score, resampling or the control-conditioned
# Gumbel; treated_vs_control.R) plus the Hermes multi-condition tests and genetic interactions. All but the z-score's
# default (total reads, CPM) read a run of the Normalization module (<project>/normalization/<run>/experiment).

cmp_normalization_runs <- function(project_root) {
  x <- ytab_analysis_runs(file.path(project_root, "normalization"))
  x[dir.exists(file.path(x$path, "experiment", "methods")), , drop = FALSE]
}

cmp_run_methods <- function(run) {
  methods <- list.dirs(file.path(run, "experiment", "methods"), full.names = FALSE, recursive = FALSE)
  methods <- methods[file.exists(file.path(run, "experiment", "methods", methods, "normalized_counts.npz")) & methods != "nonorm"]
  methods[order(methods != "betageom")]
}

# The run's samples with their sample-sheet metadata.
cmp_observations <- function(run, samples) {
  obs <- ytab_read_analysis_tsv(file.path(run, "experiment", "matrix", "observations.tsv"))
  if (!nrow(obs) || is.null(samples) || !nrow(samples)) return(obs)
  extra <- samples[, setdiff(names(samples), c(names(obs), "fastq_1", "fastq_2", "include", "warnings", "hit_file")), drop = FALSE]
  extra$sample <- samples$sample
  merged <- merge(obs, extra, by = "sample", all.x = TRUE, sort = FALSE)
  merged[match(obs$sample, merged$sample), , drop = FALSE]
}

cmp_label_columns <- function(obs) {
  skip <- c("observation_id", "sample", "display_name", "project", "dataset", "library_id", "input", "input_sha256", "layout")
  names(obs)[vapply(names(obs), function(n) !n %in% skip && length(unique(stats::na.omit(as.character(obs[[n]][nzchar(as.character(obs[[n]]))])))) >= 2L, logical(1))]
}

cmp_input_selectors <- function(prefix, runs) {
  if (!nrow(runs)) return(tags$p(class = "text-warning", "No normalization run yet. Run Normalization → Site Normalization first."))
  tagList(selectInput(paste0(prefix, "_norm_run"), "Normalization run", ytab_run_choices(runs)), uiOutput(paste0(prefix, "_method_selector")))
}

mc_methods_text <- function() ytab_methods_drawer(
  tags$p(tags$b("Questions. "), "The ANOVA asks whether a gene responds to any condition; all-pairs resampling asks which conditions differ, and in which direction. Both use only the insertion sites the gene has in the compared samples, after normalization (BGC suggested) and iterative-median recentering."),
  tags$p(tags$b("ANOVA. "), "A permutation F-test on the gene's normalized site counts, with staged permutations (up to 99,999) and Benjamini–Hochberg FDR. When one condition's replicates are much noisier than the others, ANOVA p-values are anti-conservative; the app then shows a warning, and the all-pairs results should be used."),
  tags$p(tags$b("All-pairs. "), "Site-restricted resampling for every pair of conditions, with one FDR family over all pairs (call_all_pairs). Several libraries (pools) each split into the conditions are permuted within each library (stratified)."),
  tags$p(tags$b("Plots. "), "Response correlation: Spearman correlation between conditions' responses (log2FC versus the reference), with a split reference so the shared reference does not inflate it; it needs at least two non-reference conditions. Response heatmap: log2FC of responsive genes, clustered. Correlations on responsive genes are descriptive."),
  tags$p(tags$b("Validation. "), "Pre-registered simulation (benchmark v2): type I error and FDR within their gates for both tests; rank-based U-test and ANOVA alternatives had almost no power on Hermes data."),
  ytab_reference("Adapted from the ANOVA, resampling, correlation plot and heatmap of TRANSIT (DeJesus et al. 2015, PLoS Comput Biol 11:e1004401)."),
  tags$p("Command: ", tags$code("ytab multicondition-project"), " (runs ", tags$code("ytab multicondition"), " and ", tags$code("ytab multicondition-plots"), ")."))

gi_methods_text <- function() ytab_methods_drawer(
  tags$p(tags$b("Question. "), "Does a gene's treatment response differ between two strains? Δ = log2FC(query strain B) − log2FC(reference strain A). Classes: aggravating (Δ < 0), alleviating (Δ ≥ 0 and |B| < |A|), suppressive (Δ ≥ 0 and |B| > |A|)."),
  tags$p(tags$b("Design. "), "Each strain needs at least two control/treated pairs: replicate cultures split from one pool, or several pools each split into control and treated. Control and treated cultures from one pool share clones, so comparing them site by site cancels each clone's founder abundance. When control and treated do not share clones (read correlation at shared sites below 0.5), the app warns that power is low."),
  tags$p(tags$b("Test. "), "Each control/treated pair's response is the mean over the gene's sites of log2((treated + 5)/(control + 5)). Its variance combines the replicate variance with a site-based variance, scaled by factors fitted across genes. Welch t-test with Satterthwaite degrees of freedom; Benjamini–Hochberg FDR."),
  tags$p(tags$b("Options. "), "Interval for Δ: validated, coverage 94.3–97.2%. ROPE (e.g. ±0.5): tests |Δ| > ROPE (McCarthy & Smyth 2009); validated but strict at 0.5. Holm family-wise correction: not validated (8.3% family-wise error against a 7.5% gate); use BH."),
  tags$p(tags$b("Validation. "), "Pre-registered simulation: type I error 3.1–5.6%, FDR ≤ 9.6%. A posterior-overlap test without a model of library noise made almost no calls on Hermes data."),
  ytab_reference("Adapted from the genetic-interaction method of TRANSIT (DeJesus et al. 2017, Nucleic Acids Res 45:e93)."),
  tags$p("Command: ", tags$code("ytab genetic-interactions"), "."))

comparisons_mc_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Multi-condition",
    uiOutput("mc_inputs"),
    uiOutput("mc_condition_controls"),
    tags$details(class = "ytab-more-options", tags$summary("Settings"),
      radioButtons("mc_test", "Tests", c("ANOVA and all-pairs" = "both", "ANOVA" = "anova", "All-pairs" = "all-pairs"), "both"),
      numericInput("mc_fdr", "FDR", .05, min = .001, max = .25, step = .01)),
    actionButton("mc_run", "Run multi-condition tests", class = "btn-primary"),
    tags$hr(),
    uiOutput("mc_result_selector"),
    radioButtons("mc_plot", "Plot", c("Volcano" = "volcano", "Response correlation" = "correlation", "Response heatmap" = "heatmap")),
    uiOutput("mc_pair_selector"),
    ytab_plot_downloads_ui("mc", c(anova = "Download ANOVA (TSV)", pairs = "Download all-pairs (TSV)"))),
  main = tagList(uiOutput("mc_summary"),
    ytab_plot_card("Multi-condition", uiOutput("mc_plot_frame"), description = uiOutput("mc_plot_caption")),
    ytab_table_card("Genes", DT::DTOutput("mc_table")), mc_methods_text()))

comparisons_gi_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Genetic interactions",
    uiOutput("gi_inputs"),
    uiOutput("gi_strain_controls"),
    tags$details(class = "ytab-more-options", tags$summary("Settings"),
      numericInput("gi_fdr", "FDR", .05, min = .001, max = .25, step = .01),
      numericInput("gi_rope", "ROPE: call only |Δ| demonstrably above (0 = any Δ ≠ 0; 0.5 is a common choice)", 0, min = 0, step = .1),
      selectInput("gi_correction", "Multiple-testing correction", c("Benjamini–Hochberg FDR (validated)" = "bh", "Holm family-wise (not validated)" = "holm")),
      numericInput("gi_confidence", "Interval level for Δ", .95, min = .5, max = .999, step = .01)),
    actionButton("gi_run", "Run genetic interactions", class = "btn-primary"),
    tags$hr(),
    uiOutput("gi_result_selector"),
    radioButtons("gi_plot", "Plot", c("Strain responses" = "scatter", "Interaction effects" = "forest", "Volcano" = "volcano")),
    ytab_plot_downloads_ui("gi", c(genes = "Download gene table (TSV)"))),
  main = tagList(uiOutput("gi_summary"),
    ytab_plot_card("Genetic interactions", uiOutput("gi_plot_frame"), description = uiOutput("gi_plot_caption")),
    ytab_table_card("Genes", DT::DTOutput("gi_table")), gi_methods_text()))

comparisons_module_ui <- function() navset_tab(id = "fitness_tabs",
  nav_panel("Design", value = "fitness_design", fitness_design_ui()),
  nav_panel("Treated vs Control: Run", value = "fitness_run", tvc_run_tab_ui()),
  nav_panel("Treated vs Control: Results", value = "fitness_results", tvc_results_tab_ui()),
  nav_panel("Multi-condition", value = "multicondition", comparisons_mc_ui()),
  nav_panel("Genetic Interactions", value = "gi", comparisons_gi_ui()))

mc_ggplot <- function(choice, data, pair = NULL) {
  base <- ggplot2::theme_minimal(base_size = 12)
  if (identical(choice, "correlation")) {
    d <- data$correlation; d <- d[d$gene_set == "all tested genes", , drop = FALSE]
    conditions <- unique(c(d$condition_a, d$condition_b))
    both <- rbind(d, transform(d, condition_a = condition_b, condition_b = condition_a),
                  data.frame(gene_set = d$gene_set[1], condition_a = conditions, condition_b = conditions, genes = NA, method = NA, rho = 1, ci_low = NA, ci_high = NA,
                             rho_shared_reference = NA, reference_splits = NA))
    return(ggplot2::ggplot(both, ggplot2::aes(condition_a, condition_b, fill = rho, text = sprintf("%s vs %s<br>rho %.2f (95%% CI %.2f–%.2f)<br>%s genes", condition_a, condition_b, rho, ci_low, ci_high, genes))) +
      ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho))) +
      ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1)) +
      ggplot2::labs(x = NULL, y = NULL, fill = "Spearman rho") + base)
  }
  if (identical(choice, "heatmap")) {
    m <- data$heatmap; value_columns <- grep("^log2FC_", names(m), value = TRUE)
    long <- do.call(rbind, lapply(value_columns, function(v) data.frame(gene = ifelse(is.na(m$Name) | !nzchar(m$Name), m$ORF, m$Name), ORF = m$ORF,
                                                                      condition = sub("^log2FC_", "", v), log2FC = m[[v]])))
    long$gene <- factor(long$gene, levels = rev(unique(long$gene)))
    limit <- max(abs(long$log2FC), na.rm = TRUE)
    return(ggplot2::ggplot(long, ggplot2::aes(condition, gene, fill = log2FC, text = sprintf("%s (%s)<br>%s: log2FC %.2f", gene, ORF, condition, log2FC))) +
      ggplot2::geom_tile(colour = "white", linewidth = .2) + ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-limit, limit)) +
      ggplot2::labs(x = NULL, y = NULL, fill = "log2FC vs reference") + base + ggplot2::theme(axis.text.y = ggplot2::element_blank()))
  }
  d <- data$pairs[data$pairs$pair == pair & data$pairs$status == "tested", , drop = FALSE]
  d$call_all_pairs <- factor(d$call_all_pairs, levels = c("depleted", "enriched", "not significant"))
  ggplot2::ggplot(d, ggplot2::aes(log2FC, -log10(pmax(p_value, 1e-300)), colour = call_all_pairs,
      text = sprintf("%s %s<br>log2FC %.2f<br>p %.3g; q (all pairs) %.3g<br>%d sites", ORF, ifelse(is.na(Name), "", Name), log2FC, p_value, q_all_pairs, n_sites))) +
    ggplot2::geom_point(alpha = .6, size = 1.3) + ggplot2::scale_colour_manual(values = c(depleted = "#b2182b", enriched = "#2166ac", `not significant` = "grey70"), drop = FALSE) +
    ggplot2::labs(x = sprintf("log2FC (%s)", pair), y = "-log10 p", colour = "Call (all pairs, one FDR family)") + base
}

gi_ggplot <- function(choice, genes, reference, query) {
  base <- ggplot2::theme_minimal(base_size = 12)
  d <- genes[genes$status == "tested", , drop = FALSE]
  d$call <- factor(d$call, levels = c("Aggravating", "Alleviating", "Suppressive", "No interaction"))
  colours <- c(Aggravating = "#b2182b", Alleviating = "#2166ac", Suppressive = "#1a9850", `No interaction` = "grey75")
  hover <- function(x) sprintf("%s %s<br>log2FC %s %.2f; %s %.2f<br>Δ %.2f%s<br>q %.3g; %s", x$ORF, ifelse(is.na(x$Name), "", x$Name), reference, x$log2FC_reference,
                               query, x$log2FC_query, x$delta_log2FC,
                               if ("delta_lower" %in% names(x)) sprintf(" (%.2f to %.2f)", x$delta_lower, x$delta_upper) else "", x[["p_adjusted"]] %||% x$q_bh, x$call)
  if (identical(choice, "forest")) {
    q <- d[["p_adjusted"]] %||% d$q_bh
    top <- utils::head(d[order(q, d$p_value), , drop = FALSE], 30L)
    top$gene <- factor(ifelse(is.na(top$Name) | !nzchar(top$Name), top$ORF, top$Name), levels = rev(unique(ifelse(is.na(top$Name) | !nzchar(top$Name), top$ORF, top$Name))))
    plot <- ggplot2::ggplot(top, ggplot2::aes(delta_log2FC, gene, colour = call, text = hover(top))) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::geom_point(size = 2.5) + ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
      ggplot2::labs(x = sprintf("Delta log2FC (%s - %s)", query, reference), y = NULL, colour = "Call", title = "30 genes with the smallest adjusted p-values") + base
    if ("delta_lower" %in% names(top)) plot <- plot + ggplot2::geom_errorbar(ggplot2::aes(xmin = delta_lower, xmax = delta_upper), orientation = "y", width = .3)
    return(plot)
  }
  if (identical(choice, "volcano"))
    return(ggplot2::ggplot(d, ggplot2::aes(delta_log2FC, -log10(pmax(p_value, 1e-300)), colour = call, text = hover(d))) + ggplot2::geom_point(alpha = .6, size = 1.3) +
      ggplot2::scale_colour_manual(values = colours, drop = FALSE) + ggplot2::labs(x = sprintf("Delta log2FC (%s - %s)", query, reference), y = "-log10 p", colour = "Call") + base)
  ggplot2::ggplot(d, ggplot2::aes(log2FC_reference, log2FC_query, colour = call, text = hover(d))) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") + ggplot2::geom_point(alpha = .6, size = 1.3) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::labs(x = sprintf("Treatment response in %s (log2FC)", reference), y = sprintf("Treatment response in %s (log2FC)", query), colour = "Call") + base
}

comparisons_server <- function(input, output, session, repo_root, active, active_project_path, jobs, python_bin) {
  norm_runs <- reactivePoll(4000, session,
    checkFunc = function() { x <- tryCatch(cmp_normalization_runs(active()$project_root), error = function(e) data.frame()); paste(x$path, collapse = "|") },
    valueFunc = function() tryCatch(cmp_normalization_runs(active()$project_root), error = function(e) ytab_analysis_runs("")))
  analysis_runs <- function(folder, marker) reactivePoll(3000, session,
    checkFunc = function() { x <- ytab_analysis_runs(tryCatch(file.path(active()$project_root, folder), error = function(e) ""), marker); paste(x$path, x$modified, collapse = "|") },
    valueFunc = function() ytab_analysis_runs(tryCatch(file.path(active()$project_root, folder), error = function(e) ""), marker))
  mc_runs <- analysis_runs("multicondition", "manifest.json"); gi_runs <- analysis_runs("genetic_interactions", "run_manifest.json")
  start <- function(script, args, stage, message) {
    if (jobs$job_is_running()) { showNotification("Wait for the current job to finish.", type = "warning"); return(invisible(FALSE)) }
    tryCatch({ jobs$start_job(python_bin(), c(file.path(repo_root, "scripts", "local", script), args), active_project_path(), stage, wd = repo_root)
               showNotification(message, type = "message") }, error = function(e) showNotification(conditionMessage(e), type = "error"))
  }
  tvc_server(input, output, session, repo_root, active, active_project_path, norm_runs, start)
  reference_files <- reactive({ cfg <- yaml::read_yaml(active_project_path()); resolve <- function(x) if (grepl("^/", x)) x else file.path(repo_root, x)
                                list(fasta = resolve(cfg$reference$fasta), gff = resolve(cfg$reference$gff)) })

  # Shared: normalization run and method per analysis.
  for (prefix in c("mc", "gi")) local({
    pre <- prefix
    output[[paste0(pre, "_inputs")]] <- renderUI(cmp_input_selectors(pre, norm_runs()))
    output[[paste0(pre, "_method_selector")]] <- renderUI({
      run <- input[[paste0(pre, "_norm_run")]]; req(run); methods <- cmp_run_methods(run)
      selectInput(paste0(pre, "_method"), "Normalization method", stats::setNames(methods, as.character(norm_method_label(methods))))
    })
  })
  mc_obs <- reactive({ req(input$mc_norm_run); cmp_observations(input$mc_norm_run, active()$samples) })
  gi_obs <- reactive({ req(input$gi_norm_run); cmp_observations(input$gi_norm_run, active()$samples) })

  # ---- Multi-condition -------------------------------------------------------------------------------------------
  output$mc_condition_controls <- renderUI({
    obs <- mc_obs(); req(nrow(obs) > 0); columns <- cmp_label_columns(obs)
    default <- intersect(c("treatment", "condition", "condition_label", "role"), columns)[1] %||% columns[1] %||% "role"
    tagList(selectInput("mc_condition_column", "Condition from sample sheet column", unique(c(columns, "role")), default), uiOutput("mc_condition_table"))
  })
  output$mc_condition_table <- renderUI({
    obs <- mc_obs(); column <- input$mc_condition_column; req(nrow(obs) > 0, column %in% names(obs))
    values <- as.character(obs[[column]]); values[is.na(values)] <- ""
    tagList(tags$p(class = "text-muted", "Edit a condition to regroup samples; leave it empty to leave a sample out."),
      lapply(seq_len(nrow(obs)), function(i) textInput(paste0("mc_cond_", i), obs$display_name[[i]] %||% obs$sample[[i]], values[[i]])),
      uiOutput("mc_reference_selector"),
      checkboxInput("mc_stratify", "Several pools, each split into the conditions (permute within pools)", length(unique(obs$library_id)) > 1L))
  })
  mc_assignment <- reactive({
    obs <- mc_obs(); req(nrow(obs) > 0)
    conditions <- vapply(seq_len(nrow(obs)), function(i) trimws(input[[paste0("mc_cond_", i)]] %||% ""), "")
    data.frame(observation_id = obs$observation_id, condition = conditions, role = obs$role, stringsAsFactors = FALSE)[nzchar(conditions), , drop = FALSE]
  })
  output$mc_reference_selector <- renderUI({
    a <- mc_assignment(); req(nrow(a) > 0); controls <- a$condition[a$role == "control"]
    default <- if (length(controls)) names(sort(table(controls), decreasing = TRUE))[[1L]] else unique(a$condition)[[1L]]
    selectInput("mc_reference", "Reference condition", unique(a$condition), default)
  })
  observeEvent(input$mc_run, {
    a <- mc_assignment(); ref <- reference_files()
    if (length(unique(a$condition)) < 2L) { showNotification("Assign at least two conditions.", type = "error"); return() }
    out <- ytab_new_run_dir(file.path(active()$project_root, "multicondition"), "multicondition")
    args <- c("--project-config", active_project_path(), "--normalized", file.path(input$mc_norm_run, "experiment", "methods"),
              "--normalization", input$mc_method %||% "betageom", as.vector(rbind("--sample", paste0(a$observation_id, "=", a$condition))),
              "--reference-condition", input$mc_reference %||% a$condition[[1L]], "--test", input$mc_test %||% "both",
              "--fdr", as.character(input$mc_fdr %||% .05), if (isTRUE(input$mc_stratify)) "--stratify-libraries", "--output", out)
    start("ytab_multicondition_project.py", args, "multicondition", "Multi-condition tests started; permutation tests can take several minutes.")
  }, ignoreInit = TRUE)

  mc_selected <- reactive({ x <- mc_runs(); req(nrow(x) > 0); path <- input$mc_result %||% x$path[[1L]]; req(path %in% x$path); path })
  mc_data <- reactive({ run <- mc_selected()
    list(manifest = ytab_read_manifest(run), anova = ytab_read_analysis_tsv(file.path(run, "tests", "anova_genes.tsv")),
         pairs = ytab_read_analysis_tsv(file.path(run, "tests", "pairs_genes.tsv")), correlation = ytab_read_analysis_tsv(file.path(run, "plots", "response_correlation.tsv")),
         heatmap = ytab_read_analysis_tsv(file.path(run, "plots", "heatmap_matrix.tsv"))) })
  output$mc_result_selector <- renderUI({ x <- mc_runs(); if (!nrow(x)) return(tags$p(class = "text-muted", "No multi-condition result yet.")); selectInput("mc_result", "Result", ytab_run_choices(x)) })
  output$mc_pair_selector <- renderUI({ req(identical(input$mc_plot, "volcano")); p <- mc_data()$pairs; req(nrow(p) > 0); selectInput("mc_pair", "Pair", unique(p$pair)) })
  output$mc_summary <- renderUI({
    if (!nrow(mc_runs())) return(ytab_empty_state("No multi-condition result yet", "Choose a normalization run, assign conditions and run the tests."))
    m <- mc_data()$manifest; req(m)
    cards <- c(if (!is.null(m$anova_calls)) list(ytab_metric_card("ANOVA: genes that vary", format(m$anova_calls$varies %||% 0, big.mark = ","), sprintf("of %s tested", format(sum(unlist(m$anova_calls)) - (m$anova_calls[["not tested"]] %||% 0), big.mark = ",")))),
               lapply(names(m$all_pairs_calls), function(pair) { x <- m$all_pairs_calls[[pair]]; ytab_metric_card(pair, sprintf("%s ↓ / %s ↑", x$depleted %||% 0, x$enriched %||% 0), "depleted / enriched (all-pairs)") }))
    tagList(tags$div(class = "ytab-stat-grid", cards),
      tags$p(class = "text-muted", sprintf("Conditions: %s (reference %s); normalization %s%s.", paste(unlist(m$conditions), collapse = ", "), m$reference_condition,
                                           m$normalization, if (isTRUE(m$stratified)) "; permuted within pools" else "")),
      if (!is.null(m$anova_warning)) ytab_warning_banner("Use the all-pairs results", paste0(m$anova_warning, ". ANOVA p-values are anti-conservative here.")),
      if (!identical(m$plots$status, "success") || length(m$plots$notes)) tags$p(class = "text-muted", paste(c(m$plots$message, unlist(m$plots$notes)), collapse = " ")))
  })
  output$mc_plot_caption <- renderUI(switch(input$mc_plot %||% "volcano",
    volcano = "Each point is a gene; colours are calls with one FDR family over all pairs.",
    correlation = "Spearman correlation between conditions' responses (log2FC versus the reference) over all tested genes. Needs at least two non-reference conditions.",
    heatmap = "log2FC versus the reference for responsive genes, clustered on their profiles."))
  output$mc_plot_frame <- renderUI({
    req(nrow(mc_runs()) > 0); d <- mc_data(); choice <- input$mc_plot %||% "volcano"
    if (identical(choice, "correlation") && !nrow(d$correlation)) return(tags$p(class = "text-muted", "No response correlation: it needs at least two non-reference conditions."))
    if (identical(choice, "heatmap") && !nrow(d$heatmap)) return(tags$p(class = "text-muted", "No responsive genes to show."))
    plotly::plotlyOutput("mc_plotly", height = if (identical(choice, "heatmap")) "700px" else "520px")
  })
  mc_current_plot <- function() { choice <- input$mc_plot %||% "volcano"; mc_ggplot(choice, mc_data(), input$mc_pair %||% unique(mc_data()$pairs$pair)[[1L]]) }
  output$mc_plotly <- plotly::renderPlotly(ytab_ggplotly(mc_current_plot()))
  output$mc_table <- DT::renderDT({
    d <- mc_data(); p <- d$pairs; req(nrow(p) > 0); pair <- input$mc_pair %||% unique(p$pair)[[1L]]
    p <- p[p$pair == pair & p$status == "tested", intersect(c("ORF", "Name", "n_sites", "log2FC", "p_value", "q_all_pairs", "call_all_pairs", "mean_control", "mean_treatment"), names(p)), drop = FALSE]
    if (nrow(d$anova)) p <- merge(p, d$anova[, c("ORF", "q_bh", "call")], by = "ORF", all.x = TRUE, suffixes = c("", "_anova"))
    for (column in intersect(c("log2FC", "p_value", "q_all_pairs", "mean_control", "mean_treatment", "q_bh"), names(p))) p[[column]] <- signif(p[[column]], 3)
    names(p)[names(p) == "q_bh"] <- "q_anova"; names(p)[names(p) == "call"] <- "call_anova"
    ytab_datatable(p[order(p$q_all_pairs), , drop = FALSE], page_length = 15L)
  })
  mc_name <- function(ext) paste0(active()$project_id, "_multicondition_", input$mc_plot %||% "volcano", "_", basename(mc_selected()), ext)
  for (format in c("png", "pdf")) local({ ext <- format
    output[[paste0("mc_download_", ext)]] <- downloadHandler(filename = function() mc_name(paste0(".", ext)), content = function(file) ytab_save_ggplot(file, mc_current_plot(), width = 10, height = 8)) })
  output$mc_download_anova <- downloadHandler(filename = function() mc_name("_anova.tsv"), content = function(file) file.copy(file.path(mc_selected(), "tests", "anova_genes.tsv"), file))
  output$mc_download_pairs <- downloadHandler(filename = function() mc_name("_pairs.tsv"), content = function(file) file.copy(file.path(mc_selected(), "tests", "pairs_genes.tsv"), file))

  # ---- Genetic interactions --------------------------------------------------------------------------------------
  output$gi_strain_controls <- renderUI({
    obs <- gi_obs(); req(nrow(obs) > 0); columns <- cmp_label_columns(obs)
    default <- intersect(c("background", "strain", "guessed_background"), columns)[1] %||% columns[1]
    tagList(selectInput("gi_strain_column", "Strain from sample sheet column", columns, default), uiOutput("gi_strain_pick"))
  })
  gi_assignment <- reactive({
    obs <- gi_obs(); column <- input$gi_strain_column; req(nrow(obs) > 0, column %in% names(obs))
    data.frame(observation_id = obs$observation_id, sample = obs$display_name %||% obs$sample, strain = as.character(obs[[column]]), role = obs$role,
               library = obs$library_id, stringsAsFactors = FALSE)
  })
  output$gi_strain_pick <- renderUI({
    a <- gi_assignment(); strains <- unique(a$strain[nzchar(a$strain)]); req(length(strains) >= 1L)
    tagList(selectInput("gi_reference", "Reference strain (A, e.g. wild type)", strains, strains[[1L]]),
            selectInput("gi_query", "Query strain (B, e.g. the mutant)", strains, strains[[min(2L, length(strains))]]),
            tags$details(class = "ytab-more-options", tags$summary(sprintf("Samples (%d)", nrow(a))), DT::datatable(a[, c("sample", "strain", "role", "library")],
                         rownames = FALSE, options = list(dom = "t", pageLength = 50))))
  })
  observeEvent(input$gi_run, {
    a <- gi_assignment(); ref <- reference_files()
    if (identical(input$gi_reference, input$gi_query)) { showNotification("Choose two different strains.", type = "error"); return() }
    a <- a[a$strain %in% c(input$gi_reference, input$gi_query), , drop = FALSE]
    out <- ytab_new_run_dir(file.path(active()$project_root, "genetic_interactions"), "gi")
    args <- c("--normalized", file.path(input$gi_norm_run, "experiment", "methods"), "--normalization", input$gi_method %||% "betageom",
              "--gff", ref$gff, "--reference", ref$fasta, "--output", out, "--reference-strain", input$gi_reference, "--query-strain", input$gi_query,
              as.vector(rbind("--sample", paste0(a$observation_id, "=", a$strain, ":", a$role))), "--fdr", as.character(input$gi_fdr %||% .05),
              "--rope", as.character(input$gi_rope %||% 0), "--correction", input$gi_correction %||% "bh", "--confidence", as.character(input$gi_confidence %||% .95))
    start("ytab_genetic_interactions.py", args, "genetic_interactions", "Genetic interactions started.")
  }, ignoreInit = TRUE)

  gi_selected <- reactive({ x <- gi_runs(); req(nrow(x) > 0); path <- input$gi_result %||% x$path[[1L]]; req(path %in% x$path); path })
  gi_data <- reactive(list(manifest = read_json_safely(file.path(gi_selected(), "run_manifest.json")), genes = ytab_read_analysis_tsv(file.path(gi_selected(), "genes.tsv"))))
  output$gi_result_selector <- renderUI({ x <- gi_runs(); if (!nrow(x)) return(tags$p(class = "text-muted", "No genetic-interaction result yet.")); selectInput("gi_result", "Result", ytab_run_choices(x)) })
  output$gi_summary <- renderUI({
    if (!nrow(gi_runs())) return(ytab_empty_state("No genetic-interaction result yet", "Choose a normalization run, two strains and run. Each strain needs at least two control/treated pairs."))
    m <- gi_data()$manifest; req(m); d <- m$diagnostics; calls <- d$calls
    cards <- c(list(ytab_metric_card("Genes tested", format(d$tested, big.mark = ","), sprintf("%s vs %s", d$query_strain, d$reference_strain))),
               lapply(c("Aggravating", "Alleviating", "Suppressive"), function(k) ytab_metric_card(k, format(calls[[k]] %||% 0, big.mark = ","), "genes")))
    tagList(tags$div(class = "ytab-stat-grid", cards),
      if (!is.null(d$warning) && nzchar(d$warning)) ytab_warning_banner("Low power: control and treated may not share clones", d$warning),
      if (identical(d$settings$correction %||% "bh", "holm")) ytab_warning_banner("Holm correction is not validated", "Simulated family-wise error was 8.3% against a 5% target; prefer BH."))
  })
  output$gi_plot_caption <- renderUI(switch(input$gi_plot %||% "scatter",
    scatter = "Each gene's treatment response in the two strains. Genes on the dashed diagonal respond equally; coloured genes are called interactions.",
    forest = "Δ with its interval for the 30 genes with the smallest adjusted p-values.",
    volcano = "Δ against the p-value; coloured genes are called interactions."))
  output$gi_plot_frame <- renderUI({ req(nrow(gi_runs()) > 0); plotly::plotlyOutput("gi_plotly", height = if (identical(input$gi_plot, "forest")) "700px" else "520px") })
  gi_current_plot <- function() { d <- gi_data(); gi_ggplot(input$gi_plot %||% "scatter", d$genes, d$manifest$diagnostics$reference_strain, d$manifest$diagnostics$query_strain) }
  output$gi_plotly <- plotly::renderPlotly(ytab_ggplotly(gi_current_plot()))
  output$gi_table <- DT::renderDT({
    g <- gi_data()$genes; req(nrow(g) > 0); g <- g[g$status == "tested", , drop = FALSE]
    keep <- intersect(c("ORF", "Name", "log2FC_reference", "log2FC_query", "delta_log2FC", "delta_lower", "delta_upper", "p_value", "p_adjusted", "q_bh", "call", "interaction", "sites_reference", "sites_query"), names(g))
    g <- g[, keep, drop = FALSE]; for (column in keep[vapply(g, is.numeric, logical(1))]) g[[column]] <- signif(g[[column]], 3)
    ytab_datatable(g[order(g$p_value), , drop = FALSE], page_length = 15L)
  })
  gi_name <- function(ext) paste0(active()$project_id, "_gi_", input$gi_plot %||% "scatter", "_", basename(gi_selected()), ext)
  for (format in c("png", "pdf")) local({ ext <- format
    output[[paste0("gi_download_", ext)]] <- downloadHandler(filename = function() gi_name(paste0(".", ext)), content = function(file) ytab_save_ggplot(file, gi_current_plot(), width = 10, height = 8)) })
  output$gi_download_genes <- downloadHandler(filename = function() gi_name("_genes.tsv"), content = function(file) file.copy(file.path(gi_selected(), "genes.tsv"), file))
  invisible(list(norm_runs = norm_runs, mc_runs = mc_runs, gi_runs = gi_runs, mc_data = mc_data, gi_data = gi_data, mc_assignment = mc_assignment, gi_assignment = gi_assignment))
}
