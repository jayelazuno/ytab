# Treated vs Control: three ways to call genes depleted or enriched in treated libraries relative to their controls.
# The local-abundance z-score (default) uses the Fitness run and results views (ytab fitness); resampling and the
# control-conditioned Gumbel run on a Site Normalization run (ytab treated-vs-control) and have their own views here.

tvc_method_choices <- c("Local-abundance z-score" = "zscore",
                        "Resampling" = "resampling",
                        "Control-conditioned Gumbel (depletion only)" = "gumbel")   # validation status: Methods drawer

# Normalizations each method accepts: the Gumbel needs a scalar factor (its count model uses it as the exposure);
# resampling takes the scalar methods and the beta-geometric corrections.
tvc_allowed_methods <- list(resampling = c("betageom", "aBGC", "ttr", "totreads", "nzmean"), gumbel = c("ttr", "totreads", "nzmean"))

# Normalization choices grouped by Site Normalization run, latest first; values are "<run>|<method>". `allowed`: the
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
  tags$p(tags$b("Local-abundance z-score (default). "), "The replicate libraries of each condition are averaged; each gene's log2 fold change between the treated and control means is divided by the noise expected at the gene's abundance. The noise model is fitted to the two controls of each background (running SD of control–control log2 ratios against abundance, SD = m1 + m2·A^m3), and the call threshold is the 99th percentile of control–control |z|. Default input: BGC from the latest Site Normalization run (gene totals of its normalized site counts, expressed per million of the mean library total so the 0.5 pseudocount means the same as with counts per million); total reads (counts per million) and any other Site Normalization method can be chosen. Calls: depleted, enriched or none (one call per gene from the replicate means); the per-pool z-scores are kept as a diagnostic."),
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

tvc_results_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Results",
    uiOutput("tvc_result_selector"),
    uiOutput("tvc_group_selector"),
    uiOutput("tvc_plot_selector"),
    uiOutput("tvc_gene_selector"),
    ytab_plot_downloads_ui("tvc", c(genes = "Download genes (TSV)"))),
  main = tagList(uiOutput("tvc_summary"),
    ytab_plot_card("Treated vs control", uiOutput("tvc_plot_frame"), description = uiOutput("tvc_plot_caption")),
    ytab_table_card("Genes", DT::DTOutput("tvc_table")), tvc_methods_text()))

tvc_run_tab_ui <- function() tagList(
  radioButtons("tvc_method", "Method", tvc_method_choices, "zscore"),
  conditionalPanel("input.tvc_method == 'zscore'", fitness_run_ui()),
  conditionalPanel("input.tvc_method != 'zscore'", tvc_run_ui()))

tvc_results_tab_ui <- function() tagList(
  radioButtons("tvc_results_method", "Method", tvc_method_choices, "zscore"),
  conditionalPanel("input.tvc_results_method == 'zscore'", fitness_results_ui()),
  conditionalPanel("input.tvc_results_method != 'zscore'", tvc_results_ui()))

# Gene-table x axis for the volcano: resampling has log2FC, the Gumbel the scaled read ratio of the gene.
tvc_effect_column <- function(genes) if ("log2FC" %in% names(genes)) "log2FC" else "log2_scaled_ratio"

tvc_plot_choices <- function(method) if (identical(method, "gumbel"))
  c("Volcano" = "volcano", "Observed vs expected run" = "runs", "p-value histogram" = "pvalues", "Gene sites" = "sites") else
  c("Volcano" = "volcano", "MA plot" = "ma", "p-value histogram" = "pvalues", "Gene sites" = "sites")

tvc_tested <- function(genes) genes[genes$status %in% "tested" | genes$in_family %in% c("True", TRUE), , drop = FALSE]

tvc_call_colours <- c(depleted = "#b2182b", `depleted run` = "#b2182b", enriched = "#2166ac", `not significant` = "grey70")

# One plot of a resampling or Gumbel result group. sites: the group's sites.tsv.gz rows for `gene` (gene sites only).
tvc_plot <- function(choice, genes, method, sites = NULL, gene = NULL) {
  base <- ggplot2::theme_minimal(base_size = 12)
  levels <- if (identical(method, "gumbel")) c("depleted run", "not significant") else c("depleted", "enriched", "not significant")
  d <- tvc_tested(genes); d$call <- factor(d$call, levels = levels)
  hover <- function(d) sprintf("%s %s<br>p %.3g; q %.3g", d$ORF, ifelse(is.na(d$Name), "", d$Name), d$p_value, d$q_bh)
  if (identical(choice, "ma")) {
    d <- d[is.finite(d$log2FC) & is.finite(d$mean_abundance) & d$mean_abundance > 0, , drop = FALSE]
    d$A <- log2(d$mean_abundance); d$hover <- hover(d)
    return(ggplot2::ggplot(d, ggplot2::aes(A, log2FC, colour = call, text = hover)) + ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
      ggplot2::geom_point(alpha = .6, size = 1.3) + ggplot2::scale_colour_manual(values = tvc_call_colours, drop = FALSE) +
      ggplot2::labs(x = "log2 mean normalized count per site (abundance)", y = "log2FC (treated / control)", colour = "Call") + base)
  }
  if (identical(choice, "pvalues")) {
    d <- d[is.finite(d$p_value), , drop = FALSE]
    return(ggplot2::ggplot(d, ggplot2::aes(p_value)) + ggplot2::geom_histogram(breaks = seq(0, 1, .05), fill = "#7fa7c9") +
      ggplot2::geom_hline(yintercept = nrow(d) / 20, linetype = "dashed", colour = "grey30") +
      ggplot2::labs(x = "p-value", y = "Genes") + base)
  }
  if (identical(choice, "runs")) {
    d <- d[is.finite(d$r) & is.finite(d$null_mean_run), , drop = FALSE]
    d$hover <- sprintf("%s<br>observed run %s sites; expected %.1f", hover(d), d$r, d$null_mean_run)
    return(ggplot2::ggplot(d, ggplot2::aes(null_mean_run, r, colour = call, text = hover)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40") + ggplot2::geom_jitter(width = .05, height = .15, alpha = .6, size = 1.3) +
      ggplot2::scale_colour_manual(values = tvc_call_colours, drop = FALSE) +
      ggplot2::labs(x = "Expected longest depleted run (sites; control model)", y = "Observed longest depleted run (sites)", colour = "Call") + base)
  }
  if (identical(choice, "sites")) {
    if (is.null(sites) || !nrow(sites)) return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "No sites for this gene") + ggplot2::theme_void())
    libs <- grep("^(control|treatment)_[0-9]+$", names(sites), value = TRUE)
    long <- do.call(rbind, lapply(libs, function(l) data.frame(position = sites$position, library = l, count = as.numeric(sites[[l]]),
      condition = ifelse(startsWith(l, "control"), "control", "treated"), depleted = if ("depleted" %in% names(sites)) sites$depleted %in% c("True", TRUE) else FALSE)))
    long$library <- sub("^treatment_", "treated ", sub("^control_", "control ", long$library))
    shades <- c(stats::setNames(c("grey25", "grey60", "grey80")[seq_len(sum(startsWith(unique(long$library), "control")))], sort(unique(long$library[long$condition == "control"]))),
                stats::setNames(c("#b2182b", "#ef8a62", "#fddbc7")[seq_len(sum(startsWith(unique(long$library), "treated")))], sort(unique(long$library[long$condition == "treated"]))))
    long$hover <- sprintf("%s<br>position %s<br>%s %.1f", gene, format(long$position, big.mark = ","), long$library, long$count)
    p <- ggplot2::ggplot(long, ggplot2::aes(position, count, colour = library, text = hover)) +
      ggplot2::geom_point(size = 2, alpha = .8) + ggplot2::scale_colour_manual(values = shades) +
      ggplot2::labs(x = sprintf("%s position (bp)", gene), y = if (identical(method, "gumbel")) "Reads per site" else "Normalized count per site",
                    colour = "Library") + base
    if (any(long$depleted)) p <- p + ggplot2::geom_rug(data = long[long$depleted & long$condition == "control", , drop = FALSE], sides = "b", colour = "#b2182b", inherit.aes = FALSE, ggplot2::aes(x = position))
    return(p)
  }
  tvc_ggplot(genes, method)
}

tvc_ggplot <- function(genes, method) {
  x <- tvc_effect_column(genes)
  d <- genes[genes$status %in% "tested" | genes$in_family %in% c("True", TRUE), , drop = FALSE]
  d <- d[is.finite(d[[x]]) & is.finite(d$p_value), , drop = FALSE]
  levels <- if (identical(method, "gumbel")) c("depleted run", "not significant") else c("depleted", "enriched", "not significant")
  d$call <- factor(d$call, levels = levels)
  ggplot2::ggplot(d, ggplot2::aes(.data[[x]], -log10(pmax(p_value, 1e-300)), colour = call,
      text = sprintf("%s %s<br>%s %.2f<br>p %.3g; q %.3g", ORF, ifelse(is.na(Name), "", Name), x, .data[[x]], p_value, q_bh))) +
    ggplot2::geom_point(alpha = .6, size = 1.3) +
    ggplot2::scale_colour_manual(values = c(depleted = "#b2182b", `depleted run` = "#b2182b", enriched = "#2166ac", `not significant` = "grey70"), drop = FALSE) +
    ggplot2::labs(x = if (identical(x, "log2FC")) "log2FC (treated / control)" else "log2 scaled read ratio (treated / control)", y = "-log10 p", colour = "Call") +
    ggplot2::theme_minimal(base_size = 12)
}

tvc_server <- function(input, output, session, repo_root, active, active_project_path, norm_runs, start) {
  tvc_runs <- reactivePoll(3000, session,
    checkFunc = function() { x <- ytab_analysis_runs(tryCatch(file.path(active()$project_root, "treated_vs_control"), error = function(e) "")); paste(x$path, x$modified, collapse = "|") },
    valueFunc = function() ytab_analysis_runs(tryCatch(file.path(active()$project_root, "treated_vs_control"), error = function(e) "")))
  design <- reactive({
    path <- file.path(active()$project_root, "config", "comparison_design.csv")
    if (!file.exists(path)) return(data.frame())
    d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    d[tolower(trimws(as.character(d$include))) %in% c("1", "true", "yes", "y"), , drop = FALSE]
  })

  # Z-score input: total reads (CPM) or any method of a Site Normalization run ("<run>|<method>").
  output$fitness_normalization_selector <- renderUI({
    # Default: BGC of the latest run that has it (total reads only when no run has BGC); a valid earlier choice is kept.
    choices <- c(list("Without a Site Normalization run" = c("Total reads (CPM)" = "cpm")), tvc_norm_choices(norm_runs()))
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
    if (!nrow(runs)) return(tags$p(class = "text-warning", "No normalization run yet. Run Normalization → Site Normalization first."))
    choices <- tvc_norm_choices(runs, allowed)
    if (!length(choices)) return(tags$p(class = "text-warning", sprintf("No Site Normalization run has a normalization the %s accepts (%s). Add one in Site Normalization.",
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
  method_runs <- reactive({
    x <- tvc_runs(); method <- input$tvc_results_method %||% "resampling"
    if (!nrow(x)) return(x)
    keep <- vapply(x$path, function(p) identical(ytab_read_manifest(p)$method, method), logical(1))
    x[keep, , drop = FALSE]
  })
  output$tvc_result_selector <- renderUI({
    x <- method_runs()
    if (!nrow(x)) return(tags$p(class = "text-muted", "No result for this method yet; run it from Treated vs Control: Run."))
    selectInput("tvc_result", "Result", ytab_run_choices(x))
  })
  selected <- reactive({ x <- method_runs(); req(nrow(x) > 0); path <- input$tvc_result %||% x$path[[1L]]; req(path %in% x$path); path })
  data <- reactive({ run <- selected(); list(manifest = ytab_read_manifest(run), genes = ytab_read_analysis_tsv(file.path(run, "genes.tsv"))) })
  output$tvc_group_selector <- renderUI({
    g <- data()$genes; req(nrow(g) > 0); groups <- unique(g$group)
    if (length(groups) < 2L) return(NULL)
    m <- data()$manifest$groups
    selectInput("tvc_group", "Background and condition", stats::setNames(groups, vapply(groups, function(k) paste(m[[k]]$background, m[[k]]$condition, sep = ", "), "")))
  })
  group_genes <- reactive({ g <- data()$genes; req(nrow(g) > 0); group <- input$tvc_group %||% g$group[[1L]]; g[g$group == group, , drop = FALSE] })
  output$tvc_summary <- renderUI({
    if (!nrow(method_runs())) return(ytab_empty_state("No result yet", "Choose a normalization run and comparisons in Treated vs Control: Run."))
    m <- data()$manifest; req(m)
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
  output$tvc_plot_selector <- renderUI({ m <- data()$manifest; req(m); radioButtons("tvc_plot", "Plot", tvc_plot_choices(m$method), "volcano") })
  output$tvc_gene_selector <- renderUI({
    req(identical(input$tvc_plot, "sites")); g <- tvc_tested(group_genes()); req(nrow(g) > 0)
    g <- g[order(g$q_bh), , drop = FALSE]
    selectizeInput("tvc_gene", "Gene", stats::setNames(g$ORF, ifelse(is.na(g$Name) | g$Name == g$ORF, g$ORF, paste0(g$ORF, " (", g$Name, ")"))), selected = g$ORF[[1L]])
  })
  group_sites <- reactive({
    req(identical(input$tvc_plot, "sites")); group <- input$tvc_group %||% group_genes()$group[[1L]]
    path <- file.path(selected(), group, "sites.tsv.gz")
    if (!file.exists(path)) return(data.frame())
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  })
  gene_sites <- reactive({ s <- group_sites(); gene <- input$tvc_gene %||% ""; if (!nrow(s) || !nzchar(gene)) return(s[0, , drop = FALSE])
                           s[vapply(strsplit(as.character(s$genes), "[;,]"), function(x) gene %in% x, logical(1)), , drop = FALSE] })
  output$tvc_plot_caption <- renderUI(switch(input$tvc_plot %||% "volcano",
    ma = "Each point is a gene: mean normalized count per site (abundance) against log2 fold change; colours are calls at the chosen FDR.",
    pvalues = "p-values of the tested genes. With no effect they are uniform (dashed line: expected count per bin); a peak near 0 holds the affected genes.",
    runs = "Each gene's longest run of depleted sites against the longest run expected by chance from the controls' count model; genes far above the line are depleted runs.",
    sites = if (identical(data()$manifest$method, "gumbel")) "Reads at each insertion site of the gene, per library; red ticks mark sites the test counts as depleted." else
      "Normalized count at each insertion site of the gene, per library (before recentering).",
    if (identical(data()$manifest$method, "gumbel")) "Each point is a gene tested for a depleted run of sites; the x axis is the gene's treated/control read ratio after scaling." else
      "Each point is a gene; colours are calls at the chosen FDR."))
  output$tvc_plot_frame <- renderUI({ req(nrow(group_genes()) > 0); plotly::plotlyOutput("tvc_plotly", height = "520px") })
  current_plot <- function() { choice <- input$tvc_plot %||% "volcano"
    tvc_plot(choice, group_genes(), data()$manifest$method, if (identical(choice, "sites")) gene_sites() else NULL, input$tvc_gene) }
  output$tvc_plotly <- plotly::renderPlotly(ytab_ggplotly(current_plot()))
  output$tvc_table <- DT::renderDT({
    g <- group_genes(); keep <- intersect(c("ORF", "Name", "n_sites", "n_eligible", "n_depleted", "run_span_bp", "mean_control", "mean_treatment",
                                            "log2FC", "log2_scaled_ratio", "p_value", "q_bh", "call"), names(g))
    g <- g[, keep, drop = FALSE]
    for (column in intersect(c("mean_control", "mean_treatment", "log2FC", "log2_scaled_ratio", "p_value", "q_bh"), names(g))) g[[column]] <- signif(g[[column]], 3)
    pv <- intersect(c("p_value", "q_bh"), names(g))
    for (column in pv) g[[column]] <- formatC(g[[column]], format = "g", digits = 3)
    ytab_datatable(g[order(as.numeric(g$q_bh)), , drop = FALSE], page_length = 15L)
  })
  name <- function(ext) paste0(active()$project_id, "_treated_vs_control_", basename(selected()), ext)
  for (format in c("png", "pdf")) local({ ext <- format
    output[[paste0("tvc_download_", ext)]] <- downloadHandler(filename = function() name(paste0(".", ext)), content = function(file) ytab_save_ggplot(file, current_plot(), width = 10, height = 8)) })
  output$tvc_download_genes <- downloadHandler(filename = function() name("_genes.tsv"), content = function(file) file.copy(file.path(selected(), "genes.tsv"), file))
}
