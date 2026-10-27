# Gene Explorer > Across Species > Essentiality across species: label-free HMM calls of one-to-one orthologs in two or
# more projects (`ytab cross-species`; S. cerevisiae systematic names as the hub). The command writes every table; the
# tab draws the shared and species-specific essential sets, the agreement between species and the ortholog table.

cross_species_ui <- function() ytab_two_column_layout(
  controls = ytab_control_panel("Essentiality across species",
    uiOutput("xs_project_selector"),
    helpText("One project per species, each with an essentiality run (its newest run is used). Genes are matched as one-to-one orthologs through S. cerevisiae."),
    actionButton("xs_run", "Compare species", class = "btn-primary"),
    tags$hr(),
    uiOutput("xs_result_selector"),
    radioButtons("xs_plot", "Plot", c("Shared and species-specific essential genes" = "sets", "Agreement between species" = "agreement")),
    ytab_plot_downloads_ui("xs", c(orthologs = "Download ortholog calls (TSV)", sets = "Download gene sets (TSV)", pairwise = "Download pairwise agreement (TSV)"))),
  main = tagList(
    uiOutput("xs_status"),
    ytab_plot_card("Essentiality across species", uiOutput("xs_plot_frame"), description = uiOutput("xs_plot_caption")),
    ytab_table_card("One-to-one orthologs", DT::DTOutput("xs_table"), description = "Search a gene (any species' ID or name, or the S. cerevisiae ID) to see its essentiality score and call in every species."),
    ytab_methods_drawer(
      tags$p(tags$b("Orthologs. "), "Each project's genes are mapped to S. cerevisiae systematic names with the ortholog table named in its project configuration (reference: orthology_file); an S. cerevisiae project is the hub itself. Only one-to-one orthologs are kept, and only genes scored by the HMM in every chosen species, so all species are compared on the same genes."),
      tags$p(tags$b("Calls. "), "Each species' label-free HMM score (P_ES_run, the posterior probability of a run of essential windows in the coding sequence) and call (essential when P_ES_run ≥ 0.5) from the project's newest essentiality run."),
      tags$p(tags$b("Agreement. "), "For each pair of species: Cohen's κ of the essential calls (agreement beyond chance; 1 = identical calls, 0 = chance) and the Spearman correlation of the scores. Differences combine biology (species-specific essentiality) with differences in library density and transposon; compare with the libraries' QC before reading a species-specific set as biology."),
      ytab_reference("Cohen J (1960). A coefficient of agreement for nominal scales. Educ Psychol Meas 20:37–46."),
      tags$p("Command: ", tags$code("ytab cross-species"), ". Results are written to ", tags$code("output/cross_species"), " and shown in every project they include."))))

# Species names in italics: HTML for the interactive plot, plotmath for downloads.
xs_species_labels <- function(species, html = TRUE) if (html) paste0("<i>", species, "</i>") else parse(text = sprintf('bolditalic("%s")', species))

# Bars: genes per set of species in which they are essential (most common first); below, the species of each set.
# Species-specific sets take the species' colour, genes essential in all species teal, others grey.
xs_sets_ggplot <- function(inter, species, html = TRUE, top = 15L) {
  d <- utils::head(inter[order(-inter$genes), , drop = FALSE], top); d$set <- factor(seq_len(nrow(d)))
  d$group <- ifelse(d$n_species == 0, "none", ifelse(d$n_species == length(species), "all", ifelse(d$n_species == 1, d$essential_in, "some")))
  fills <- c(all = YTAB_FOCAL, some = YTAB_COMPARATOR, none = "grey90", stats::setNames(ytab_species_colour(species), species))
  ymax <- max(d$genes); step <- ymax * .09; rows <- -step * seq_along(species)
  dots <- do.call(rbind, lapply(seq_along(species), function(i) data.frame(set = d$set, y = rows[[i]], on = d[[species[[i]]]] == 1, species = species[[i]])))
  links <- do.call(rbind, lapply(seq_len(nrow(d)), function(k) { on <- rows[d[k, species] == 1]
    if (length(on) > 1) data.frame(set = d$set[[k]], ymin = min(on), ymax = max(on)) }))
  d$tip <- sprintf("Essential in: %s<br>%s genes", d$essential_in, format(d$genes, big.mark = ","))
  ticks <- pretty(c(0, ymax)); ticks <- ticks[ticks <= ymax * 1.05]
  labels <- if (html) c(format(ticks, big.mark = ","), xs_species_labels(species)) else c(parse(text = sprintf('bold("%s")', format(ticks, big.mark = ","))), xs_species_labels(species, FALSE))
  p <- ggplot2::ggplot(d, ggplot2::aes(set, genes)) + ggplot2::geom_col(ggplot2::aes(fill = group, text = tip), width = .7) +
    # Placed above the bar by position, not vjust: plotly ignores vjust and would centre the count on the bar top.
    ggplot2::geom_text(ggplot2::aes(y = genes + ymax * .045, label = format(genes, big.mark = ",")), size = 3.4, fontface = "bold") +
    ggplot2::geom_point(data = dots, ggplot2::aes(set, y, colour = on), size = 3.2, inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey85"), guide = "none") +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::scale_y_continuous(breaks = c(ticks, rows), labels = labels, expand = ggplot2::expansion(mult = c(.02, .1))) +
    ggplot2::scale_x_discrete(labels = NULL) +
    ggplot2::labs(x = "Set of species in which the genes are essential", y = "One-to-one orthologs") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  if (!is.null(links) && nrow(links)) p <- p + ggplot2::geom_segment(data = links, ggplot2::aes(x = set, xend = set, y = ymin, yend = ymax), linewidth = .8, inherit.aes = FALSE)
  p
}

# Pairwise agreement: Cohen's kappa of the calls (fill) with the Spearman correlation of the scores, both triangles.
xs_agreement_ggplot <- function(pairwise, species, html = TRUE) {
  both <- rbind(pairwise, transform(pairwise, species_a = species_b, species_b = pairwise$species_a))
  both$species_a <- factor(both$species_a, species); both$species_b <- factor(both$species_b, rev(species))
  both$label <- sprintf("κ %.2f\nρ %.2f", both$kappa_ES, both$spearman_P_ES_run)
  both$tip <- sprintf("%s vs %s<br>%s genes; essential in both %s, in either %s<br>Cohen's κ %.3f; Spearman ρ %.3f", both$species_a, both$species_b,
                      format(both$genes, big.mark = ","), format(both$both_ES, big.mark = ","), format(both$either_ES, big.mark = ","), both$kappa_ES, both$spearman_P_ES_run)
  ggplot2::ggplot(both, ggplot2::aes(species_a, species_b, fill = kappa_ES, text = tip)) + ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = label), fontface = "bold", size = 3.8) +
    ggplot2::scale_fill_gradient(low = "white", high = YTAB_FOCAL, limits = c(0, 1), name = "Cohen's κ") +
    ggplot2::scale_x_discrete(labels = xs_species_labels(species, html)) + ggplot2::scale_y_discrete(labels = xs_species_labels(rev(species), html)) +
    ggplot2::coord_equal() + ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank())
}

cross_species_server <- function(input, output, session, repo_root, active, active_project_path, jobs, python_bin) {
  # Comparisons span projects, so they are kept in one place (output/cross_species) and every project they include shows them.
  runs_dir <- file.path(repo_root, "output", "cross_species")
  all_runs <- reactivePoll(4000, session,
    checkFunc = function() { x <- ytab_analysis_runs(runs_dir); paste(x$path, x$modified, collapse = "|") },
    valueFunc = function() ytab_analysis_runs(runs_dir))
  runs <- reactive({
    x <- all_runs(); if (!nrow(x)) return(x)
    keep <- vapply(x$path, function(run) active()$project_id %in% vapply((ytab_read_manifest(run) %||% list())$projects %||% list(), function(p) p$project %||% "", ""), TRUE)
    x[keep, , drop = FALSE]
  })
  # projects with an essentiality run, newest run first
  candidates <- reactive({
    x <- discover_ytab_projects(repo_root, include_tests = FALSE); if (!nrow(x)) return(x)
    x$essentiality <- vapply(dirname(dirname(x$project_config)), function(root) { r <- ytab_analysis_runs(file.path(root, "essentiality"))
      if (nrow(r)) as.numeric(r$modified[[1L]]) else NA_real_ }, 0)
    x[!is.na(x$essentiality), , drop = FALSE]
  })
  output$xs_project_selector <- renderUI({
    x <- candidates()
    if (nrow(x) < 2L) return(tags$p(class = "text-muted", "Needs at least two projects of different species with an essentiality run (Essentiality → Run)."))
    # the active project, then for every other species the project with the newest essentiality run
    x <- x[order(x$project_id != active()$project_id, -x$essentiality), , drop = FALSE]
    chosen <- x$project_config[!duplicated(x$species)]
    checkboxGroupInput("xs_projects", "Projects", stats::setNames(x$project_config, sprintf("%s (%s)", x$display_name, x$species)), selected = chosen)
  })
  observeEvent(input$xs_run, {
    if (jobs$job_is_running()) { showNotification("Wait for the current job to finish.", type = "warning"); return() }
    chosen <- input$xs_projects %||% character(); species <- candidates()$species[match(chosen, candidates()$project_config)]
    if (length(chosen) < 2L) { showNotification("Choose at least two projects.", type = "error"); return() }
    if (anyDuplicated(species)) { showNotification("Choose one project per species.", type = "error"); return() }
    args <- c(file.path(repo_root, "scripts", "local", "ytab_cross_species.py"), as.vector(rbind("--project-config", chosen)),
              "--output", ytab_new_run_dir(runs_dir, "cross_species"))
    tryCatch({
      jobs$start_job(python_bin(), args, active_project_path(), "cross_species", wd = repo_root)
      showNotification("Comparing species. The result appears here when it finishes.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error"))
  }, ignoreInit = TRUE)
  output$xs_result_selector <- renderUI({
    x <- runs(); if (!nrow(x)) return(NULL)
    selectInput("xs_result", "Result", ytab_run_choices(x))
  })
  selected_run <- reactive({ req(nrow(runs()) > 0); path <- input$xs_result %||% runs()$path[[1L]]; req(path %in% runs()$path); path })
  data <- reactive({
    run <- selected_run(); m <- ytab_read_manifest(run)
    list(manifest = m, species = unlist(m$species), orthologs = ytab_read_analysis_tsv(file.path(run, "orthologs.tsv")),
         intersections = ytab_read_analysis_tsv(file.path(run, "intersections.tsv")), pairwise = ytab_read_analysis_tsv(file.path(run, "pairwise.tsv")))
  })
  output$xs_status <- renderUI({
    if (!nrow(runs())) return(ytab_empty_state("No comparison yet", "Choose one project per species and compare. Each project needs an essentiality run."))
    m <- data()$manifest
    tags$p(class = "text-muted", sprintf("%s one-to-one orthologs scored in %s; essential in all: %s; in none: %s.", format(m$orthologs, big.mark = ","),
                                         paste(unlist(m$species), collapse = ", "), format(m$essential_in_all, big.mark = ","), format(m$essential_in_none, big.mark = ",")))
  })
  plot_for <- function(html) {
    d <- data(); validate(need(nrow(d$intersections) > 0, "No result tables in this run."))
    if (identical(input$xs_plot, "agreement")) xs_agreement_ggplot(d$pairwise, d$species, html) else xs_sets_ggplot(d$intersections, d$species, html)
  }
  output$xs_plot_caption <- renderUI(if (identical(input$xs_plot, "agreement"))
    "Agreement of the essential calls between each pair of species (Cohen's κ, fill) and the correlation of their HMM scores (Spearman ρ), on the same one-to-one orthologs."
    else "Number of one-to-one orthologs essential in each set of species (dots below each bar; the 15 largest sets): bars in a species' colour are essential only in that species, teal in all species.")
  output$xs_plot_frame <- renderUI({ req(nrow(runs()) > 0); plotly::plotlyOutput("xs_plotly", height = "560px") })
  output$xs_plotly <- plotly::renderPlotly(ytab_ggplotly(plot_for(TRUE)))
  output$xs_table <- DT::renderDT({
    o <- data()$orthologs; req(nrow(o) > 0)
    for (s in data()$species) { o[[paste(s, "ES")]] <- ifelse(o[[paste(s, "ES")]] == 1, "essential", "")
      o[[paste(s, "name")]] <- ytab_display_gene_names(o[[paste(s, "gene")]], o[[paste(s, "name")]], repo_root) }
    names(o)[names(o) == "scer_id"] <- "S. cerevisiae ID"
    ytab_datatable(o, page_length = 10L)
  })
  name <- function(suffix) paste0(active()$project_id, "_cross_species_", basename(selected_run()), suffix)
  for (format in c("png", "pdf")) local({
    ext <- format
    output[[paste0("xs_download_", ext)]] <- downloadHandler(filename = function() name(paste0("_", input$xs_plot %||% "sets", ".", ext)),
      content = function(file) ytab_save_ggplot(file, plot_for(FALSE), width = if (identical(input$xs_plot, "agreement")) 7 else 10, height = if (identical(input$xs_plot, "agreement")) 6 else 6.5))
  })
  tables <- c(orthologs = "orthologs.tsv", sets = "intersections.tsv", pairwise = "pairwise.tsv")
  for (table in names(tables)) local({
    key <- table; file_name <- tables[[key]]
    output[[paste0("xs_download_", key)]] <- downloadHandler(filename = function() name(paste0("_", file_name)),
      content = function(file) file.copy(file.path(selected_run(), file_name), file))
  })
  invisible(list(runs = runs, selected_run = selected_run, data = data, candidates = candidates))
}
