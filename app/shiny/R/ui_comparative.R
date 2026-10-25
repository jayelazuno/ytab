comparative_metric_choices <- c(
  "Essentiality prediction" = "essentiality",
  "Fitness score/call" = "fitness",
  "Insertion summary" = "insertion"
)

comparative_controls_card <- function(title, ...) {
  tags$aside(class = "ytab-comparative-controls ytab-control-panel", panel_card(title, ...))
}

comparative_result_card <- function(title, ...) {
  tags$section(class = "ytab-comparative-results ytab-release-card", panel_card(title, ...))
}

comparative_layout <- function(controls, results) {
  tags$div(
    class = "row ytab-comparative-layout",
    tags$div(class = "col-sm-4 col-md-3", controls),
    tags$div(class = "col-sm-8 col-md-9", results)
  )
}

comparative_single_species_ui <- function() {
  comparative_layout(
    comparative_controls_card(
      "Single species controls",
      selectInput("comparative_single_species", "Species", choices = character()),
      selectInput("comparative_single_project", "Project", choices = character()),
      textInput("comparative_single_gene", "Gene name or ID", ""),
      selectInput("comparative_single_metric", "Metric",
                  choices = comparative_metric_choices),
      actionButton("comparative_single_search", "Search", class = "btn-primary")
    ),
    comparative_result_card(
      "Single Species View",
      tags$p(class = "ytab-stage-purpose",
             "Search for one gene to view available Tn-seq results for one species."),
      uiOutput("comparative_single_state"),
      uiOutput("comparative_single_gene_card"),
      DT::DTOutput("comparative_single_identifiers"),
      DT::DTOutput("comparative_single_results")
    )
  )
}

comparative_species_checkbox_ui <- function() {
  tagList(
    checkboxGroupInput("comparative_species_set", "Species",
                       choices = character(), selected = character())
  )
}

comparative_comparison_ui <- function() {
  comparative_layout(
    comparative_controls_card(
      "Comparative controls",
      textInput("comparative_gene", "Gene name or ID", ""),
      uiOutput("comparative_species_checks"),
      selectInput("comparative_metric", "Metric",
                  choices = c("Essentiality prediction" = "essentiality",
                              "Fitness call" = "fitness_call",
                              "Fitness score" = "fitness_score",
                              "Insertion summary" = "insertion")),
      actionButton("comparative_compare", "Compare", class = "btn-primary")
    ),
    comparative_result_card(
      "Comparative View",
      tags$p(class = "ytab-stage-purpose",
             "Compare orthology mapping and available YTAB outputs across supported species."),
      uiOutput("comparative_state"),
      uiOutput("comparative_orthogroup_summary"),
      DT::DTOutput("comparative_ortholog_table"),
      DT::DTOutput("comparative_availability_table"),
      DT::DTOutput("comparative_result_table")
    )
  )
}

comparative_gene_group_ui <- function() {
  comparative_layout(
    comparative_controls_card(
      "Gene group controls",
      fileInput("comparative_group_file", "Upload CSV gene group",
                accept = c(".csv", "text/csv")),
      textAreaInput("comparative_group_text", "Paste gene list", "",
                    rows = 6, placeholder = "One gene ID or name per line"),
      selectInput("comparative_group_species", "Species", choices = character()),
      checkboxInput("comparative_group_cross_species",
                    "Include cross-species ortholog analysis", TRUE),
      selectInput("comparative_group_metric", "Metric",
                  choices = comparative_metric_choices),
      actionButton("comparative_group_analyze", "Analyze group",
                   class = "btn-primary")
    ),
    comparative_result_card(
      "Gene Group Analysis",
      tags$p(class = "ytab-stage-purpose",
             "Map a gene list to orthogroups and summarize available YTAB outputs."),
      uiOutput("comparative_group_state"),
      uiOutput("comparative_group_summary"),
      DT::DTOutput("comparative_group_unmapped"),
      DT::DTOutput("comparative_group_availability"),
      DT::DTOutput("comparative_group_results")
    )
  )
}

comparative_ui <- function() {
  navset_tab(
    id = "comparative_tabs",
    nav_panel("Single Species View", value = "single_species",
              comparative_single_species_ui()),
    nav_panel("Comparative View", value = "comparative",
              comparative_comparison_ui()),
    nav_panel("Gene Group Analysis", value = "gene_group",
              comparative_gene_group_ui()),
    nav_panel("Essentiality across species", value = "essentiality_species",
              cross_species_ui())
  )
}

# Rows whose names or IDs contain the query (any case), or whose ID is one of `ids` (the gene's orthologs, matched
# exactly; a GFF prefix such as "gene:" is ignored).
comparative_row_matches <- function(data, query, ids = character()) {
  if (!nrow(data)) return(rep(FALSE, 0L))
  query <- tolower(trimws(as.character(query %||% "")))
  if (!nzchar(query) && !length(ids)) return(rep(FALSE, nrow(data)))
  fields <- intersect(c(
    "Standard name", "Common name", "Sc ortholog", "Sc std name",
    "feature_id", "standard_name", "common_name", "gene_id", "gene",
    "Feature ID", "Gene", "ORF", "Name"
  ), names(data))
  if (!length(fields)) return(rep(FALSE, nrow(data)))
  bare <- function(x) sub("^[a-z]+:", "", tolower(as.character(x)))
  ids <- bare(ids[nzchar(ids)])
  Reduce(`|`, lapply(data[, fields, drop = FALSE], function(value)
    (nzchar(query) & grepl(query, tolower(as.character(value)), fixed = TRUE)) | bare(value) %in% ids))
}

comparative_display_rows <- function(data, query = "", ids = character()) {
  if (!nrow(data)) return(data.frame())
  rows <- if (nzchar(trimws(query)) || length(ids)) data[comparative_row_matches(data, query, ids), ,
                                         drop = FALSE] else data
  if (!nrow(rows)) return(data.frame())
  keep <- intersect(c(
    "ytab_species", "ytab_project_id", "ytab_result_type", "ytab_target_tag",
    "ytab_analysis_id", "Standard name", "Common name", "Sc ortholog",
    "RF - G4 - ess. for FPR 0.100", "ORF", "Name", "hmm_call", "hmm_score", "hybrid_score", "agreement", "feature_id", "standard_name",
    "common_name", "final_call", "mean_log2FC", "max_abs_z", "sample",
    "total_reads", "total_hits", "percent_features_hit"
  ), names(rows))
  if (!length(keep)) keep <- names(rows)[seq_len(min(8L, ncol(rows)))]
  rows[seq_len(min(100L, nrow(rows))), keep, drop = FALSE]
}

comparative_load_metric <- function(project_config, metric) {
  if (!nzchar(project_config %||% "")) return(data.frame())
  if (metric %in% c("essentiality"))
    return(load_project_essentiality_results(project_config))
  if (metric %in% c("fitness", "fitness_call", "fitness_score"))
    return(load_project_fitness_results(project_config))
  load_project_summary_stats(project_config)
}

comparative_availability_display <- function(availability) {
  if (!is.data.frame(availability) || !nrow(availability)) return(data.frame())
  data.frame(
    Species = as.character(availability$label %||% availability$species),
    `YTAB project outputs` = as.character(availability$status %||% ""),
    Projects = as.integer(availability$project_count %||% 0L),
    `Essentiality results` = as.integer(availability$classifier_projects %||% 0L),
    `Fitness results` = as.integer(availability$fitness_projects %||% 0L),
    `Insertion summaries` = as.integer(availability$summary_projects %||% 0L),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

comparative_server <- function(input, output, session, repo_root, active = function() list()) {
  keep_or <- function(value, fallback) if (length(value) && nzchar(value[[1]])) value else fallback
  species_manifest <- reactive(load_comparative_species_manifest(repo_root))
  project_inventory <- reactive(list_ytab_projects_by_species(repo_root))
  gene_lookup <- reactive(load_comparative_gene_lookup(repo_root))
  # other IDs of the same gene (e.g. C. glabrata CAGL0 IDs in the orthogroups, GWK60 IDs in projects), from the annotation lookup
  id_aliases <- reactive({ x <- tryCatch(ytab_load_glabrata_annotation_lookup(repo_root), error = function(e) data.frame())
    if (all(c("cagl_id", "gwk60_id_clean") %in% names(x))) x[nzchar(x$cagl_id) & nzchar(x$gwk60_id_clean), c("cagl_id", "gwk60_id_clean")] else data.frame(cagl_id = character(), gwk60_id_clean = character()) })
  with_aliases <- function(ids) unique(c(ids, id_aliases()$gwk60_id_clean[id_aliases()$cagl_id %in% ids]))
  orthology <- reactive(load_comparative_orthogroups(repo_root))
  availability <- reactive(summarize_species_project_availability(repo_root))

  observe({
    species <- species_manifest()
    enabled <- species[species$enabled, , drop = FALSE]
    choices <- as.list(setNames(enabled$species, enabled$label))
    # default: the open project's species, else the first enabled species
    current <- tolower(as.character(tryCatch(active()$species, error = function(e) "") %||% ""))
    first <- if (current %in% enabled$species) current else if (length(choices)) unname(unlist(choices, use.names = FALSE))[[1]] else ""
    selected <- keep_or(isolate(input$comparative_single_species), first)
    updateSelectInput(session, "comparative_single_species",
                      choices = choices, selected = selected)
    updateSelectInput(session, "comparative_group_species",
                      choices = choices,
                      selected = keep_or(isolate(input$comparative_group_species), selected))
  })

  observe({
    projects <- project_inventory()
    species <- input$comparative_single_species %||% ""
    subset <- projects[projects$species == species, , drop = FALSE]
    choices <- if (nrow(subset))
      as.list(setNames(subset$project_config, subset$display_name)) else list()
    # default: the open project when it is of this species
    open_config <- tryCatch(normalizePath(as.character(active()$project_config %||% ""), mustWork = FALSE), error = function(e) "")
    configs <- vapply(choices, function(x) normalizePath(x, mustWork = FALSE), "")
    selected <- if (nzchar(open_config) && open_config %in% configs) choices[[match(open_config, configs)]] else if (length(choices)) choices[[1]] else character()
    updateSelectInput(session, "comparative_single_project", choices = choices, selected = selected)
  })

  # rendered with its choices: an empty group rendered after the observer's update would otherwise replace them
  output$comparative_species_checks <- renderUI({
    species <- species_manifest(); enabled <- species[species$enabled, , drop = FALSE]
    checkboxGroupInput("comparative_species_set", "Species", choices = stats::setNames(enabled$species, enabled$label),
                       selected = isolate(input$comparative_species_set) %||% enabled$species)
  })

  output$comparative_single_state <- renderUI({
    if ((input$comparative_single_search %||% 0L) == 0L)
      return(tags$p(class = "text-muted",
                    "Enter a gene name or ID, choose a metric, then click Search."))
    species <- input$comparative_single_species %||% ""
    projects <- project_inventory()
    if (!any(projects$species == species))
      return(tags$div(class = "alert alert-warning",
                      "No YTAB project outputs are available for this species yet."))
    NULL
  })

  single_result <- eventReactive(input$comparative_single_search, {
    query <- input$comparative_single_gene %||% ""
    project_config <- input$comparative_single_project %||% ""
    metric <- input$comparative_single_metric %||% "essentiality"
    data <- comparative_load_metric(project_config, metric)
    list(query = query, data = comparative_display_rows(data, query),
         raw = data, metric = metric)
  }, ignoreInit = TRUE)

  output$comparative_single_gene_card <- renderUI({
    result <- single_result()
    req(result)
    ids <- gene_lookup()[comparative_query_matches(
      gene_lookup(), result$query, c("gene_id", "gene_name", "expression_id")
    ), , drop = FALSE]
    tags$div(
      class = "ytab-result-card",
      tags$h4(if (nzchar(result$query)) result$query else "No gene query"),
      tags$p(if (nrow(ids)) "Orthology identifiers are available." else
        "No orthology lookup match is available for this query."),
      tags$p(sprintf("Matching YTAB feature rows: %d", nrow(result$data)))
    )
  })

  output$comparative_single_identifiers <- DT::renderDT({
    result <- single_result()
    req(result)
    ids <- gene_lookup()[comparative_query_matches(
      gene_lookup(), result$query, c("gene_id", "gene_name", "expression_id")
    ), , drop = FALSE]
    if (!nrow(ids)) return(NULL)
    DT::datatable(ids, rownames = FALSE, selection = "none",
                  options = list(pageLength = 8, scrollX = TRUE))
  })

  output$comparative_single_results <- DT::renderDT({
    result <- single_result()
    req(result)
    if (!nrow(result$data)) return(NULL)
    DT::datatable(result$data, rownames = FALSE, selection = "none",
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  comparative_result <- eventReactive(input$comparative_compare, {
    query <- input$comparative_gene %||% ""
    species <- input$comparative_species_set %||% character()
    mapped <- map_gene_across_species(query, orthology(), gene_lookup())
    projects <- project_inventory()
    rows <- lapply(species, function(key) {
      subset <- projects[projects$species == key, , drop = FALSE]
      if (!nrow(subset)) return(data.frame())
      data <- comparative_load_metric(subset$project_config[[1]],
                                      input$comparative_metric %||% "essentiality")
      # the gene itself (by name or ID) and its orthologs in this species (orthogroup members)
      comparative_display_rows(data, query, with_aliases(as.character(mapped$gene_id[mapped$species == key])))
    })
    rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
    list(query = query, species = species, mapped = mapped,
         results = if (length(rows)) do.call(rbind, rows) else data.frame())
  }, ignoreInit = TRUE)

  output$comparative_state <- renderUI({
    if ((input$comparative_compare %||% 0L) == 0L)
      return(tags$p(class = "text-muted",
                    "Search for a gene to compare orthology and available YTAB outputs."))
    result <- comparative_result()
    available_species <- unique(project_inventory()$species)
    selected_available <- intersect(result$species, available_species)
    if (length(selected_available) <= 1L)
      tags$div(class = "alert alert-info",
               "Only one species currently has YTAB project results loaded.")
  })

  output$comparative_orthogroup_summary <- renderUI({
    result <- comparative_result()
    req(result)
    groups <- unique(result$mapped$orthogroup_id)
    tags$div(
      class = "ytab-result-card",
      tags$h4("Orthogroup match summary"),
      tags$p(sprintf("Mapped orthogroups: %d", length(groups))),
      tags$p(sprintf("Ortholog rows: %d", nrow(result$mapped)))
    )
  })

  output$comparative_ortholog_table <- DT::renderDT({
    result <- comparative_result()
    req(result)
    if (!nrow(result$mapped)) return(NULL)
    DT::datatable(result$mapped, rownames = FALSE, selection = "none",
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$comparative_availability_table <- DT::renderDT({
    DT::datatable(comparative_availability_display(availability()),
                  rownames = FALSE, selection = "none",
                  options = list(dom = "t", scrollX = TRUE))
  })

  output$comparative_result_table <- DT::renderDT({
    result <- comparative_result()
    req(result)
    if (is.null(result$results) || !nrow(result$results)) return(NULL)
    DT::datatable(result$results, rownames = FALSE, selection = "none",
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  group_genes <- eventReactive(input$comparative_group_analyze, {
    pasted <- unlist(strsplit(input$comparative_group_text %||% "", "\\s|,|;"))
    pasted <- pasted[nzchar(trimws(pasted))]
    uploaded <- character()
    file <- input$comparative_group_file
    if (!is.null(file$datapath) && file.exists(file$datapath)) {
      table <- tryCatch(read.csv(file$datapath, stringsAsFactors = FALSE,
                                 check.names = FALSE), error = function(e) data.frame())
      if (nrow(table) && ncol(table)) uploaded <- as.character(table[[1]])
    }
    unique(trimws(c(pasted, uploaded)))
  }, ignoreInit = TRUE)

  output$comparative_group_state <- renderUI({
    if ((input$comparative_group_analyze %||% 0L) == 0L)
      tags$p(class = "text-muted",
             "Upload or paste a gene list to summarize mapped orthogroups.")
  })

  output$comparative_group_summary <- renderUI({
    genes <- group_genes()
    req(genes)
    maps <- lapply(genes, find_orthogroup_for_gene,
                   orthology = orthology(), gene_lookup = gene_lookup())
    mapped <- sum(vapply(maps, nrow, integer(1)) > 0L)
    tags$div(
      class = "ytab-stat-grid",
      tags$div(tags$b(length(genes)), "Parsed genes"),
      tags$div(tags$b(mapped), "Genes with orthogroup mapping"),
      tags$div(tags$b(length(genes) - mapped), "Unmapped genes")
    )
  })

  output$comparative_group_unmapped <- DT::renderDT({
    genes <- group_genes()
    req(genes)
    mapped <- vapply(genes, function(gene)
      nrow(find_orthogroup_for_gene(gene, orthology(), gene_lookup())) > 0L,
      logical(1))
    data <- data.frame(gene = genes[!mapped], stringsAsFactors = FALSE)
    if (!nrow(data)) return(NULL)
    DT::datatable(data, rownames = FALSE, selection = "none")
  })

  output$comparative_group_availability <- DT::renderDT({
    DT::datatable(comparative_availability_display(availability()),
                  rownames = FALSE, selection = "none",
                  options = list(dom = "t", scrollX = TRUE))
  })

  output$comparative_group_results <- DT::renderDT({
    genes <- group_genes()
    req(genes)
    species <- input$comparative_group_species %||% ""
    projects <- project_inventory()
    subset <- projects[projects$species == species, , drop = FALSE]
    if (!nrow(subset)) return(NULL)
    data <- comparative_load_metric(subset$project_config[[1]],
                                    input$comparative_group_metric %||% "essentiality")
    rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L,
                   lapply(genes, function(gene) comparative_display_rows(data, gene)))
    rows <- if (length(rows)) do.call(rbind, rows) else data.frame()
    if (is.null(rows) || !nrow(rows)) return(NULL)
    DT::datatable(rows, rownames = FALSE, selection = "none",
                  options = list(pageLength = 10, scrollX = TRUE))
  })
}
