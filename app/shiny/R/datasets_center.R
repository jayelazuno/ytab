# Datasets page: the public Tn-seq studies YTAB has reprocessed, from resources/datasets/catalog.tsv and libraries.tsv
# (written by `ytab datasets`); the app only reads them. Each row is one study in one project.

datasets_catalog <- function(repo_root) {
  read <- function(name) {
    path <- file.path(repo_root, "resources", "datasets", name)
    if (!file.exists(path)) return(data.frame())
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", na.strings = character(), colClasses = "character")
  }
  catalog <- read("catalog.tsv"); libraries <- read("libraries.tsv")
  if (nrow(catalog)) {
    # One row per study, not per project: a study whose strains need different references is processed as separate
    # projects (a project has one reference), but it is still one study and is listed once.
    catalog$id <- catalog$study
    catalog$open_config <- ifelse(file.exists(file.path(repo_root, catalog$project_config)), catalog$project_config, "")
  }
  list(catalog = catalog, libraries = libraries)
}

datasets_number <- function(x, digits = 1) { v <- suppressWarnings(as.numeric(x)); ifelse(is.na(v), "", formatC(v, format = "f", digits = digits)) }

# Search: catalogue rows whose fields, or whose libraries' IDs, names and conditions, contain every query word.
datasets_search <- function(data, words) {
  x <- data$catalog; empty <- data.frame(kind = character(), id = character(), title = character(), section = character(), snippet = character(), stringsAsFactors = FALSE)
  if (!nrow(x) || !length(words)) return(empty)
  libs <- data$libraries
  text <- vapply(seq_len(nrow(x)), function(i) {
    l <- libs[libs$study == x$study[[i]] & libs$project == x$project[[i]], intersect(c("sample", "display_name", "condition", "background", "notes"), names(libs)), drop = FALSE]
    tolower(paste(c(unlist(x[i, setdiff(names(x), c("id", "open_config"))]), unlist(l)), collapse = " "))
  }, "")
  hit <- vapply(text, function(t) all(vapply(words, grepl, logical(1), x = t, fixed = TRUE)), logical(1))
  if (!any(hit)) return(empty)
  x <- x[hit, , drop = FALSE]
  data.frame(kind = "dataset", id = x$id, title = sprintf("%s: %s, %s", x$citation, x$species, x$strain), section = "Dataset",
             snippet = sprintf("%s. %s libraries; %s. %s", x$transposon, x$runs_included, x$accession, x$design), stringsAsFactors = FALSE)
}

datasets_modal <- function() modalDialog(title = "Public datasets", size = "xl", easyClose = TRUE, footer = modalButton("Close"),
  tags$div(class = "ytab-datasets",
    tags$p(class = "text-muted", "Published Hermes and Ac/Ds libraries reprocessed with YTAB from raw reads. Search any field, including drug, strain or run ID; select a row for its conditions, QC, results and libraries. Sites per kb: median insertion density of the analysed libraries; top site: the largest share of a library's reads at one site."),
    DT::DTOutput("datasets_table"), uiOutput("datasets_detail")))

# Catalog rows of one study, in the order they appear.
datasets_studies <- function(x) if (!nrow(x)) list() else split(x, factor(x$study, levels = unique(x$study)))

# One study's rows as a single row: its projects joined, its libraries summed, its conditions and results unioned.
datasets_merge_group <- function(g) {
  if (nrow(g) == 1) return(g)
  join <- function(column) paste(unique(unlist(strsplit(g[[column]], "; ", fixed = TRUE))), collapse = "; ")
  g$runs <- as.character(sum(suppressWarnings(as.integer(g$runs)), na.rm = TRUE))
  g$runs_included <- as.character(sum(suppressWarnings(as.integer(g$runs_included)), na.rm = TRUE))
  g$conditions <- join("conditions"); g$in_app <- join("in_app"); g$backgrounds <- join("backgrounds")
  g$project <- paste(g$project, collapse = " ")
  g[1, , drop = FALSE]
}

datasets_table_data <- function(x, libraries) {
  if (!nrow(x)) return(data.frame())
  groups <- datasets_studies(x)
  x <- do.call(rbind, lapply(groups, datasets_merge_group))
  study <- ifelse(nzchar(x$doi), sprintf('<a href="https://doi.org/%s" target="_blank" rel="noopener">%s</a>', x$doi, htmltools::htmlEscape(x$citation)), htmltools::htmlEscape(x$citation))
  n_conditions <- vapply(strsplit(x$conditions, "; ", fixed = TRUE), length, 1L)
  short <- function(v) paste(unique(sub(" \\(.*", "", sub("^Treated vs control.*", "Treated vs control", strsplit(v, "; ", fixed = TRUE)[[1]]))), collapse = ", ")
  # hidden column: conditions, strains and every run ID and name, so the table's search finds them
  hidden <- vapply(seq_len(nrow(x)), function(i) {
    l <- libraries[libraries$study == x$study[[i]] & libraries$project == x$project[[i]], , drop = FALSE]
    paste(x$conditions[[i]], x$backgrounds[[i]], x$design[[i]], x$project[[i]], paste(l$sample, l$display_name, collapse = " "))
  }, "")
  data.frame(Study = study, Species = sprintf("<i>%s</i>", htmltools::htmlEscape(x$species)), Strain = htmltools::htmlEscape(x$strain),
             Transposon = htmltools::htmlEscape(x$transposon), Libraries = suppressWarnings(as.integer(x$runs_included)),
             Conditions = n_conditions, Accession = htmltools::htmlEscape(x$accession),
             `Sites per kb` = suppressWarnings(as.numeric(datasets_number(x$median_sites_per_kb))),
             `Top site (%)` = suppressWarnings(round(100 * as.numeric(x$max_top_site_share), 1)),
             `In the app` = vapply(x$in_app, short, ""), search = hidden, check.names = FALSE, stringsAsFactors = FALSE)
}

datasets_detail_ui <- function(row, libraries) {
  field <- function(label, value) if (inherits(value, "shiny.tag") || (length(value) == 1L && !is.na(value) && nzchar(value))) tags$tr(tags$th(label), tags$td(value))
  link <- function(url, text) tags$a(href = url, target = "_blank", rel = "noopener", text)
  projects <- strsplit(row$project, " ")[[1]]
  libs <- libraries[libraries$study == row$study & libraries$project %in% projects, , drop = FALSE]
  share <- suppressWarnings(as.numeric(row$max_top_site_share))
  tags$div(class = "ytab-dataset-detail",
    tags$h4(row$title),
    tags$p(class = "text-muted", sprintf("%s (%s) %s.", row$citation, row$year, row$journal), " ",
           if (nzchar(row$doi)) link(paste0("https://doi.org/", row$doi), paste0("doi:", row$doi)), " ",
           if (nzchar(row$pmid)) link(paste0("https://pubmed.ncbi.nlm.nih.gov/", row$pmid), paste0("PMID ", row$pmid))),
    tags$table(class = "table table-sm ytab-dataset-fields", tags$tbody(
      field("Species", tags$i(row$species)), field("Strain", row$strain), field("Transposon", row$transposon), field("Read design", row$read_design),
      field("Accession", row$accession), field("Design", row$design),
      field("Libraries", sprintf("%s runs (%s analysed): %s", row$runs, row$runs_included, row$roles)),
      field("Strains in the project", row$backgrounds), field("Conditions", row$conditions),
      field("Density", if (nzchar(row$median_sites_per_kb)) sprintf("median %s insertion sites per kb across the analysed libraries", datasets_number(row$median_sites_per_kb))),
      field("Top site", if (!is.na(share)) sprintf("up to %.1f%% of a library's reads at one site (%s)", 100 * share, row$max_top_site)),
      field("Results in the app", row$in_app),
      field(if (length(projects) > 1) "Projects" else "Project",
            if (length(projects) > 1)
              tagList(lapply(projects, function(p) tags$div(tags$code(p), " - ",
                sprintf("%d libraries", sum(libs$project == p & nzchar(libs$library))))),
                tags$p(class = "text-muted", "One study, processed as separate projects because its strains are mapped to different reference genomes; a project has one reference."))
            else tags$code(row$project)),
      field("Experimental design", "Which runs make each library and which analyses use it: Methods > Experimental designs of the public datasets, or ytab datasets --find RUN"))),
    tags$div(class = "ytab-actions",
      if (nzchar(row$open_config))
        tagList(lapply(seq_along(projects), function(i) actionButton(paste0("datasets_open_project_", i),
          if (length(projects) > 1) paste("Open", projects[[i]]) else "Open project", class = "btn-primary ytab-button")))
      else tags$span(class = "text-muted", "This project is not on this computer; download the YTAB projects bundle to open it.")),
    if (nrow(libs)) tagList(tags$h5("Libraries"), DT::DTOutput("datasets_libraries")))
}

datasets_server <- function(input, output, session, repo_root, data, open_project = NULL) {
  # The table shows one row per study, so a clicked row index must be read against the grouped rows, not the catalog.
  grouped <- do.call(rbind, lapply(datasets_studies(data$catalog), datasets_merge_group))
  selected <- reactiveVal("")
  open_datasets <- function(id = "") { selected(id); showModal(datasets_modal()) }
  observeEvent(input$open_datasets, open_datasets(), ignoreInit = TRUE)
  observeEvent(input$open_datasets_nav, open_datasets(), ignoreInit = TRUE)
  output$datasets_table <- DT::renderDT({
    x <- grouped; req(nrow(x)); pre <- match(isolate(selected()), x$id)
    DT::datatable(datasets_table_data(data$catalog, data$libraries), escape = FALSE, rownames = FALSE, filter = "none",
      selection = list(mode = "single", selected = if (!is.na(pre)) pre else NULL, target = "row"), class = "compact stripe hover",
      options = list(pageLength = 15, lengthMenu = c(15, 50), scrollX = TRUE, autoWidth = FALSE,
                     search = list(search = ""), dom = "ftip", columnDefs = list(list(className = "dt-right", targets = c(4, 5, 7, 8)), list(visible = FALSE, targets = 10))))
  })
  observeEvent(input$datasets_table_rows_selected, { i <- input$datasets_table_rows_selected; if (length(i)) selected(grouped$id[[i]]) }, ignoreInit = TRUE)
  current <- reactive({ x <- data$catalog; g <- x[x$id == selected(), , drop = FALSE]; if (!nrow(g)) NULL else datasets_merge_group(g) })
  output$datasets_detail <- renderUI({ row <- current(); if (is.null(row)) return(NULL); datasets_detail_ui(as.list(row), data$libraries) })
  current_projects <- reactive({ row <- current(); if (is.null(row)) character() else strsplit(row$project, " ")[[1]] })
  output$datasets_libraries <- DT::renderDT({
    row <- current(); req(row); l <- data$libraries; l <- l[l$study == row$study & l$project == row$project, , drop = FALSE]
    t <- data.frame(Run = l$sample, Name = l$display_name, Role = l$role, Condition = l$condition, Strain = l$background, Pool = l$pool,
                    Library = if ("library" %in% names(l)) ifelse(nzchar(l$library), l$library, "not analysed") else "",
                    `Used in` = if ("used_in" %in% names(l)) gsub("; ", "<br>", htmltools::htmlEscape(l$used_in)) else "", Reads = suppressWarnings(as.numeric(l$input_reads)),
                    `Mapped (%)` = suppressWarnings(as.numeric(datasets_number(l$percent_mapped))), `Insertion sites` = suppressWarnings(as.numeric(l$insertion_sites)),
                    `Sites per kb` = suppressWarnings(as.numeric(datasets_number(l$sites_per_kb))),
                    `Top site share (%)` = suppressWarnings(round(100 * as.numeric(l$top_site_share), 1)), check.names = FALSE)
    t <- DT::datatable(t, rownames = FALSE, filter = "top", selection = "none", class = "compact stripe hover", escape = setdiff(seq_along(t), which(names(t) == "Used in")),
                  options = list(pageLength = 10, scrollX = TRUE, dom = "ftip"))
    DT::formatRound(t, c("Reads", "Insertion sites"), digits = 0)
  })
  for (i in 1:4) local({ k <- i
    observeEvent(input[[paste0("datasets_open_project_", k)]], {
      p <- current_projects(); req(length(p) >= k, is.function(open_project)); removeModal()
      open_project(normalizePath(file.path(repo_root, "output/projects", p[[k]], "config/project.yaml"), winslash = "/"))
    }, ignoreInit = TRUE)
  })
  list(open = open_datasets)
}
