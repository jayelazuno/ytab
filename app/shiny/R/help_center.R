# Help center: methods pages, tutorial, release notes and search. Content is Markdown in docs/ (docs/methods/index.tsv
# lists the methods pages); the app only renders it. Images use paths relative to their page (e.g. img/qc.png), served
# under the resource path "ytab-docs".

help_pages <- function(repo_root) {
  index <- file.path(repo_root, "docs", "methods", "index.tsv")
  pages <- if (file.exists(index)) utils::read.delim(index, stringsAsFactors = FALSE) else data.frame(order = integer(), section = character(), title = character(), file = character())
  extra <- data.frame(order = c(0L, 999L), section = c("Getting started", "About"), title = c("Tutorial: Gale 2020 from project to results", "Release notes"),
                      file = c("docs/methods/tutorial.md", "docs/RELEASE_NOTES.md"), stringsAsFactors = FALSE)
  pages <- rbind(extra[1, ], pages, extra[2, ])
  pages <- pages[file.exists(file.path(repo_root, pages$file)), , drop = FALSE]
  pages$id <- make.names(sub("\\.md$", "", basename(pages$file)), unique = TRUE)
  rownames(pages) <- NULL
  pages
}

help_markdown <- function(repo_root, file) paste(readLines(file.path(repo_root, file), warn = FALSE, encoding = "UTF-8"), collapse = "\n")

help_render <- function(repo_root, file) {
  html <- commonmark::markdown_html(help_markdown(repo_root, file), extensions = TRUE)
  sub_dir <- sub("^docs/?", "", dirname(file)); base <- if (nzchar(sub_dir)) paste0("ytab-docs/", sub_dir) else "ytab-docs"
  # relative image paths resolve against the page's folder
  html <- gsub('<img src="(?!https?:|/|ytab-docs/)([^"]+)"', paste0('<img class="ytab-help-image" src="', base, '/\\1"'), html, perl = TRUE)
  HTML(html)
}

# Plain text of a page for search (Markdown syntax stripped).
help_plain_text <- function(repo_root, file) {
  text <- help_markdown(repo_root, file)
  text <- gsub("```[\\s\\S]*?```", " ", text, perl = TRUE)
  text <- gsub("!\\[[^]]*\\]\\([^)]*\\)", " ", text)
  text <- gsub("\\[([^]]*)\\]\\([^)]*\\)", "\\1", text)
  text <- gsub("[#*`|>_]", " ", text)
  gsub("[ \t]+", " ", text)
}

# Tabs of the workspace and words that lead to them (for search).
help_tab_index <- data.frame(
  value = c("project", "run", "quality_control", "normalization", "fitness", "essentiality", "gene_explorer", "screen_coverage"),
  title = c("Project", "Run", "Quality Control", "Normalization", "Condition Comparisons", "Essentiality", "Gene Explorer", "Screen Planning"),
  keywords = c("samples sample sheet provenance project setup",
               "mapping fastq hit files summarytable alignment preprocessing",
               "qc library quality saturation midlc sequence bias insertions density mapping statistics",
               "normalization bgc beta-geometric ttr total reads quantile zinb midlc scaling",
               "fitness treated control z-score resampling gumbel depleted enriched ma plot volcano top hits comparison design",
               "essential essentiality hmm hybrid posterior chromosome track classifier",
               "gene domain insertions genome browser igv go enrichment across species ortholog",
               "screen planning coverage colonies poisson mutant identification sanger arbitrary pcr"),
  stringsAsFactors = FALSE)

# Search: tabs whose name or keywords match, then pages ranked by the number of query words they contain.
help_search <- function(repo_root, query, pages = help_pages(repo_root), max_hits = 12L) {
  words <- unique(tolower(strsplit(trimws(query %||% ""), "\\s+")[[1]])); words <- words[nchar(words) >= 2L]
  empty <- data.frame(kind = character(), id = character(), title = character(), section = character(), snippet = character(), stringsAsFactors = FALSE)
  if (!length(words)) return(empty)
  tab_text <- tolower(paste(help_tab_index$title, help_tab_index$keywords))
  tab_hits <- help_tab_index[vapply(tab_text, function(t) all(vapply(words, grepl, logical(1), x = t, fixed = TRUE)), logical(1)), , drop = FALSE]
  tabs <- data.frame(kind = rep("tab", nrow(tab_hits)), id = tab_hits$value, title = tab_hits$title, section = rep("Go to tab", nrow(tab_hits)),
                     snippet = rep("", nrow(tab_hits)), stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(pages)), function(i) {
    text <- help_plain_text(repo_root, pages$file[[i]]); lower <- tolower(paste(pages$title[[i]], text))
    found <- vapply(words, grepl, logical(1), x = lower, fixed = TRUE)
    if (!any(found)) return(NULL)
    at <- regexpr(words[found][[1]], tolower(text), fixed = TRUE)
    snippet <- if (at > 0) paste0("…", trimws(substr(text, max(1L, at - 80L), at + 160L)), "…") else ""
    data.frame(kind = "page", id = pages$id[[i]], title = pages$title[[i]], section = pages$section[[i]], snippet = gsub("\\s+", " ", snippet),
               score = sum(found) + 2 * all(found), stringsAsFactors = FALSE)
  })
  found <- do.call(rbind, Filter(Negate(is.null), rows))
  if (!is.null(found) && nrow(found)) found <- found[order(-found$score), setdiff(names(found), "score"), drop = FALSE] else found <- empty
  utils::head(rbind(tabs, found), max_hits)
}

help_center_modal <- function(pages, selected) {
  sections <- unique(pages$section)
  nav <- lapply(sections, function(s) {
    rows <- pages[pages$section == s, , drop = FALSE]
    tags$div(class = "ytab-help-nav-section", tags$h5(s),
             lapply(seq_len(nrow(rows)), function(i) actionLink(paste0("help_page_", rows$id[[i]]), rows$title[[i]],
                                                                 class = if (identical(rows$id[[i]], selected)) "ytab-help-link active" else "ytab-help-link")))
  })
  modalDialog(title = "Methods and help", size = "l", easyClose = TRUE, footer = modalButton("Close"),
    tags$div(class = "ytab-help-center", tags$nav(class = "ytab-help-nav", nav), tags$article(class = "ytab-help-page", uiOutput("help_page_body"))))
}

help_search_modal <- function() modalDialog(title = "Search YTAB", size = "m", easyClose = TRUE, footer = modalButton("Close"),
  textInput("help_search_query", NULL, placeholder = "e.g. BGC, essential, top hits, colonies, IGV", width = "100%"),
  uiOutput("help_search_results"))

help_center_server <- function(input, output, session, repo_root, go_to = NULL) {
  if (dir.exists(file.path(repo_root, "docs"))) addResourcePath("ytab-docs", file.path(repo_root, "docs"))
  pages <- help_pages(repo_root); current <- reactiveVal(if (nrow(pages)) pages$id[[1]] else "")
  open_page <- function(id) { current(id); showModal(help_center_modal(pages, id)) }
  output$help_page_body <- renderUI({ id <- current(); req(nzchar(id)); help_render(repo_root, pages$file[match(id, pages$id)]) })
  for (i in seq_len(nrow(pages))) local({ id <- pages$id[[i]]; observeEvent(input[[paste0("help_page_", id)]], open_page(id), ignoreInit = TRUE) })
  first_method <- function() pages$id[match("Building a pool", pages$section)] %||% pages$id[[1]]
  observeEvent(input$open_help_methods, open_page(first_method()), ignoreInit = TRUE)
  observeEvent(input$open_help_methods_nav, open_page(first_method()), ignoreInit = TRUE)
  observeEvent(input$open_help_tutorial, open_page("tutorial"), ignoreInit = TRUE)
  observeEvent(input$open_help_release, open_page("RELEASE_NOTES"), ignoreInit = TRUE)
  observeEvent(input$open_search, showModal(help_search_modal()), ignoreInit = TRUE)
  observeEvent(input$open_search_nav, showModal(help_search_modal()), ignoreInit = TRUE)
  output$help_search_results <- renderUI({
    hits <- help_search(repo_root, input$help_search_query %||% "", pages)
    if (!nzchar(trimws(input$help_search_query %||% ""))) return(tags$p(class = "text-muted", "Search the methods, the tutorial, the release notes and the tabs."))
    if (!nrow(hits)) return(tags$p(class = "text-muted", "Nothing found."))
    tags$ul(class = "ytab-search-results", lapply(seq_len(nrow(hits)), function(i) tags$li(
      tags$a(href = "#", class = "ytab-search-hit", `data-kind` = hits$kind[[i]], `data-id` = hits$id[[i]], hits$title[[i]]),
      tags$span(class = "ytab-search-section", hits$section[[i]]), if (nzchar(hits$snippet[[i]])) tags$p(class = "ytab-search-snippet", hits$snippet[[i]]))),
      tags$script(HTML("$('.ytab-search-hit').on('click', function(e){ e.preventDefault(); Shiny.setInputValue('help_search_open', {kind: $(this).data('kind'), id: $(this).data('id')}, {priority: 'event'}); });")))
  })
  observeEvent(input$help_search_open, {
    hit <- input$help_search_open
    if (identical(hit$kind, "tab") && is.function(go_to)) { removeModal(); go_to(hit$id) } else open_page(hit$id)
  }, ignoreInit = TRUE)
  invisible(pages)
}
