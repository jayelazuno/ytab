gene_domain_explorer_server <- function(input, output, session, active,
                                        active_project_path, repo_root,
                                        python_bin, log_text) {
  search_state <- reactiveVal(list(status = "idle", candidates = list()))
  selected_gene <- reactiveVal("")
  plot_state <- reactiveVal(NULL)
  track_data <- reactiveVal(data.frame())

  ordered_tracks <- function(data) gene_domain_order_tracks(data, active()$samples)
  refresh_preset_choices <- function(data) {
    data <- ordered_tracks(data)
    choices <- gene_domain_preset_choices(data)
    selected <- input$gene_domain_track_preset %||% "all"
    if (!selected %in% unname(choices)) selected <- "all"
    updateSelectInput(session, "gene_domain_track_preset", choices = choices, selected = selected)
    selected
  }
  sync_sample_selector <- function(preset = input$gene_domain_track_preset %||% "all") {
    data <- ordered_tracks(track_data())
    rows <- gene_domain_preset_track_rows(preset, data, input$gene_domain_samples %||% character())
    if (identical(preset, "custom")) {
      rows <- data
      selected <- intersect(input$gene_domain_samples %||% character(), rows$sample)
    } else {
      selected <- rows$sample
    }
    choices <- if (nrow(rows)) as.list(setNames(rows$sample, rows$track_name)) else list()
    updateSelectizeInput(session, "gene_domain_samples",
                         choices = choices, selected = selected, server = TRUE)
  }

  output$gene_domain_project_indicator <- renderUI({
    req(active())
    p <- active()
    tags$p(class = "text-muted", "Project: ", tags$b(p$project_id),
           " · Species: ", tags$b(p$species))
  })

  refresh_tracks <- function() {
    req(active_project_path(), active())
    data <- tryCatch(
      gene_domain_tracks(python_bin, repo_root, active_project_path(),
                         input$gene_domain_track_source %||% "raw"),
      error = function(e) data.frame()
    )
    data <- ordered_tracks(data)
    track_data(data)
    preset <- refresh_preset_choices(data)
    sync_sample_selector(preset)
  }

  observeEvent(active_project_path(), refresh_tracks(), ignoreInit = FALSE)
  observeEvent(input$gene_domain_track_source, refresh_tracks(), ignoreInit = TRUE)
  observeEvent(input$gene_domain_track_preset, {
    sync_sample_selector(input$gene_domain_track_preset %||% "all")
  }, ignoreInit = TRUE)

  output$gene_domain_preset_message <- renderUI({
    req(active())
    data <- ordered_tracks(track_data())
    if (!nrow(data)) return(tags$p(class = "text-muted", "No insertion tracks are available."))
    pairs <- vapply(as.character(1:4), function(pool)
      any(data$role == "parent" & data$pool == pool) &&
        any(data$role == "treated" & data$pool == pool), logical(1))
    if (!any(pairs))
      tags$p(class = "text-muted", "Matched pool presets require parent/treated and pool metadata.")
    else NULL
  })

  observeEvent(input$gene_domain_search, {
    req(active_project_path())
    query <- trimws(input$gene_domain_query %||% "")
    if (!nzchar(query)) {
      search_state(list(status = "no_query", candidates = list()))
      selected_gene("")
      return()
    }
    result <- tryCatch(
      gene_domain_json(python_bin, repo_root, active_project_path(),
                       c("--query", query)),
      error = function(e) list(status = "failure", error = conditionMessage(e))
    )
    search_state(result)
    candidates <- result$candidates %||% list()
    if (identical(result$status, "resolved") && length(candidates) == 1L)
      selected_gene(as.character(candidates[[1]]$gene_id %||% ""))
    else
      selected_gene("")
    log_text(paste(c(result$command %||% "", result$stdout %||% "",
                     result$stderr %||% "", paste("Exit status:", result$exit_status %||% "")),
                   collapse = "\n"))
  }, ignoreInit = TRUE)

  output$gene_domain_status <- renderUI({
    state <- search_state()
    status <- as.character(state$status %||% "idle")
    if (identical(status, "idle"))
      return(tags$p(class = "text-muted", "Enter a gene name or ID, then click Search gene."))
    if (identical(status, "no_query"))
      return(tags$p(class = "ytab-warning", "Enter a gene query first."))
    if (identical(status, "no_match"))
      return(tags$p(class = "ytab-warning", "No matching gene found."))
    if (identical(status, "ambiguous"))
      return(tags$p(class = "text-muted", "Multiple genes matched. Select one row before generating a figure."))
    if (identical(status, "failure"))
      return(tags$div(class = "alert alert-danger", state$error %||% "Gene search failed."))
    tags$p(class = "text-muted", "Gene resolved. Configure tracks and generate a figure.")
  })

  output$gene_domain_candidates <- DT::renderDT({
    data <- gene_domain_candidate_table(search_state()$candidates)
    if (!nrow(data)) return(NULL)
    DT::datatable(data, rownames = FALSE, selection = "single",
                  options = list(pageLength = 8, scrollX = TRUE))
  })

  observeEvent(input$gene_domain_candidates_rows_selected, {
    data <- gene_domain_candidate_table(search_state()$candidates)
    ix <- input$gene_domain_candidates_rows_selected
    if (length(ix) == 1L && nrow(data) >= ix)
      selected_gene(as.character(data$gene_id[[ix]]))
  }, ignoreInit = TRUE)

  output$gene_domain_gene_summary <- renderUI({
    req(active())
    gene_id <- selected_gene()
    if (!nzchar(gene_id)) return(NULL)
    data <- gene_domain_candidate_table(search_state()$candidates)
    row <- data[data$gene_id == gene_id, , drop = FALSE]
    if (!nrow(row)) return(tags$p("Selected gene: ", tags$b(gene_id)))
    tags$dl(
      class = "ytab-meta",
      tags$dt("Selected gene"), tags$dd(tags$b(row$display_name[[1]]), " · ", row$gene_id[[1]]),
      tags$dt("Location"), tags$dd(sprintf("%s:%s-%s (%s)", row$chromosome[[1]], row$start[[1]], row$end[[1]], row$strand[[1]])),
      tags$dt("Product"), tags$dd(row$product[[1]] %||% "")
    )
  })

  observeEvent(input$gene_domain_generate, {
    req(active_project_path())
    gene_id <- selected_gene()
    if (!nzchar(gene_id)) {
      showNotification("Search and resolve a gene before generating a figure.", type = "warning")
      return()
    }
    samples <- input$gene_domain_samples %||% character()
    args <- c(
      "--gene", gene_id,
      "--track-source", input$gene_domain_track_source %||% "raw",
      "--track-preset", input$gene_domain_track_preset %||% "all",
      "--samples", if (length(samples)) paste(samples, collapse = ",") else "all",
      "--flank-bp", as.character(as.integer(input$gene_domain_flank_bp %||% 1000L)),
      "--width-px", as.character(as.integer(input$gene_domain_width_px %||% 1800L)),
      "--label-mode", input$gene_domain_label_mode %||% "full",
      if (!isTRUE(input$gene_domain_show_site_counts)) "--hide-site-counts",
      if (!isTRUE(input$gene_domain_show_domains)) "--hide-domains",
      if (!isTRUE(input$gene_domain_show_direction)) "--hide-direction"
    )
    result <- tryCatch(
      gene_domain_json(python_bin, repo_root, active_project_path(), args),
      error = function(e) list(status = "failure", error = conditionMessage(e))
    )
    plot_state(result)
    log_text(paste(c(result$command %||% "", result$stdout %||% "",
                     result$stderr %||% "", paste("Exit status:", result$exit_status %||% "")),
                   collapse = "\n"))
    if (identical(result$status, "failure")) {
      showNotification(result$error %||% "Figure generation failed.", type = "error", duration = NULL)
    } else {
      showNotification("Gene insertion figure is ready.", type = "message")
    }
  }, ignoreInit = TRUE)

  output$gene_domain_domain_status <- renderUI({
    manifest <- plot_state()$manifest %||% NULL
    if (is.null(manifest)) return(NULL)
    if (isTRUE(manifest$domains_available))
      tags$p(class = "text-muted", "Domain source: ", manifest$domain_source)
    else
      tags$p(class = "text-muted", "No domain annotation available for this gene.")
  })

  output$gene_domain_figure_ui <- renderUI({
    state <- plot_state()
    if (is.null(state)) {
      return(ytab_empty_state(
        "Search for a gene, choose tracks, then click Generate figure.",
        "The figure preview and downloads appear here after a successful run."
      ))
    }
    if (identical(state$status, "failure")) {
      return(ytab_result_card(
        "Gene Explorer could not generate a figure.",
        status = "warning",
        tags$p("Reason: ", state$error %||% "Figure generation failed."),
        tags$details(class = "ytab-technical-details",
                     tags$summary("Show technical details"),
                     tags$div(class = "ytab-technical-console",
                              tags$pre(paste(c(state$command %||% "", state$stdout %||% "", state$stderr %||% ""),
                                             collapse = "\n"))))
      ))
    }
    if (!gene_domain_valid_output(state)) {
      return(ytab_empty_state(
        "No valid Gene Explorer output is available yet.",
        "Search for a gene, choose tracks, then click Generate figure."
      ))
    }
    ytab_static_image_output_card(
      "Generated: Gene/domain insertion figure",
      "gene_domain_figure",
      description = "Python-generated PNG from existing insertion tracks; displayed as a static image with preserved aspect ratio.",
      height = input$gene_domain_display_height %||% "large",
      downloads = tagList(downloadButton("gene_domain_download_png", "Download PNG"),
                          downloadButton("gene_domain_download_table", "Download table"))
    )
  })

  output$gene_domain_figure <- renderImage({
    manifest <- plot_state()$manifest %||% NULL
    if (is.null(manifest)) return(NULL)
    src <- as.character(manifest$figure_path %||% "")
    if (!nzchar(src) || !file.exists(src)) return(NULL)
    list(src = src, contentType = "image/png", alt = "Gene/domain insertion figure")
  }, deleteFile = FALSE)

  output$gene_domain_insertions <- DT::renderDT({
    manifest <- plot_state()$manifest %||% NULL
    if (is.null(manifest)) return(NULL)
    data <- gene_domain_read_insertions(manifest)
    if (!nrow(data)) return(DT::datatable(data.frame(Message = "No insertions in the selected region."), rownames = FALSE))
    DT::datatable(data, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$gene_domain_technical <- renderText({
    state <- plot_state() %||% search_state()
    paste(c(state$command %||% "", state$stdout %||% "", state$stderr %||% ""), collapse = "\n")
  })

  output$gene_domain_download_png <- downloadHandler(
    filename = function() basename(as.character(plot_state()$manifest$figure_path %||% "gene_domain_insertions.png")),
    content = function(file) {
      validate(need(gene_domain_valid_output(plot_state()), "No Gene Explorer PNG is available for download."))
      file.copy(as.character(plot_state()$manifest$figure_path), file, overwrite = TRUE)
    }
  )

  output$gene_domain_download_table <- downloadHandler(
    filename = function() basename(as.character(plot_state()$manifest$table_path %||% "gene_domain_insertions.csv")),
    content = function(file) {
      validate(need(gene_domain_valid_output(plot_state()), "No Gene Explorer table is available for download."))
      file.copy(as.character(plot_state()$manifest$table_path), file, overwrite = TRUE)
    }
  )
}
