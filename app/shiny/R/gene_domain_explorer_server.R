gene_domain_explorer_server <- function(input, output, session, active,
                                        active_project_path, repo_root,
                                        python_bin, log_text,
                                        projects = reactive(data.frame()), reference_state = NULL) {
  if (is.null(reference_state)) reference_state <- ytab_reference_state()
  ytab_reference_server("gene_domain", input, output, session, repo_root, active, projects, reference_state,
                        display_modes = FALSE)
  search_state <- reactiveVal(list(status = "idle", candidates = list()))
  selected_gene <- reactiveVal("")
  plot_state <- reactiveVal(NULL)
  track_data <- reactiveVal(data.frame())
  custom_track_selection <- reactiveVal(character())

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
      prior <- custom_track_selection()
      selected <- if (length(prior)) intersect(prior, rows$sample) else
        intersect(input$gene_domain_samples %||% character(), rows$sample)
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
    tags$p(class = "text-muted", "Dataset: ", tags$b(project_display_label(p)),
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
  observeEvent(input$gene_domain_samples, {
    if (identical(input$gene_domain_track_preset %||% "all", "custom")) {
      custom_track_selection(input$gene_domain_samples %||% character())
    }
  }, ignoreInit = TRUE)

  output$gene_domain_preset_message <- renderUI({
    req(active())
    data <- ordered_tracks(track_data())
    if (!nrow(data)) return(tags$p(class = "text-muted", "No insertion tracks are available."))
    preset <- input$gene_domain_track_preset %||% "all"
    rows <- gene_domain_preset_track_rows(preset, data, input$gene_domain_samples %||% character())
    if (!identical(preset, "custom") && !nrow(rows))
      return(tags$p(class = "ytab-warning", "No tracks match this preset for the current dataset."))
    if (identical(preset, "matched_pairs") && !nrow(rows))
      return(tags$p(class = "ytab-warning", "No matched parent/treated pool pairs were detected. Use All tracks or Custom."))
    if (!length(gene_domain_matched_pools(data)))
      tags$p(class = "text-muted", "Matched pool presets require parent/treated and pool metadata.")
    else
      tags$p(class = "text-muted", "Choose which insertion tracks are displayed for the selected gene.")
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
      return(tags$p(class = "text-muted", "Enter a gene name or ID, then click Search gene. The gene and its figure appear together."))
    if (identical(status, "no_query"))
      return(tags$p(class = "ytab-warning", "Enter a gene query first."))
    if (identical(status, "no_match"))
      return(tags$p(class = "ytab-warning", "No matching gene found."))
    if (identical(status, "ambiguous"))
      return(tags$p(class = "text-muted", "Multiple genes matched. Choose one under Matching genes."))
    if (identical(status, "failure"))
      return(tags$div(class = "alert alert-danger", state$error %||% "Gene search failed."))
    tags$p(class = "text-muted", "Gene resolved; the figure is below.")
  })

  output$gene_domain_candidates <- DT::renderDT({
    data <- gene_domain_candidate_table(search_state()$candidates)
    if (!nrow(data)) return(NULL)
    data <- ytab_join_glabrata_display(data, repo_root, "gene_id")
    if ("chromosome" %in% names(data)) {
      data$chromosome_display <- qc_plot_chromosome_display(data$chromosome, active())
    }
    mapped <- data.frame(`Gene ID` = data$cagl_display_id,
                         `Gene name` = data$gene_display_name,
                         `S. cerevisiae orthology` = data$cg_to_sc_relationship_display,
                         data[, intersect(c("chromosome_display", "start", "end", "strand",
                                             "product", "match_type"), names(data)),
                              drop = FALSE],
                         ytab_glabrata_gene_detail_columns(data),
                         check.names = FALSE)
    names(mapped)[names(mapped) == "chromosome_display"] <- "Chromosome"
    data <- mapped
    ytab_gene_details_datatable(data, page_length = 8, filter = "none",
                                selection = "single",
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
    row <- ytab_join_glabrata_display(row, repo_root, "gene_id")
    chromosome_display <- qc_plot_chromosome_display(row$chromosome[[1]], active())
    tags$dl(
      class = "ytab-meta",
      tags$dt("Gene ID"), tags$dd(row$cagl_display_id[[1]]),
      tags$dt("Gene name"), tags$dd(tags$b(row$gene_display_name[[1]])),
      tags$dt("S. cerevisiae orthology"), tags$dd(row$cg_to_sc_relationship_display[[1]] %||% ""),
      tags$dt("Location"), tags$dd(sprintf("%s:%s-%s (%s)", chromosome_display, row$start[[1]], row$end[[1]], row$strand[[1]])),
      tags$dt("Product"), tags$dd(row$product[[1]] %||% "")
    )
  })

  # The figure is the answer to the search, not a second step: it is drawn as soon as a gene resolves, and redrawn
  # when the tracks, the style or a display option changes. Debounced because several of those inputs move together
  # (a preset change rewrites the sample list), and each draw is a Python process.
  figure_request <- debounce(reactive(list(
    gene = selected_gene(),
    samples = input$gene_domain_samples %||% character(),
    style = input$gene_domain_track_style %||% "sites",
    source = input$gene_domain_track_source %||% "raw",
    preset = input$gene_domain_track_preset %||% "all",
    flank = input$gene_domain_flank_bp %||% 1000L,
    width = input$gene_domain_width_px %||% 1800L,
    label_mode = input$gene_domain_label_mode %||% "full",
    site_counts = isTRUE(input$gene_domain_show_site_counts),
    domains = isTRUE(input$gene_domain_show_domains),
    direction = isTRUE(input$gene_domain_show_direction),
    references = reference_state() %||% character()
  )), 800)

  observeEvent(figure_request(), {
    request <- figure_request()
    if (!nzchar(request$gene %||% "")) { plot_state(NULL); return() }
    req(active_project_path())
    if (!length(request$samples)) {
      plot_state(list(status = "failure", error = "No tracks match this preset for the current dataset."))
      return()
    }
    args <- c(
      "--gene", request$gene,
      "--track-source", request$source,
      "--track-preset", request$preset,
      "--samples", paste(request$samples, collapse = ","),
      "--track-style", request$style,
      "--flank-bp", as.character(as.integer(request$flank)),
      "--width-px", as.character(as.integer(request$width)),
      "--label-mode", request$label_mode,
      if (!request$site_counts) "--hide-site-counts",
      if (!request$domains) "--hide-domains",
      if (!request$direction) "--hide-direction",
      # Each chosen project draws its own copy of the gene below this one's tracks. Build the flag/value pairs only
      # when there is at least one: rbind() of a flag and an empty vector still yields the bare flag, which the
      # parser rejects as a missing argument.
      if (length(request$references)) as.character(rbind("--reference-config", request$references))
    )
    result <- tryCatch(
      gene_domain_json(python_bin, repo_root, active_project_path(), args),
      error = function(e) list(status = "failure", error = conditionMessage(e))
    )
    plot_state(result)
    log_text(paste(c(result$command %||% "", result$stdout %||% "",
                     result$stderr %||% "", paste("Exit status:", result$exit_status %||% "")),
                   collapse = "\n"))
    if (identical(result$status, "failure"))
      showNotification(result$error %||% "Figure generation failed.", type = "error", duration = NULL)
  }, ignoreInit = FALSE)

  output$gene_domain_style_note <- renderUI({
    if (identical(input$gene_domain_track_style %||% "sites", "reads"))
      tags$p(class = "text-muted", "Each insertion site is a stem whose height is the reads at it, log-scaled ",
             "within its own track, so what compares between tracks of different depth is where the reads fall, ",
             "not how many. An essential gene reads as a gap over the coding sequence.")
    else
      tags$p(class = "text-muted", "Each insertion site is a tick of equal height, its width scaled by the reads ",
             "at it: the map of what was hit, with the site count per track in the label.")
  })

  output$gene_domain_domain_status <- renderUI({
    manifest <- plot_state()$manifest %||% NULL
    if (is.null(manifest)) return(NULL)
    if (isTRUE(manifest$domains_available))
      tags$p(class = "text-muted", "Domain source: ", manifest$domain_source)
    else
      tags$p(class = "text-muted", "No domain annotation available for this gene.")
  })

  has_figure <- reactive(gene_domain_valid_output(plot_state()))
  output$gene_domain_has_figure <- reactive(has_figure())
  # The flag drives a conditionalPanel, which Shiny would otherwise never evaluate while the panel is hidden.
  outputOptions(output, "gene_domain_has_figure", suspendWhenHidden = FALSE)

  output$gene_domain_figure_message <- renderUI({
    state <- plot_state()
    if (is.null(state)) {
      return(ytab_empty_state(
        "Search for a gene to see its insertions.",
        "The figure is drawn as soon as the gene resolves, and redrawn when you change the figure style, the tracks or another dataset to compare with."
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
    if (!has_figure())
      return(ytab_empty_state(
        "No valid Gene Explorer output is available yet.",
        "Search for a gene; if this persists, see the technical details below."
      ))
    NULL
  })

  output$gene_domain_figure_caption <- renderUI({
    req(has_figure())
    tags$p(class = "text-muted",
           if (identical(input$gene_domain_display_height %||% "fit", "fit"))
             "Scaled to fit the window; choose a preview height under Display, or download the PNG, to see it at full size."
           else "Shown at the chosen preview height; scroll the panel to see the whole figure.")
  })

  output$gene_domain_figure <- renderImage({
    manifest <- plot_state()$manifest %||% NULL
    if (is.null(manifest)) return(NULL)
    src <- as.character(manifest$figure_path %||% "")
    if (!nzchar(src) || !file.exists(src)) return(NULL)
    # A tall figure (one lane per track, plus any comparison datasets) overflows a fixed-height panel, so the
    # default caps it against the viewport and the browser scales it down whole rather than making the user scroll.
    choice <- input$gene_domain_display_height %||% "fit"
    list(src = src, contentType = "image/png", alt = "Gene/domain insertion figure",
         class = if (identical(choice, "fit")) "ytab-image-fit" else NULL,
         style = if (identical(choice, "fit")) NULL else paste0("max-height:", ytab_plot_height_px(choice), ";"))
  }, deleteFile = FALSE)

  output$gene_domain_insertions <- DT::renderDT({
    manifest <- plot_state()$manifest %||% NULL
    if (is.null(manifest)) return(NULL)
    data <- gene_domain_read_insertions(manifest)
    if (!nrow(data)) return(DT::datatable(data.frame(Message = "No insertions in the selected region."), rownames = FALSE))
    if ("chromosome" %in% names(data)) {
      data$chromosome <- qc_plot_chromosome_display(data$chromosome, active())
      names(data)[names(data) == "chromosome"] <- "Chromosome"
    }
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
      data <- gene_domain_read_insertions(plot_state()$manifest)
      if ("chromosome" %in% names(data)) {
        data$chromosome <- qc_plot_chromosome_display(data$chromosome, active())
        names(data)[names(data) == "chromosome"] <- "Chromosome"
      }
      write.csv(data, file, row.names = FALSE)
    }
  )
}
