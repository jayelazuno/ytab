essentiality_action_label <- function(mode, cached, force, preview_label, run_label,
                                      cached_label, rerun_label) {
  if (identical(mode, "preview")) return(preview_label)
  if (isTRUE(force)) return(rerun_label)
  if (isTRUE(cached)) cached_label else run_label
}

essentiality_relative_path <- function(path, project_root) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- paste0(normalizePath(project_root, winslash = "/", mustWork = TRUE), "/")
  if (startsWith(path, root)) substring(path, nchar(root) + 1L) else basename(path)
}

essentiality_table_row_count <- function(path) {
  if (!nzchar(path %||% "") || !file.exists(path)) return(NA_integer_)
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) return(0L)
  header <- if (startsWith(trimws(lines[[1]]), "RDF")) 2L else 1L
  max(0L, length(lines) - header)
}

essentiality_server <- function(input, output, session, active, active_project_path,
                                repo_root, jobs, status_tick, go_to, python_bin,
                                log_text) {
  parent_samples <- reactiveVal(character())
  selected_target_value <- reactiveVal(NA_real_)
  selected_target_tag <- reactiveVal("")
  target_mode <- reactiveVal("recommended")
  final_target <- reactiveVal("")
  initialized_project <- reactiveVal("")
  selected_result_tag <- reactiveVal("")
  completed_seen <- reactiveVal("")
  pending_final_target <- reactiveVal("")

  set_selected_target <- function(tag, mode = "existing") {
    if (!essentiality_valid_tag(tag)) return(FALSE)
    selected_target_tag(tag)
    selected_target_value(essentiality_tag_value(tag))
    target_mode(mode)
    TRUE
  }

  observeEvent(active_project_path(), {
    path <- active_project_path()
    if (identical(path, initialized_project())) return()
    detected <- detect_essentiality_parent_samples(active())
    parent_samples(detected)
    recommendation <- essentiality_recommendation(active()$project_root)
    available <- discover_essentiality_targets(active()$project_root)
    initial <- if (recommendation$available) recommendation$tag else if (nrow(available)) available$Tag[[nrow(available)]] else ""
    if (nzchar(initial)) set_selected_target(initial, if (recommendation$available) "recommended" else "existing")
    else {
      selected_target_tag("")
      selected_target_value(NA_real_)
      target_mode("recommended")
    }
    final_target(essentiality_final_target(active()$project_root))
    selected_result_tag("")
    initialized_project(path)
  }, ignoreInit = FALSE)

  job_progress_server("normalization_job", jobs, active, compact = FALSE,
                      expected_stage = "sample_normalization")
  job_progress_server("summary_normalized_job", jobs, active, compact = FALSE,
                      expected_stage = "summary_normalized")
  job_progress_server("combine_job", jobs, active, compact = FALSE,
                      expected_stage = "combined_hits")
  job_progress_server("combined_summary_job", jobs, active, compact = FALSE,
                      expected_stage = "summary_combined")
  job_progress_server("classifier_job", jobs, active, compact = FALSE,
                      expected_stage = "classifier")

  observe({
    if (!jobs$job_is_running()) return()
    invalidateLater(1000, session)
    jobs$poll_job()
    progress <- jobs$last_progress()
    if (is.null(progress) || progress$status %in% c("queued", "starting", "running")) return()
    key <- paste(progress$job_id %||% "", progress$status %||% "", sep = ":")
    if (!identical(key, completed_seen())) {
      completed_seen(key)
      status_tick(status_tick() + 1L)
      final_target(essentiality_final_target(active()$project_root))
    }
  })

  available_targets <- reactive({
    status_tick()
    discover_essentiality_targets(active()$project_root)
  })
  recommendation <- reactive({
    status_tick()
    essentiality_recommendation(active()$project_root)
  })
  target_summary <- reactive({
    status_tick()
    essentiality_target_summary_data(active()$project_root)
  })
  target_evaluation <- reactive({
    status_tick()
    essentiality_target_evaluation_data(active()$project_root)
  })
  current_target <- reactive(list(value = selected_target_value(), tag = selected_target_tag()))
  stage_states <- reactive({
    progress <- jobs$current_progress()
    essentiality_stage_state(active(), parent_samples(), selected_target_tag(),
                             progress, repo_root)
  })

  output$essentiality_parent_summary <- renderUI({
    NULL
  })
  output$essentiality_parent_table <- DT::renderDT({
    review <- essentiality_parent_review(active(), parent_samples())
    if (!nrow(review)) return(NULL)
    DT::datatable(review, rownames = FALSE, selection = "none",
                  options = list(dom = "t", paging = FALSE, scrollX = TRUE))
  })
  output$essentiality_parent_override <- renderUI({
    inventory <- essentiality_parent_inventory(active())
    candidates <- inventory$Sample[inventory$Included & !inventory$Explicit_treated]
    tagList(
      tags$details(open = if (!length(parent_samples())) NA else NULL,
                   tags$summary("Override parent selection"),
                   selectizeInput("essentiality_parent_override_select",
                                  "Temporary parent libraries",
                                  choices = as.list(setNames(candidates, candidates)),
                                  selected = parent_samples(), multiple = TRUE),
                   helpText("This session-only selection does not alter sample_sheet.csv. Libraries explicitly marked treated are unavailable."))
    )
  })
  observeEvent(input$essentiality_parent_override_select, {
    selected <- detect_essentiality_parent_samples(active(),
                                                   input$essentiality_parent_override_select %||% character())
    parent_samples(selected)
  }, ignoreInit = TRUE)

  selected_target_notice <- function(target_tag = selected_target_tag(),
                                     target_value = selected_target_value(),
                                     include_action = TRUE) {
    if (!nzchar(target_tag %||% "")) {
      return(tags$div(
        class = "alert alert-warning",
        tags$p("Choose a normalization target before continuing."),
        if (isTRUE(include_action))
          actionButton("essentiality_go_normalize_choose",
                       "Go to Normalize & Choose Target", class = "btn-primary")
      ))
    }
    tagList(
      tags$p(tags$b(format_target_label(target_value, target_tag, "user"))),
      tags$details(
        tags$summary("Technical details"),
        tags$p(format_target_label(target_value, target_tag, "technical"))
      )
    )
  }
  observeEvent(input$essentiality_go_normalize_choose,
               go_to("essentiality", "essentiality_tabs", "normalize"),
               ignoreInit = TRUE)

  register_stage_job_disclosure <- function(prefix, expected_stage, module_id) {
    history_id <- paste0(prefix, "_job_history")
    output[[paste0(prefix, "_job_disclosure")]] <- renderUI({
      status_tick()
      progress <- jobs$current_progress()
      active_match <- essentiality_job_matches(progress, expected_stage) &&
        (progress$status %||% "") %in% c("queued", "starting", "running")
      history <- essentiality_stage_job_history(active()$project_root, expected_stage)
      history_ui <- if (nrow(history)) tags$details(
        tags$summary("View stage job history"),
        DT::DTOutput(history_id)
      ) else NULL
      if (active_match)
        return(tagList(job_progress_ui(module_id), history_ui))
      if (!nrow(history)) return(NULL)
      last <- history[1L, , drop = FALSE]
      tagList(
        tags$details(
          tags$summary("Last stage job"),
          tags$dl(
            class = "ytab-meta",
            tags$dt("Time"), tags$dd(last$Time),
            tags$dt("Stage"), tags$dd(last$Stage),
            tags$dt("Target"), tags$dd(last$Target),
            tags$dt("Mode"), tags$dd(last$Mode),
            tags$dt("Status"), tags$dd(last$Status),
            tags$dt("Elapsed"), tags$dd(last$Elapsed)
          )
        ),
        history_ui
      )
    })
    output[[history_id]] <- DT::renderDT({
      history <- essentiality_stage_job_history(active()$project_root, expected_stage)
      if (!nrow(history)) return(NULL)
      DT::datatable(
        history, rownames = FALSE, selection = "none",
        options = list(dom = "t", pageLength = 10, scrollX = TRUE)
      )
    })
  }
  register_stage_job_disclosure(
    "normalization", "sample_normalization", "normalization_job"
  )
  register_stage_job_disclosure(
    "summary_normalized", "summary_normalized", "summary_normalized_job"
  )
  register_stage_job_disclosure("combine", "combined_hits", "combine_job")
  register_stage_job_disclosure(
    "combined_summary", "summary_combined", "combined_summary_job"
  )
  register_stage_job_disclosure("classifier", "classifier", "classifier_job")

  manual_normalization <- reactive(parse_essentiality_targets(input$normalization_manual_targets %||% ""))
  output$normalization_target_preview <- DT::renderDT({
    preview <- manual_normalization()$table
    DT::datatable(preview, rownames = FALSE, selection = "none",
                  options = list(dom = "t", paging = FALSE))
  })
  normalization_target_argument <- reactive({
    if (identical(input$normalization_target_mode %||% "auto", "auto")) return("auto")
    parsed <- manual_normalization()
    validate(need(parsed$valid, "Enter at least one valid, unique positive target."))
    parsed$argument
  })
  output$normalization_prerun_summary <- renderUI({
    target <- tryCatch(normalization_target_argument(), error = function(e) "")
    tags$dl(class = "ytab-execution-summary",
            tags$dt("Parents"), tags$dd(length(parent_samples())),
            tags$dt("Target mode"), tags$dd(if (identical(target, "auto")) "Auto scan" else "Manual target list"),
            tags$dt("Minimum site retention"), tags$dd(input$normalization_min_retention %||% .95),
            tags$dt("Upsampling"), tags$dd("Never"),
            tags$dt("Output"), tags$dd(tags$code("sample_normalization/")))
  })
  normalization_cached <- reactive({
    if (!length(parent_samples())) return(FALSE)
    if (identical(input$normalization_target_mode %||% "auto", "auto")) {
      rec <- recommendation()
      return(rec$available && all(essentiality_normalized_inputs(active()$project_root,
                                                                 rec$tag, parent_samples())$Available))
    }
    parsed <- manual_normalization()
    parsed$valid && all(vapply(parsed$tags, function(tag)
      all(essentiality_normalized_inputs(active()$project_root, tag, parent_samples())$Available),
      logical(1)))
  })
  output$normalization_run_panel <- renderUI({
    if (normalization_cached()) return(NULL)
    content <- essentiality_normalization_configuration_ui()
    essentiality_run_panel_ui(FALSE, "", "Configure parent normalization", content)
  })
  output$normalization_rerun_panel <- renderUI({
    if (!normalization_cached()) return(NULL)
    tagList(tags$h4("Rerun or reconfigure"), essentiality_normalization_configuration_ui())
  })
  output$normalization_action <- renderUI({
    mode <- input$normalization_execution_mode %||% "preview"
    force <- FALSE
    valid <- length(parent_samples()) > 0L &&
      (identical(input$normalization_target_mode %||% "auto", "auto") || manual_normalization()$valid)
    actionButton("run_sample_normalization",
                 essentiality_action_label(mode, normalization_cached(), force,
                                           "Preview normalization scan", "Run normalization scan",
                                           "Use cached normalization", "Rerun normalization"),
                 class = "btn-primary",
                 disabled = if (!valid || jobs$job_is_running()) "disabled" else NULL)
  })

  evaluate_target_argument <- reactive({
    source <- input$summary_normalized_target_source %||% "all"
    if (source == "all") return("all")
    if (source == "recommended") {
      rec <- recommendation()
      validate(need(rec$available, "No recommended target is available."))
      return(rec$tag)
    }
    if (source == "existing") {
      tags <- input$summary_normalized_existing_targets %||% character()
      validate(need(length(tags), "Select at least one existing target."))
      return(paste(tags, collapse = ","))
    }
    parsed <- parse_essentiality_targets(input$summary_normalized_manual_targets %||% "")
    validate(need(parsed$valid, "Enter at least one valid positive target."))
    parsed$argument
  })
  observe({
    choices <- essentiality_target_choices(available_targets())
    selected <- intersect(input$summary_normalized_existing_targets %||% character(),
                          unname(unlist(choices, use.names = FALSE)))
    updateSelectizeInput(session, "summary_normalized_existing_targets",
                         choices = choices, selected = selected, server = TRUE)
  })
  output$summary_normalized_prerun_summary <- renderUI({
    argument <- tryCatch(evaluate_target_argument(), error = function(e) "")
    tags$dl(class = "ytab-execution-summary",
            tags$dt("Targets"), tags$dd(if (nzchar(argument)) argument else "Not selected"),
            tags$dt("Parent libraries"), tags$dd(length(parent_samples())),
            tags$dt("Expected SummaryTable jobs"),
            tags$dd({
              count <- if (argument == "all") nrow(available_targets()) else length(strsplit(argument, ",", fixed = TRUE)[[1]])
              count * length(parent_samples())
            }),
            tags$dt("Site-retention threshold"),
            tags$dd(recommendation()$data$retention_threshold %||% .95),
            tags$dt("Feature-retention threshold"),
            tags$dd(input$summary_normalized_min_retention %||% .95))
  })
  output$summary_normalized_action <- renderUI({
    mode <- input$summary_normalized_execution_mode %||% "preview"
    force <- FALSE
    valid <- length(parent_samples()) > 0L &&
      !inherits(try(evaluate_target_argument(), silent = TRUE), "try-error")
    cached <- nrow(target_summary()) > 0L
    actionButton("run_summary_normalized",
                 essentiality_action_label(mode, cached, force,
                                           "Preview feature-level evaluation",
                                           "Run feature-level evaluation",
                                           "Use cached target evaluation",
                                           "Rerun feature-level evaluation"),
                 class = "btn-primary",
                 disabled = if (!valid || jobs$job_is_running()) "disabled" else NULL)
  })
  output$summary_normalized_run_panel <- renderUI({
    complete <- identical(essentiality_feature_evaluation_state(active()$project_root),
                          "feature_evaluation_complete")
    if (complete) return(NULL)
    content <- essentiality_evaluation_configuration_ui()
    essentiality_run_panel_ui(FALSE, "", "Run feature-level target evaluation", content)
  })
  output$summary_normalized_rerun_panel <- renderUI({
    if (!identical(essentiality_feature_evaluation_state(active()$project_root),
                   "feature_evaluation_complete")) return(NULL)
    tagList(tags$h4("Rerun or reconfigure"), essentiality_evaluation_configuration_ui())
  })

  launch_essentiality_stage <- function(stage, target, mode, force, selected_items,
                                        min_site = .95, min_feature = .95,
                                        auto_min = NA_real_, auto_max = NA_real_,
                                        auto_step = 5, seed = 0L,
                                        keep_going = FALSE) {
    if (jobs$job_is_running()) {
      showNotification(paste("A heavy job is already active:", jobs$stage()),
                       type = "warning")
      return(FALSE)
    }
    plan <- build_essentiality_command(
      repo_root, active_project_path(), stage, target = target,
      parents = parent_samples(), execution_mode = mode, force = force,
      threads = active()$threads, min_site_retention = min_site,
      min_feature_retention = min_feature, auto_min = auto_min,
      auto_max = auto_max, auto_step = auto_step, seed = seed,
      keep_going = keep_going
    )
    metadata <- list(
      execution_mode = mode, dry_run = identical(mode, "preview"),
      force = identical(mode, "run") && isTRUE(force),
      target = target, target_tag = selected_target_tag(),
      parent_samples = unname(parent_samples()),
      parent_count = length(parent_samples())
    )
    jobs$start_job(python_bin(), plan$full, active_project_path(), stage,
                   wd = repo_root, selected_items = selected_items,
                   metadata = metadata)
    TRUE
  }

  observeEvent(input$run_sample_normalization, {
    target <- normalization_target_argument()
    items <- if (identical(target, "auto")) "auto recommendation" else
      parse_essentiality_targets(target)$tags
    launch_essentiality_stage(
      "sample_normalization", target,
      input$normalization_execution_mode %||% "preview",
      FALSE, items,
      min_site = input$normalization_min_retention %||% .95,
      auto_min = input$normalization_auto_min,
      auto_max = input$normalization_auto_max,
      auto_step = input$normalization_auto_step %||% 5,
      keep_going = TRUE
    )
  }, ignoreInit = TRUE)

  observeEvent(input$run_summary_normalized, {
    target <- evaluate_target_argument()
    tags <- if (identical(target, "all")) available_targets()$Tag else strsplit(target, ",", fixed = TRUE)[[1]]
    items <- as.vector(outer(tags, parent_samples(), paste, sep = " × "))
    launch_essentiality_stage(
      "summary_normalized", target,
      input$summary_normalized_execution_mode %||% "preview",
      FALSE, items,
      min_feature = input$summary_normalized_min_retention %||% .95,
      keep_going = TRUE
    )
  }, ignoreInit = TRUE)

  output$normalization_result_summary <- renderUI({
    NULL
  })

  output$target_evaluation_empty <- renderUI({
    if (!nrow(target_summary())) tags$p(class = "text-muted",
      "Run at least one normalization target before feature-level evaluation.")
  })
  output$target_evaluation_status_summary <- renderUI({
    NULL
  })
  output$target_comparison_ui <- renderUI({
    data <- target_summary()
    visual <- essentiality_target_visual_state(data, recommendation())
    if (identical(visual$mode, "empty"))
      return(tags$div(class = "ytab-empty-state", visual$message))
    if (identical(visual$mode, "single")) {
      view <- essentiality_target_summary_view(data)
      row <- view[1L, , drop = FALSE]
      cells <- lapply(names(row), function(name)
        tags$div(tags$b(as.character(row[[name]][[1L]])), name))
      return(tagList(
        tags$div(class = "ytab-single-target-card", cells),
        tags$p(class = "text-muted", visual$message)
      ))
    }
    if (identical(visual$mode, "compact"))
      return(tagList(
        DT::DTOutput("normalization_target_summary_table"),
        selectInput("normalization_target_plot_choice", "Plot",
                    choices = c("Combined retention" = "combined",
                                "Site retention" = "site",
                                "Feature retention" = "feature")),
        uiOutput("normalization_target_selected_plot")
      ))
    tagList(
      DT::DTOutput("normalization_target_summary_table"),
      selectInput("normalization_target_plot_choice", "Plot",
                  choices = c("Combined retention" = "combined",
                              "Site retention" = "site",
                              "Feature retention" = "feature")),
      uiOutput("normalization_target_selected_plot")
    )
  })
  output$normalization_target_summary_table <- DT::renderDT({
    view <- essentiality_target_summary_view(target_summary())
    if (nrow(view) < 2L) return(NULL)
    DT::datatable(view, rownames = FALSE, filter = "top", selection = "none",
                  options = list(scrollX = TRUE, pageLength = 10))
  })
  retention_plot <- function(kind = c("site", "feature", "combined")) {
    kind <- match.arg(kind)
    data <- target_summary()
    if (!nrow(data) || !"target" %in% names(data)) {
      qc_plot_empty("No target-retention values are available.")
      return()
    }
    old <- qc_plot_par(mar = c(5, 5, 3, 1))
    on.exit(par(old), add = TRUE)
    target <- suppressWarnings(as.numeric(data$target))
    site_col <- intersect(c("min_hit_site_retention_fraction",
                            "mean_hit_site_retention_fraction"), names(data))
    feature_col <- intersect(c("min_feature_retention_fraction",
                               "mean_feature_retention_fraction"), names(data))
    site <- if (length(site_col)) suppressWarnings(as.numeric(data[[site_col[[1]]]])) else rep(NA, nrow(data))
    feature <- if (length(feature_col)) suppressWarnings(as.numeric(data[[feature_col[[1]]]])) else rep(NA, nrow(data))
    if (kind == "site") {
      plot(target, site, type = "b", pch = 16, col = "black", lwd = qc_plot_lwd,
           xlab = "Target", ylab = "Site retention", main = "Site retention by target",
           ylim = range(c(site, 0, 1), na.rm = TRUE))
    } else if (kind == "feature") {
      plot(target, feature, type = "b", pch = 16, col = "black", lwd = qc_plot_lwd,
           xlab = "Target", ylab = "Feature retention",
           main = "Feature retention by target",
           ylim = range(c(feature, 0, 1), na.rm = TRUE))
    } else {
      plot(target, site, type = "b", pch = 16, col = "black", lwd = qc_plot_lwd,
           xlab = "Target", ylab = "Retention", main = "Combined target-retention view",
           ylim = range(c(site, feature, 0, 1), na.rm = TRUE))
      lines(target, feature, type = "b", pch = 17, col = "black", lwd = qc_plot_lwd)
      legend("bottomleft", c("Site", "Feature"), col = c("black", "black"),
             pch = c(16, 17), lty = 1, lwd = qc_plot_lwd, bty = "n", text.font = 2)
    }
    rec <- normalize_recommendation_state(recommendation())
    if (rec$available) {
      index <- which.min(abs(target - rec$value))
      if (length(index) && is.finite(target[[index]]))
        text(target[[index]], if (kind == "feature") feature[[index]] else site[[index]],
             "recommended", pos = 3, cex = .8)
    }
  }
  output$normalization_site_retention_plot <- renderPlot(retention_plot("site"))
  output$normalization_feature_retention_plot <- renderPlot(retention_plot("feature"))
  output$normalization_combined_retention_plot <- renderPlot(retention_plot("combined"))
  output$normalization_target_selected_plot <- renderUI({
    choice <- input$normalization_target_plot_choice %||% "combined"
    output_id <- switch(choice, site = "normalization_site_retention_plot",
                        feature = "normalization_feature_retention_plot",
                        "normalization_combined_retention_plot")
    plotOutput(output_id, height = "320px")
  })

  output$normalization_recommendation_card <- renderUI({
    rec <- normalize_recommendation_state(recommendation())
    if (!rec$available) return(tags$div(class = "ytab-result-card ytab-result-empty",
      tags$h4("Recommended target"), tags$p("No recommendation is available.")))
    summary <- target_summary()
    row <- if (nrow(summary) && "target_tag" %in% names(summary))
      summary[summary$target_tag == rec$tag, , drop = FALSE] else data.frame()
    value <- function(column, fallback = "Not available")
      if (nrow(row) && column %in% names(row)) row[[column]][[1]] else fallback
    hover <- paste(
      paste("Recommendation type:", rec$label),
      paste("Site retention:", value("min_hit_site_retention_fraction",
                                     rec$site_retention)),
      paste("Feature retention:", value("min_feature_retention_fraction",
                                        rec$feature_retention)),
      paste("Parent libraries passing:", value("sample_count",
                                               rec$parents_passing)),
      paste("Reason:", rec$reason),
      sep = "\n"
    )
    tagList(
      tags$div(class = "ytab-result-card ytab-result-matching", title = hover,
               tags$div(class = "ytab-result-heading",
                        tags$h4("Recommended target"),
                        tags$span(class = "ytab-result-badge", rec$tag)),
               tags$p(format_target_label(rec$value, rec$tag, "user")),
               if (rec$preliminary) tags$p(class = "ytab-warning",
                 "Preliminary recommendation; confirm with feature-level evaluation."),
               actionButton("use_recommended_target", "Use recommended target",
                            class = "btn-secondary"),
               NULL)
    )
  })
  observeEvent(input$use_recommended_target, {
    rec <- recommendation()
    if (rec$available) set_selected_target(rec$tag, "recommended")
  }, ignoreInit = TRUE)
  observeEvent(input$continue_to_combine, {
    req(nzchar(selected_target_tag()))
    go_to("essentiality", "essentiality_tabs", "combine")
  }, ignoreInit = TRUE)

  selected_evaluated_passes <- reactive({
    tag <- input$essentiality_target_choice %||% ""
    data <- target_summary()
    if (!nzchar(tag) || !nrow(data) || !"target_tag" %in% names(data)) return(TRUE)
    row <- data[data$target_tag == tag, , drop = FALSE]
    if (!nrow(row)) return(TRUE)
    site <- suppressWarnings(as.numeric(row$min_hit_site_retention_fraction %||% NA))
    feature <- suppressWarnings(as.numeric(row$min_feature_retention_fraction %||% NA))
    site_threshold <- suppressWarnings(as.numeric(recommendation()$data$retention_threshold %||% .95))
    feature_threshold <- suppressWarnings(as.numeric(input$summary_normalized_min_retention %||% .95))
    is.finite(site) && site >= site_threshold && is.finite(feature) && feature >= feature_threshold
  })
  output$essentiality_target_selection <- renderUI({
    data <- target_summary()
    if (!nrow(data) || !"target_tag" %in% names(data)) return(NULL)
    tags <- unique(as.character(data$target_tag))
    values <- vapply(tags, essentiality_tag_value, numeric(1))
    choices <- as.list(setNames(tags, vapply(tags, function(tag)
      format_target_label(essentiality_tag_value(tag), tag, "user"), character(1))))
    tagList(
      tags$p(if (nzchar(selected_target_tag()))
        paste("Selected normalization target:",
              essentiality_numeric_text(selected_target_value())) else
        "No normalization target has been selected."),
      tags$details(
        tags$summary("Choose another target"),
        selectInput("essentiality_target_choice", "Available targets",
                    choices = choices,
                    selected = if (selected_target_tag() %in% tags)
                      selected_target_tag() else tags[[1]]),
        if (!selected_evaluated_passes())
          tagList(tags$div(class = "alert alert-warning",
                           "This target did not meet the configured retention criteria and may remove insertion or feature information."),
                  checkboxInput("essentiality_failed_target_ack",
                                "I acknowledge the retention warning.", FALSE)),
        actionButton("use_selected_evaluated_target", "Use selected target",
                     disabled = if (!selected_evaluated_passes() &&
                                   !isTRUE(input$essentiality_failed_target_ack)) "disabled" else NULL)
      ),
      if (nzchar(selected_target_tag()))
        actionButton("continue_to_combine", "Continue to Combine Parents",
                     class = "btn-primary")
    )
  })
  observeEvent(input$use_selected_evaluated_target, {
    if (!selected_evaluated_passes() && !isTRUE(input$essentiality_failed_target_ack)) {
      showNotification("Acknowledge the retention warning before selecting this target.",
                       type = "warning")
      return()
    }
    if (set_selected_target(input$essentiality_target_choice %||% "", "existing"))
      showNotification("Essentiality target updated.")
  }, ignoreInit = TRUE)

  combine_inputs <- reactive({
    if (!nzchar(selected_target_tag()) || !length(parent_samples())) return(data.frame())
    essentiality_normalized_inputs(active()$project_root, selected_target_tag(),
                                   parent_samples())
  })
  output$combine_prerequisite_summary <- renderUI({
    if (!nzchar(selected_target_tag()))
      return(selected_target_notice())
    inputs <- combine_inputs()
    missing <- if (nrow(inputs)) inputs$Sample[!inputs$Available] else parent_samples()
    tagList(
      selected_target_notice(include_action = FALSE),
      tags$dl(class = "ytab-execution-summary",
              tags$dt("Parent libraries"), tags$dd(length(parent_samples())),
              tags$dt("Normalized hit files available"),
              tags$dd(if (nrow(inputs)) sum(inputs$Available) else 0L),
              tags$dt("Missing inputs"),
              tags$dd(if (length(missing)) paste(missing, collapse = ", ") else "None")),
      if (length(missing)) tags$div(class = "alert alert-warning",
        "Select a target with normalized outputs for every parent library.",
        actionButton("combine_back_to_normalize",
                     "Go to Normalize & Choose Target"))
    )
  })
  output$combine_parent_inputs <- DT::renderDT({
    inputs <- combine_inputs()
    if (!nrow(inputs)) return(NULL)
    visible <- data.frame(Sample = inputs$Sample,
                          `Normalized hit file` = vapply(inputs$Path, essentiality_relative_path,
                                                         character(1), active()$project_root),
                          Available = ifelse(inputs$Available, "Yes", "No"),
                          check.names = FALSE)
    DT::datatable(visible, rownames = FALSE, selection = "none",
                  options = list(dom = "t", paging = FALSE, scrollX = TRUE))
  })
  observeEvent(input$combine_back_to_normalize,
               go_to("essentiality", "essentiality_tabs", "normalize"),
               ignoreInit = TRUE)
  combined_hits_path <- reactive(if (nzchar(selected_target_tag()))
    file.path(active()$project_root, "combined_hits", selected_target_tag(),
              paste0("combined_parent_hits.", selected_target_tag(), ".txt")) else "")
  combined_hits_cached <- reactive(nzchar(combined_hits_path()) &&
    file.exists(combined_hits_path()) && file.info(combined_hits_path())$size > 0)
  output$combine_run_panel <- renderUI({
    if (combined_hits_cached()) return(NULL)
    content <- tagList(
      uiOutput("combine_prerun_configuration"),
      essentiality_execution_controls("combine")
    )
    essentiality_run_panel_ui(FALSE, "", "Configure parent combination", content)
  })
  output$combine_rerun_panel <- renderUI({
    if (!combined_hits_cached()) return(NULL)
    tagList(tags$h4("Recombine parent libraries"),uiOutput("combine_prerun_configuration"),essentiality_execution_controls("combine"))
  })
  output$combine_prerun_configuration <- renderUI({
    inputs <- combine_inputs()
    tags$div(
      class = "ytab-config-summary",
      tags$span(tags$b(if (nzchar(selected_target_tag()))
        format_target_label(selected_target_value(), selected_target_tag(), "user")
        else "Not selected")),
      tags$span(tags$b(if (nrow(inputs)) sum(inputs$Available) else 0L),
                " parent inputs available")
    )
  })
  output$combine_action <- renderUI({
    mode <- input$combine_execution_mode %||% "preview"
    force <- FALSE
    ready <- nrow(combine_inputs()) > 0L && all(combine_inputs()$Available)
    actionButton("run_combine_hits",
                 essentiality_action_label(mode, combined_hits_cached(), force,
                                           "Preview parent combination",
                                           "Combine parent libraries",
                                           "Use cached combined library",
                                           "Recombine parent libraries"),
                 class = "btn-primary",
                 disabled = if (!ready || jobs$job_is_running()) "disabled" else NULL)
  })
  observeEvent(input$run_combine_hits, {
    req(nzchar(selected_target_tag()), nrow(combine_inputs()), all(combine_inputs()$Available))
    launch_essentiality_stage("combined_hits", selected_target_tag(),
                              input$combine_execution_mode %||% "preview",
                              FALSE, selected_target_tag())
  }, ignoreInit = TRUE)
  output$combined_hits_status <- renderUI({
    status_tick()
    path <- combined_hits_path()
    if (!nzchar(path) || !file.exists(path)) return(NULL)
    tags$div(
      class = "ytab-result-card ytab-result-matching",
      tags$div(
        class = "ytab-result-heading",
        tags$h4("Combined parent library complete"),
        tags$span(class = "ytab-result-badge", "Ready")
      ),
      tags$p(tags$b("Normalization target: "),
             format_target_label(selected_target_value(), selected_target_tag(), "user"))
    )
  })
  output$combined_hits_result <- renderUI({
    status_tick()
    path <- combined_hits_path()
    if (!nzchar(path) || !file.exists(path)) return(tags$p(class = "text-muted",
      "Select a target with normalized outputs for every parent library."))
    manifest <- essentiality_manifest(active()$project_root, "combine",
                                      selected_target_tag()) %||% list()
    tagList(
      tags$div(
        class = "ytab-result-heading",
        tags$h4("Combined parent library complete"),
        tags$span(class = "ytab-result-badge",
                  if (identical(manifest$status %||% "", "skipped"))
                    "Cached result" else "Complete")
      ),
      tags$dl(class = "ytab-meta",
            tags$dt("Normalization target"),
            tags$dd(format_target_label(selected_target_value(),
                                        selected_target_tag(), "user")),
            tags$dt("Parent libraries combined"),
            tags$dd(length(manifest$selected_samples %||% parent_samples())),
            tags$dt("Combined insertion sites"),
            tags$dd(manifest$total_combined_sites %||% "Not recorded"),
            tags$dt("Combined reads"),
            tags$dd(manifest$total_combined_reads %||% "Not recorded"),
            tags$dt("Combined hit-file path"),
            tags$dd(tags$code(essentiality_relative_path(path, active()$project_root))),
            tags$dt("Manifest"),
            tags$dd(tags$code(file.path("manifests", "combined_hits",
                                        paste0(selected_target_tag(), ".combined_hits_manifest.json")))))
    )
  })
  observeEvent(input$combine_continue_summary,
               go_to("essentiality", "essentiality_tabs", "combined_summary"),
               ignoreInit = TRUE)

  combined_feature_path <- reactive(if (nzchar(selected_target_tag()))
    file.path(active()$project_root, "summary_combined", selected_target_tag(),
              paste0("combined_feature_table.", selected_target_tag(), ".txt")) else "")
  combined_feature_cached <- reactive(nzchar(combined_feature_path()) &&
    file.exists(combined_feature_path()) && file.info(combined_feature_path())$size > 0)
  output$combined_summary_prerequisite <- renderUI({
    if (!nzchar(selected_target_tag()))
      return(selected_target_notice())
    path <- combined_hits_path()
    reference <- active()$reference %||% list()
    annotation <- c(reference$feature_table, reference$gff, reference$gtf)
    annotation <- annotation[nzchar(as.character(annotation))]
    annotation_ready <- length(annotation) && any(file.exists(vapply(annotation,
      ytab_resolve_path, character(1), repo_root)))
    tagList(
      selected_target_notice(include_action = FALSE),
      tags$dl(class = "ytab-execution-summary",
              tags$dt("Combined hit file"),
              tags$dd(if (nzchar(path) && file.exists(path))
                tags$code(essentiality_relative_path(path, active()$project_root)) else "Missing"),
              tags$dt("Annotation/reference status"),
              tags$dd(if (annotation_ready) "Available" else "Missing"),
              tags$dt("Expected feature-table output"),
              tags$dd(if (nzchar(selected_target_tag()))
                tags$code(file.path("summary_combined", selected_target_tag(),
                                    paste0("combined_feature_table.",
                                           selected_target_tag(), ".txt"))) else "Target required")),
      if (!combined_hits_cached()) tags$div(class = "alert alert-warning",
        "Create the combined normalized parent hit file first.")
    )
  })
  output$combined_summary_action <- renderUI({
    mode <- input$combined_summary_execution_mode %||% "preview"
    force <- FALSE
    actionButton("run_summary_combined",
                 essentiality_action_label(mode, combined_feature_cached(), force,
                                           "Preview combined SummaryTable",
                                           "Build combined feature table",
                                           "Use cached feature table",
                                           "Rebuild combined feature table"),
                 class = "btn-primary",
                 disabled = if (!combined_hits_cached() || jobs$job_is_running())
                   "disabled" else NULL)
  })
  output$combined_summary_run_panel <- renderUI({
    if (combined_feature_cached()) return(NULL)
    content <- essentiality_execution_controls("combined_summary")
    essentiality_run_panel_ui(FALSE, "", "Configure combined SummaryTable", content)
  })
  output$combined_summary_rerun_panel <- renderUI({
    if (!combined_feature_cached()) return(NULL)
    tagList(tags$h4("Rebuild combined feature table"),essentiality_execution_controls("combined_summary"))
  })
  observeEvent(input$run_summary_combined, {
    req(combined_hits_cached())
    launch_essentiality_stage("summary_combined", selected_target_tag(),
                              input$combined_summary_execution_mode %||% "preview",
                              FALSE,
                              selected_target_tag())
  }, ignoreInit = TRUE)
  output$combined_summary_status <- renderUI({
    status_tick()
    path <- combined_feature_path()
    if (!nzchar(path) || !file.exists(path)) return(NULL)
    tags$div(
      class = "ytab-result-card ytab-result-matching",
      tags$div(
        class = "ytab-result-heading",
        tags$h4("Combined parent feature table complete"),
        tags$span(class = "ytab-result-badge", "Ready")
      ),
      tags$p(tags$b("Features: "), essentiality_table_row_count(path))
    )
  })
  output$combined_summary_result <- renderUI({
    status_tick()
    path <- combined_feature_path()
    if (!nzchar(path) || !file.exists(path)) return(tags$p(class = "text-muted",
      "Create the combined normalized parent hit file first."))
    manifest <- essentiality_manifest(active()$project_root, "combined_summary",
                                      selected_target_tag()) %||% list()
    tagList(
      tags$div(
        class = "ytab-result-heading",
        tags$h4("Combined parent feature table complete"),
        tags$span(class = "ytab-result-badge", "Complete")
      ),
      tags$dl(class = "ytab-meta",
            tags$dt("Normalization target"),
            tags$dd(format_target_label(selected_target_value(),
                                        selected_target_tag(), "user")),
            tags$dt("Feature count"), tags$dd(essentiality_table_row_count(path)),
            tags$dt("Stable feature table"),
            tags$dd(tags$code(essentiality_relative_path(path, active()$project_root))),
            tags$dt("Summary statistics"),
            tags$dd(if (nzchar(manifest$summary_stats_file %||% ""))
              tags$code(essentiality_relative_path(manifest$summary_stats_file,
                                                   active()$project_root)) else "Not available"),
            tags$dt("Completion status"), tags$dd(manifest$status %||% "complete"))
    )
  })
  observeEvent(input$combined_summary_continue_classifier,
               go_to("essentiality", "essentiality_tabs", "classifier"),
               ignoreInit = TRUE)

  classifier_target <- reactive({
    tag <- selected_target_tag()
    validate(need(essentiality_valid_tag(tag),
                  "Choose a normalization target before continuing."))
    list(value = selected_target_value(), tag = tag,
         source = "Selected normalization target")
  })
  output$classifier_recommended_target <- renderUI({
    if (!nzchar(selected_target_tag()))
      return(selected_target_notice())
    feature <- file.path(
      active()$project_root, "summary_combined", selected_target_tag(),
      paste0("combined_feature_table.", selected_target_tag(), ".txt")
    )
    tagList(
      selected_target_notice(include_action = FALSE),
      tags$dl(
        class = "ytab-meta",
        tags$dt("Combined feature input"),
        tags$dd(if (file.exists(feature))
          tags$code(essentiality_relative_path(feature, active()$project_root)) else
          "Missing")
      )
    )
  })
  classifier_resources <- reactive({
    target <- tryCatch(classifier_target(), error = function(e) list(tag = ""))
    essentiality_classifier_resources(active(), repo_root, target$tag %||% "")
  })
  output$classifier_species_support <- renderUI({
    resources <- classifier_resources()
    if (!isTRUE(attr(resources, "supported_species")))
      tags$div(class = "alert alert-danger",
               paste("Classifier resources are currently unavailable for species",
                     shQuote(active()$species), ". The current classifier supports glabrata table mode."))
  })
  output$classifier_resource_table <- DT::renderDT({
    data <- classifier_resources()
    visible <- data[, c("Resource", "Status", "Requirement"), drop = FALSE]
    DT::datatable(visible, rownames = FALSE, selection = "none",
                  options = list(dom = "t", paging = FALSE))
  })
  output$classifier_readiness_summary <- renderUI({
    resources <- classifier_resources()
    available <- sum(resources$Status == "Available")
    required <- sum(resources$Requirement == "Required")
    missing <- resources$Resource[
      resources$Requirement == "Required" & resources$Status != "Available"
    ]
    tags$div(
      class = if (length(missing))
        "ytab-readiness-line is-blocked" else "ytab-readiness-line is-ready",
      tags$b(sprintf("Required resources: %d of %d available", available, required)),
      if (length(missing))
        tags$span(paste("Missing:", paste(missing, collapse = ", "))) else
        tags$span("Classifier inputs are ready.")
    )
  })
  classifier_seed <- reactive({
    text <- trimws(input$classifier_seed %||% "0")
    value <- suppressWarnings(as.integer(text))
    validate(need(nzchar(text) && !is.na(value) &&
                  identical(as.character(value), text),
                  "Random seed must be an integer."))
    value
  })
  classifier_prediction_path <- reactive({
    target <- tryCatch(classifier_target(), error = function(e) list(tag = ""))
    if (!nzchar(target$tag %||% "")) return("")
    file.path(active()$project_root, "classifier", target$tag,
              paste0("essentiality_predictions.", target$tag, ".csv"))
  })
  classifier_cached <- reactive({
    path <- classifier_prediction_path()
    target <- tryCatch(classifier_target(), error = function(e) list(tag = ""))
    seed <- tryCatch(classifier_seed(), error = function(e) NA_integer_)
    nzchar(path) && nzchar(target$tag %||% "") && !is.na(seed) &&
      essentiality_classifier_cache_matches(active()$project_root,
                                             target$tag, seed)
  })
  output$classifier_run_panel <- renderUI({
    if (classifier_cached()) return(NULL)
    content <- essentiality_classifier_configuration_ui()
    essentiality_run_panel_ui(FALSE, "", "Configure classifier", content)
  })
  output$classifier_rerun_panel <- renderUI({
    if (!classifier_cached()) return(NULL)
    tagList(tags$h4("Rerun classifier"),essentiality_classifier_configuration_ui())
  })
  classifier_ready <- reactive({
    target <- tryCatch(classifier_target(), error = function(e) NULL)
    resources <- classifier_resources()
    !is.null(target) && isTRUE(attr(resources, "supported_species")) &&
      all(resources$Status == "Available")
  })
  output$classifier_prerun_summary <- renderUI({
    target <- tryCatch(classifier_target(), error = function(e) list(value = NA,
                                                                    tag = ""))
    feature <- if (nzchar(target$tag)) file.path(active()$project_root,
                                                 "summary_combined", target$tag,
                                                 paste0("combined_feature_table.",
                                                        target$tag, ".txt")) else ""
    missing <- classifier_resources()$Resource[classifier_resources()$Status == "Missing"]
    tags$dl(class = "ytab-execution-summary",
            tags$dt("Normalization target"), tags$dd(if (nzchar(target$tag))
              format_target_label(target$value, target$tag, "user") else
              "Not resolved"),
            tags$dt("Combined feature table"),
            tags$dd(if (nzchar(feature) && file.exists(feature))
              tags$code(essentiality_relative_path(feature, active()$project_root)) else "Missing"),
            tags$dt("Feature count"),
            tags$dd(if (nzchar(feature)) essentiality_table_row_count(feature) else "Unavailable"),
            tags$dt("Random seed"),
            tags$dd(tryCatch(classifier_seed(), error = function(e) "Invalid")),
            tags$dt("Resource readiness"),
            tags$dd(if (length(missing)) paste("Missing:", paste(missing, collapse = ", ")) else "Ready"),
            tags$dt("Matching cache"), tags$dd(if (classifier_cached()) "Available" else "No"),
            tags$dt("Expected predictions file"),
            tags$dd(if (nzchar(target$tag)) tags$code(file.path("classifier", target$tag,
              paste0("essentiality_predictions.", target$tag, ".csv"))) else "Target required"))
  })
  output$classifier_action <- renderUI({
    mode <- input$classifier_execution_mode %||% "preview"
    force <- FALSE
    seed_valid <- !inherits(try(classifier_seed(), silent = TRUE), "try-error")
    actionButton("run_classifier",
                 essentiality_action_label(mode, classifier_cached(), force,
                                           "Preview classifier", "Run classifier",
                                           "Use cached predictions", "Rerun classifier"),
                 class = "btn-primary",
                 disabled = if (!classifier_ready() || !seed_valid ||
                               jobs$job_is_running()) "disabled" else NULL)
  })
  observeEvent(input$run_classifier, {
    req(classifier_ready())
    target <- classifier_target()
    set_selected_target(target$tag, "selected")
    launch_essentiality_stage("classifier", target$tag,
                              input$classifier_execution_mode %||% "preview",
                              FALSE, target$tag,
                              seed = classifier_seed())
  }, ignoreInit = TRUE)
  output$classifier_result_summary <- renderUI({
    status_tick()
    target <- tryCatch(classifier_target(), error = function(e) NULL)
    if (is.null(target)) return(tags$p(class = "text-muted",
      "Build a combined parent feature table before running the classifier."))
    path <- classifier_prediction_path()
    if (!file.exists(path)) return(tags$p(class = "text-muted",
      "No matching classifier result is available for this target."))
    result <- choose_essentiality_result(
      discover_classifier_results(active()$project_root, target$tag, final_target()),
      selected_target = target$tag, final_target = final_target())
    summary <- essentiality_result_summary(result)
    tagList(
      tags$div(
        class = "ytab-result-heading",
        tags$h4("Essentiality classifier complete"),
        tags$span(class = "ytab-result-badge",
                  if (classifier_cached()) "Matching result" else "Complete")
      ),
      tags$div(
        class = "ytab-stat-grid",
        tags$div(tags$b(format_target_label(target$value, target$tag, "badge")),
                 "Target"),
        tags$div(tags$b(summary$total), "Features"),
        lapply(names(summary$counts), function(label)
          tags$div(tags$b(unname(summary$counts[[label]])),
                   paste("Predicted", label))),
        tags$div(tags$b(summary$excluded), "Excluded")
      )
    )
  })

  classifier_results <- reactive({
    status_tick()
    discover_classifier_results(active()$project_root, selected_target_tag(),
                                final_target())
  })
  current_result <- reactive(choose_essentiality_result(
    classifier_results(), selected_result_tag(), selected_target_tag(),
    final_target()))
  results_view <- reactive(
    essentiality_results_view_contract(
      input$essentiality_results_view %||% "overview"
    )$view
  )
  output$essentiality_smoke_warning <- renderUI({
    if (essentiality_smoke_project(active()))
      tags$div(class = "alert alert-warning",
               "This project uses a shallow FASTQ subset for software validation. Essentiality predictions are not final biological conclusions.")
  })
  output$classifier_results_empty <- renderUI({
    if (!length(classifier_results())) tags$p(class = "text-muted",
      "No essentiality classifier result is available for this project.")
  })
  output$classifier_result_selector <- renderUI({
    results <- classifier_results()
    if (length(results) <= 1L) return(NULL)
    selectInput("classifier_result_choice", "Classifier result",
                choices = essentiality_result_labels(results,
                  recommendation()$tag, final_target()),
                selected = current_result()$tag %||% "")
  })
  observeEvent(input$classifier_result_choice,
               selected_result_tag(input$classifier_result_choice %||% ""),
               ignoreInit = TRUE)
  output$classifier_result_state <- renderUI({
    result <- current_result()
    if (is.null(result)) return(NULL)
    tags$div(class = if (grepl("matching|final", result$classification))
      "ytab-result-card ytab-result-matching" else "ytab-result-card ytab-result-historical",
      tags$div(class = "ytab-result-heading",
               tags$h4(tools::toTitleCase(result$classification)),
               tags$span(class = "ytab-result-badge", result$status)),
      tags$dl(class = "ytab-meta",
              tags$dt("Normalization target"),
              tags$dd(format_target_label(result$value, result$tag, "user")),
              tags$dt("Prediction rows"), tags$dd(result$rows),
              tags$dt("Selected normalization target"),
              tags$dd(if (identical(result$tag, selected_target_tag())) "Yes" else "No"),
              tags$dt("Saved final classifier target"),
              tags$dd(if (nzchar(final_target()))
                format_target_label(essentiality_tag_value(final_target()),
                                    final_target(), "badge") else "Not set")),
      tags$details(
        tags$summary("Technical details"),
        tags$p(format_target_label(result$value, result$tag, "technical"))
      ))
  })
  output$selected_target_results_summary <- renderUI({
    if (!nzchar(selected_target_tag()))
      return(selected_target_notice())
    selected_target_notice(include_action = FALSE)
  })
  output$essentiality_results_view_body <- renderUI({
    view <- results_view()
    if (identical(view, "overview")) {
      result <- current_result()
      columns <- if (is.null(result)) list(label = "") else
        essentiality_prediction_columns(read_essentiality_predictions(result$path))
      return(essentiality_results_overview_ui(nzchar(columns$label %||% "")))
    }
    if (identical(view, "predictions"))
      return(essentiality_results_predictions_ui())
    if (identical(view, "visualizations")) {
      result <- current_result()
      columns <- if (is.null(result)) list(label = "", score = "", metric = "") else
        essentiality_prediction_columns(read_essentiality_predictions(result$path))
      return(essentiality_results_visualizations_ui(columns))
    }
    if (identical(view, "provenance"))
      return(tagList(
        uiOutput("classifier_provenance"),
        essentiality_technical_details("classifier_results_technical")
      ))
    uiOutput("essentiality_download_cards")
  })
  output$essentiality_visualizations_body <- renderUI({
    result <- current_result()
    columns <- if (is.null(result)) list(label = "", score = "", metric = "") else
      essentiality_prediction_columns(read_essentiality_predictions(result$path))
    essentiality_results_visualizations_ui(columns)
  })
  output$essentiality_visualization_selector <- renderUI({
    result <- current_result()
    columns <- if (is.null(result)) list(label = "", score = "", metric = "") else
      essentiality_prediction_columns(read_essentiality_predictions(result$path))
    choices <- c()
    if (nzchar(columns$label %||% "")) choices <- c(choices, "Prediction-label distribution" = "label")
    if (nzchar(columns$score %||% "")) choices <- c(choices, "Score distribution" = "score")
    if (nzchar(columns$metric %||% "") && nzchar(columns$label %||% "")) choices <- c(choices, "Feature metric versus prediction" = "metric")
    plots <- essentiality_generated_plot_inventory(active())
    if (nrow(plots)) choices <- c(choices, setNames(paste0("generated:", seq_len(nrow(plots))),
                                                    paste("Generated:", plots$target, plots$title)))
    if (!length(choices)) return(tags$p(class = "text-muted", "No classifier visualizations are available."))
    selectInput("essentiality_visualization_choice", "Plot", choices = choices)
  })
  output$essentiality_selected_visualization <- renderUI({
    choice <- input$essentiality_visualization_choice %||% "label"
    if (startsWith(choice, "generated:")) {
      plots <- essentiality_generated_plot_inventory(active())
      index <- suppressWarnings(as.integer(sub("^generated:", "", choice)))
      if (is.na(index) || index < 1L || index > nrow(plots))
        return(tags$p(class = "text-muted", "Selected classifier plot is unavailable."))
      return(tags$article(
        class = "ytab-plot-card",
        tags$h4(paste(plots$target[[index]], plots$title[[index]])),
        tags$div(class = "ytab-diagnostic-preview",
          tags$img(src = plots$served_url[[index]], alt = plots$filename[[index]],
                   loading = "lazy", style = "width:100%;max-height:560px;object-fit:contain")
        ),
        tags$details(tags$summary("Show filename"), tags$code(plots$filename[[index]]))
      ))
    }
    output_id <- switch(choice, score = "classifier_score_plot",
                        metric = "classifier_metric_plot", "classifier_visual_label_plot")
    plotOutput(output_id, height = "360px")
  })
  result_data <- reactive({
    result <- current_result()
    if (is.null(result)) data.frame() else read_essentiality_predictions(result$path)
  })
  result_filtered <- reactive(filter_essentiality_results(
    result_data(), input$classifier_gene_search %||% "",
    input$classifier_label_filter %||% "All",
    input$classifier_inclusion_filter %||% "All"))
  observe({
    data <- result_data()
    columns <- essentiality_prediction_columns(data)
    labels <- if (nrow(data) && nzchar(columns$label))
      sort(unique(as.character(data[[columns$label]]))) else character()
    labels <- labels[!is.na(labels) & nzchar(labels)]
    choices <- as.list(c("All", labels))
    names(choices) <- c("All labels", labels)
    updateSelectInput(session, "classifier_label_filter", choices = choices,
                      selected = if ((input$classifier_label_filter %||% "All") %in%
                                    c("All", labels))
                        input$classifier_label_filter %||% "All" else "All")
  })
  output$classifier_summary_cards <- renderUI({
    result <- current_result()
    if (is.null(result)) return(NULL)
    summary <- essentiality_result_summary(result)
    label_counts <- paste(
      paste(names(summary$counts), unname(summary$counts), sep = ": "),
      collapse = "; "
    )
    tags$div(class = "ytab-stat-grid",
             title = paste("Label counts:", label_counts,
                           "\nExcluded features:", summary$excluded),
             tags$div(tags$b(summary$total), "Features"),
             tags$div(tags$b(format_target_label(result$value, result$tag, "badge")),
                      "Target"),
             tags$div(tags$b(if (identical(result$tag, final_target()))
               "Final target" else "Not final"), "Final target status"))
  })
  output$classifier_filtered_count <- renderUI({
    if (nrow(result_data())) tags$p(sprintf("Showing %d of %d features.",
                                            nrow(result_filtered()),
                                            nrow(result_data())))
  })
  output$classifier_predictions_table <- DT::renderDT({
    data <- essentiality_visible_results(result_filtered())
    if (!nrow(data)) return(NULL)
    DT::datatable(data, rownames = FALSE, filter = "none",
                  options = list(scrollX = TRUE, pageLength = 15,
                                 deferRender = TRUE, processing = TRUE))
  })
  classifier_label_plot <- function() {
    result <- current_result()
    req(!is.null(result))
    counts <- essentiality_result_summary(result)$counts
    if (!length(counts)) {
      qc_plot_empty("No classifier label column is available.")
    } else {
      old <- qc_plot_par(mar = c(5, 5, 3, 1))
      on.exit(par(old), add = TRUE)
      barplot(counts, col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
              ylab = "Features", main = "Prediction-label distribution")
    }
  }
  output$classifier_label_plot <- renderPlot({
    req(identical(results_view(), "overview"))
    classifier_label_plot()
  })
  output$classifier_visual_label_plot <- renderPlot({
    classifier_label_plot()
  })
  output$classifier_score_plot <- renderPlot({
    data <- result_data()
    columns <- essentiality_prediction_columns(data)
    if (!nrow(data) || !nzchar(columns$score)) {
      qc_plot_empty("No numeric prediction score is available.")
      return()
    }
    values <- suppressWarnings(as.numeric(data[[columns$score]]))
    values <- values[is.finite(values)]
    if (!length(values)) {
      qc_plot_empty("No numeric prediction score is available.")
    } else {
      old <- qc_plot_par(mar = c(5, 5, 3, 1))
      on.exit(par(old), add = TRUE)
      hist(values, col = qc_plot_fill, border = qc_plot_border, lwd = qc_plot_lwd,
           xlab = columns$score, main = "Score distribution")
    }
  })
  output$classifier_metric_plot <- renderPlot({
    data <- result_data()
    columns <- essentiality_prediction_columns(data)
    if (!nrow(data) || !nzchar(columns$metric) || !nzchar(columns$label)) {
      qc_plot_empty("Suitable feature-metric columns are unavailable.")
      return()
    }
    values <- suppressWarnings(as.numeric(data[[columns$metric]]))
    labels <- as.factor(data[[columns$label]])
    keep <- is.finite(values) & !is.na(labels)
    if (!any(keep)) {
      qc_plot_empty("Suitable feature-metric columns are unavailable.")
    } else {
      old <- qc_plot_par(mar = c(5, 5, 3, 1))
      on.exit(par(old), add = TRUE)
      boxplot(values[keep] ~ labels[keep], col = qc_plot_fill, border = qc_plot_border,
              lwd = qc_plot_lwd, xlab = "Classifier label", ylab = columns$metric,
              main = "Feature metric versus prediction")
    }
  })
  output$classifier_provenance <- renderUI({
    result <- current_result()
    if (is.null(result)) return(NULL)
    provenance <- essentiality_result_provenance(active()$project_root, result,
                                                 parent_samples())
    resource_count <- length(provenance$resources)
    tags$section(
      class = "ytab-review-section",
      tags$h4("Target provenance"),
      tags$dl(class = "ytab-meta",
              tags$dt("Normalization target"),
              tags$dd(format_target_label(provenance$target_value,
                                          provenance$target_tag, "user")),
              tags$dt("Recommendation type"), tags$dd(provenance$recommendation_type),
              tags$dt("Minimum site retention"), tags$dd(provenance$minimum_site_retention),
              tags$dt("Minimum feature retention"), tags$dd(provenance$minimum_feature_retention),
              tags$dt("Parent libraries"), tags$dd(provenance$parent_libraries),
              tags$dt("Combined sites"), tags$dd(provenance$combined_sites),
              tags$dt("Combined reads"), tags$dd(provenance$combined_reads),
              tags$dt("Combined feature count"), tags$dd(provenance$combined_feature_count),
              tags$dt("Random seed"), tags$dd(provenance$random_seed),
              tags$dt("Classifier resource hash groups"), tags$dd(resource_count),
              tags$dt("Final-target status"),
              tags$dd(if (provenance$final) "Final" else "Not final")),
      tags$details(
        tags$summary("Technical provenance"),
        tags$div(
          class = "ytab-technical-console",
          tags$p(format_target_label(provenance$target_value,
                                     provenance$target_tag, "technical")),
          tags$p(tags$b("Predictions: "), tags$code(result$path)),
          tags$p(tags$b("Classifier manifest: "),
                 tags$code(result$manifest_path %||% "Not available")),
          tags$p(tags$b("Input SHA-256: "),
                 tags$code(result$manifest$input_sha256 %||% "Not recorded")),
          tags$p(tags$b("Classifier SHA-256: "),
                 tags$code(result$manifest$classifier_script_sha256 %||%
                             "Not recorded")),
          tags$p(tags$b("Runner SHA-256: "),
                 tags$code(result$manifest$classifier_runner_sha256 %||%
                             "Not recorded"))
        )
      )
    )
  })
  output$classifier_final_target_action <- renderUI({
    result <- current_result()
    if (is.null(result) || result$status == "failed") return(NULL)
    tags$div(class = "ytab-result-card",
             tags$p(tags$b("Selected normalization target: "),
                    format_target_label(result$value, result$tag, "badge")),
             tags$p(tags$b("Saved final classifier target: "),
                    if (nzchar(final_target()))
                      format_target_label(essentiality_tag_value(final_target()),
                                          final_target(), "badge") else "Not set"),
             actionButton("set_final_classifier_target",
                          "Set this as final classifier target",
                          class = "btn-primary"))
  })
  show_final_target_modal <- function(result) {
    if (is.null(result)) return(invisible(FALSE))
    pending_final_target(result$tag)
    selected_result_tag(result$tag)
    showModal(modalDialog(
      title = paste("Set", result$tag, "as the final classifier target?"),
      tags$p("This records the reviewed target for project exports and downstream annotation. It does not rerun the classifier."),
      footer = tagList(modalButton("Cancel"),
                       actionButton("confirm_final_classifier_target",
                                    "Set final target", class = "btn-primary"))
    ))
    invisible(TRUE)
  }
  observeEvent(input$set_final_classifier_target, {
    result <- current_result()
    req(!is.null(result))
    show_final_target_modal(result)
  }, ignoreInit = TRUE)
  observeEvent(input$classifier_view_predictions, {
    target <- tryCatch(classifier_target(), error = function(e) NULL)
    if (!is.null(target)) selected_result_tag(target$tag)
    go_to("essentiality", "essentiality_tabs", "results")
  }, ignoreInit = TRUE)
  observeEvent(input$classifier_set_final_target, {
    target <- tryCatch(classifier_target(), error = function(e) NULL)
    req(!is.null(target))
    result <- choose_essentiality_result(
      classifier_results(), selected = target$tag,
      selected_target = target$tag, final_target = final_target()
    )
    req(!is.null(result), identical(result$tag, target$tag))
    show_final_target_modal(result)
  }, ignoreInit = TRUE)
  observeEvent(input$confirm_final_classifier_target, {
    tag <- pending_final_target()
    result <- choose_essentiality_result(
      classifier_results(), selected = tag,
      selected_target = tag, final_target = final_target()
    )
    req(!is.null(result))
    removeModal()
    args <- c(file.path(repo_root, "scripts", "local", "ytab_run_classifier.py"),
              "--project-config", active_project_path(), "--target", result$tag,
              "--save-final-target-only")
    command_result <- tryCatch(run_process_sync(python_bin(), args, wd = repo_root),
                               error = function(e) list(status = 1L, stdout = "",
                                                        stderr = conditionMessage(e)))
    log_text(paste(c(format_command_for_display(python_bin(), args),
                     command_result$stdout, command_result$stderr,
                     paste("Exit status:", command_result$status)), collapse = "\n"))
    if (identical(as.integer(command_result$status), 0L)) {
      final_target(essentiality_final_target(active()$project_root))
      pending_final_target("")
      status_tick(status_tick() + 1L)
      showNotification(paste(result$tag, "is now the final classifier target."))
    } else {
      pending_final_target("")
      showNotification("The final classifier target could not be saved.",
                       type = "error", duration = NULL)
    }
  }, ignoreInit = TRUE)

  technical_ui <- function(stage, manifest = NULL, inputs = character(),
                           outputs = character(), warnings = character()) {
    progress <- jobs$last_progress()
    job <- jobs$current_job()
    matches <- identical(job$stage %||% "", stage)
    failed <- matches && identical(progress$status %||% "", "failed")
    tags$details(
      class = "ytab-technical-details", open = if (failed) NA else NULL,
      tags$summary("Technical details"),
      tags$div(class = "ytab-technical-console",
               tags$p(tags$b("Exact command: "),
                      tags$code(if (matches) jobs$command() else
                        paste(manifest$command_run %||% "No command recorded",
                              collapse = " "))),
               tags$p(tags$b("Progress file: "),
                      tags$code(if (matches) job$progress_file %||% "Not available" else "Not available")),
               tags$p(tags$b("Manifest: "),
                      tags$code(manifest$manifest_path %||% "Not available")),
               if (length(inputs)) tags$p(tags$b("Inputs: "),
                 tags$code(paste(inputs, collapse = "; "))),
               if (length(outputs)) tags$p(tags$b("Outputs: "),
                 tags$code(paste(outputs, collapse = "; "))),
               if (length(warnings)) tags$p(tags$b("Warnings: "),
                 paste(warnings, collapse = "; ")),
               tags$details(tags$summary("stdout"),
                            tags$pre(if (matches) paste(jobs$current_stdout(),
                                                      collapse = "\n") else "No active stdout.")),
               tags$details(tags$summary("stderr"),
                            tags$pre(if (matches) paste(jobs$current_stderr(),
                                                      collapse = "\n") else "No active stderr.")))
    )
  }
  output$normalization_technical <- renderUI({
    tag <- if (recommendation()$available) recommendation()$tag else selected_target_tag()
    manifest <- if (nzchar(tag)) essentiality_manifest(active()$project_root,
                                                       "normalize", tag) else list()
    manifest$manifest_path <- if (nzchar(tag)) file.path("manifests", "sample_normalization",
      paste0(tag, ".sample_normalization_manifest.json")) else ""
    technical_ui("sample_normalization", manifest,
                 manifest$input_hit_files %||% character(),
                 manifest$normalized_hit_files %||% character(),
                 manifest$warnings %||% character())
  })
  output$summary_normalized_technical <- renderUI({
    technical_ui("summary_normalized", list(),
                 if (nrow(target_evaluation())) target_evaluation()$normalized_hits_file %||% character() else character(),
                 c(file.path("sample_normalization", "normalization_target_evaluation.csv"),
                   file.path("sample_normalization", "normalization_target_summary.csv")))
  })
  output$combine_technical <- renderUI({
    manifest <- essentiality_manifest(active()$project_root, "combine",
                                      selected_target_tag()) %||% list()
    manifest$manifest_path <- if (nzchar(selected_target_tag())) file.path("manifests",
      "combined_hits", paste0(selected_target_tag(),
                             ".combined_hits_manifest.json")) else ""
    technical_ui("combined_hits", manifest,
                 manifest$input_normalized_hit_files %||% character(),
                 manifest$combined_hits_file %||% character(),
                 manifest$warnings %||% character())
  })
  output$combined_summary_technical <- renderUI({
    manifest <- essentiality_manifest(active()$project_root,
                                      "combined_summary",
                                      selected_target_tag()) %||% list()
    manifest$manifest_path <- if (nzchar(selected_target_tag())) file.path("manifests",
      "summary_combined", paste0(selected_target_tag(),
                                ".summary_combined_manifest.json")) else ""
    technical_ui("summary_combined", manifest,
                 manifest$input_combined_hits_file %||% character(),
                 manifest$stable_combined_feature_table %||% character(),
                 manifest$warnings %||% character())
  })
  output$classifier_technical <- renderUI({
    target <- tryCatch(classifier_target(), error = function(e) list(tag = ""))
    manifest <- if (nzchar(target$tag)) essentiality_manifest(active()$project_root,
      "classifier", target$tag) %||% list() else list()
    manifest$manifest_path <- if (nzchar(target$tag)) file.path("manifests",
      "classifier", paste0(target$tag, ".classifier_manifest.json")) else ""
    technical_ui("classifier", manifest,
                 manifest$input_combined_feature_table %||% character(),
                 manifest$stable_prediction_table %||% character(),
                 manifest$warnings %||% character())
  })
  output$classifier_results_technical <- renderUI({
    result <- current_result()
    if (is.null(result)) return(tags$details(class = "ytab-technical-details",
      tags$summary("Technical details"), tags$p("No classifier result is available.")))
    technical_ui("classifier", result$manifest %||% list(),
                 result$manifest$input_combined_feature_table %||% character(),
                 result$path, result$manifest$warnings %||% character())
  })

  essentiality_download <- function(path, message) {
    if (!nzchar(path %||% "") || !file.exists(path) || dir.exists(path)) {
      showNotification(message, type = "warning")
      validate(need(FALSE, message))
    }
    path
  }
  current_downloads <- reactive({
    result <- current_result()
    if (is.null(result)) list() else essentiality_download_paths(active()$project_root,
                                                                 result)
  })
  output$essentiality_download_cards <- renderUI({
    req(identical(results_view(), "downloads"))
    paths <- current_downloads()
    availability <- essentiality_download_availability(paths)
    cards <- list(
      list("Predictions CSV", "download_essentiality_predictions", "predictions",
           "Stable classifier predictions"),
      list("Filtered predictions", "download_essentiality_filtered", "predictions",
           "Current prediction filters"),
      list("Combined parent feature table", "download_essentiality_feature_table",
           "combined_feature", "Classifier input table"),
      list("Target evaluation", "download_essentiality_evaluation",
           "target_evaluation", "Stored target-retention evaluation"),
      list("Recommendation record", "download_essentiality_recommendation",
           "normalization_recommendation", "Stored normalization recommendation"),
      list("Classifier manifest", "download_essentiality_manifest",
           "classifier_manifest", "Cache and resource provenance"),
      list("Final-target record", "download_essentiality_final_target",
           "final_target", "Explicitly saved project target")
    )
    tags$div(
      class = "ytab-download-grid",
      lapply(cards, function(card) {
        available <- isTRUE(unname(availability[card[[3L]]]))
        tags$section(
          class = paste("ytab-download-card", if (!available) "is-disabled" else ""),
          tags$h4(card[[1L]]),
          tags$p(card[[4L]]),
          if (available)
            downloadButton(card[[2L]], paste("Download", card[[1L]])) else
            tags$button(
              type = "button", class = "btn btn-default", disabled = "disabled",
              paste(card[[1L]], "unavailable")
            )
        )
      })
    )
  })
  output$download_essentiality_predictions <- downloadHandler(
    filename = function() basename(current_downloads()$predictions %||%
                                     "essentiality_predictions.csv"),
    content = function(file) file.copy(essentiality_download(
      current_downloads()$predictions %||% "",
      "No essentiality predictions are available."), file, overwrite = TRUE))
  output$download_essentiality_filtered <- downloadHandler(
    filename = function() paste0(current_result()$tag %||% "essentiality",
                                  ".filtered_predictions.csv"),
    content = function(file) write.csv(result_filtered(), file, row.names = FALSE))
  output$download_essentiality_feature_table <- downloadHandler(
    filename = function() basename(current_downloads()$combined_feature %||%
                                     "combined_feature_table.txt"),
    content = function(file) file.copy(essentiality_download(
      current_downloads()$combined_feature %||% "",
      "No combined parent feature table is available."), file, overwrite = TRUE))
  output$download_essentiality_evaluation <- downloadHandler(
    filename = function() "normalization_target_evaluation.csv",
    content = function(file) file.copy(essentiality_download(
      current_downloads()$target_evaluation %||% "",
      "No target-evaluation table is available."), file, overwrite = TRUE))
  output$download_essentiality_recommendation <- downloadHandler(
    filename = function() basename(current_downloads()$normalization_recommendation %||%
                                     "normalization_recommendation.json"),
    content = function(file) file.copy(essentiality_download(
      current_downloads()$normalization_recommendation %||% "",
      "No normalization recommendation is available."), file, overwrite = TRUE))
  output$download_essentiality_manifest <- downloadHandler(
    filename = function() paste0(current_result()$tag %||% "classifier",
                                  ".classifier_manifest.json"),
    content = function(file) file.copy(essentiality_download(
      current_downloads()$classifier_manifest %||% "",
      "No classifier manifest is available."), file, overwrite = TRUE))
  output$download_essentiality_final_target <- downloadHandler(
    filename = function() "final_classifier_target.txt",
    content = function(file) file.copy(essentiality_download(
      current_downloads()$final_target %||% "",
      "No final classifier target has been saved."), file, overwrite = TRUE))

  invisible(list(
    parent_samples = parent_samples,
    selected_target = selected_target_value,
    selected_target_tag = selected_target_tag,
    target_mode = target_mode,
    final_target = final_target,
    stage_states = stage_states,
    available_targets = available_targets,
    recommendation = recommendation
  ))
}
