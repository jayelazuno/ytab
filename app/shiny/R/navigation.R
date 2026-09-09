navigation_state_defaults <- function(top = "overview") {
  list(
    active_top_tab = top,
    active_project_data_subtab = "project",
    active_preprocessing_subtab = "samples",
    active_qc_subtab = "mapping_qc",
    active_essentiality_subtab = "normalize",
    active_fitness_subtab = "fitness_design",
    active_go_analysis_subtab = "gene_set",
    active_comparative_subtab = "single_species",
    active_results_view = "overview"
  )
}

navigation_top_values <- c(
  "overview", "project_data", "pipeline", "preprocessing", "quality_control",
  "essentiality", "fitness", "go_analysis", "gene_domain_explorer", "comparative_species",
  "genome_browser", "results_exports", "logs_provenance", "help"
)

navigation_top_nested_ids <- c(
  project_data = "project_data_tabs",
  preprocessing = "preprocessing_tabs",
  quality_control = "qc_tabs",
  essentiality = "essentiality_tabs",
  fitness = "fitness_tabs",
  go_analysis = "go_analysis_tabs",
  comparative_species = "comparative_tabs"
)

navigation_nested_id_for_top <- function(top) {
  top <- as.character(top %||% "")
  if (!length(top) || !nzchar(top[[1L]])) return("")
  if (!top[[1L]] %in% names(navigation_top_nested_ids)) return("")
  unname(navigation_top_nested_ids[[top[[1L]]]])
}

navigation_nested_fields <- c(
  project_data_tabs = "active_project_data_subtab",
  preprocessing_tabs = "active_preprocessing_subtab",
  qc_tabs = "active_qc_subtab",
  essentiality_tabs = "active_essentiality_subtab",
  fitness_tabs = "active_fitness_subtab",
  go_analysis_tabs = "active_go_analysis_subtab",
  comparative_tabs = "active_comparative_subtab",
  essentiality_results_view = "active_results_view"
)

normalize_navigation_state <- function(value = NULL, top = NULL) {
  defaults <- navigation_state_defaults(top %||% "overview")
  if (!is.list(value)) return(defaults)
  for (name in intersect(names(defaults), names(value))) {
    candidate <- as.character(value[[name]] %||% "")
    if (length(candidate) && nzchar(candidate[[1L]])) defaults[[name]] <- candidate[[1L]]
  }
  if (identical(defaults$active_essentiality_subtab, "evaluate"))
    defaults$active_essentiality_subtab <- "normalize"
  if (!is.null(top) && nzchar(top)) defaults$active_top_tab <- top
  defaults
}

navigation_state_from_url <- function(value = NULL) {
  state <- normalize_navigation_state()
  if (!is.list(value)) return(state)
  top <- as.character(value$top %||% "")
  nested_id <- as.character(value$nested_id %||% "")
  nested_value <- as.character(value$nested_value %||% "")
  results_view <- as.character(value$results_view %||% "")
  if (length(top) && top[[1L]] %in% navigation_top_values)
    state <- record_navigation_selection(state, "workspace_tabs", top[[1L]])
  if (length(nested_id) && nested_id[[1L]] %in% names(navigation_nested_fields) &&
      length(nested_value) && nzchar(nested_value[[1L]]))
    state <- record_navigation_selection(state, nested_id[[1L]], nested_value[[1L]])
  if (length(results_view) && nzchar(results_view[[1L]]))
    state <- record_navigation_selection(state, "essentiality_results_view", results_view[[1L]])
  state
}

navigation_target_from_hash <- function(hash = "") {
  hash <- as.character(hash %||% "")
  if (!length(hash) || !nzchar(hash[[1L]])) return(NULL)
  hash <- sub("^#", "", hash[[1L]])
  if (!nzchar(hash)) return(NULL)
  parts <- strsplit(hash, "&", fixed = TRUE)[[1L]]
  params <- list()
  for (part in parts) {
    pair <- strsplit(part, "=", fixed = TRUE)[[1L]]
    if (length(pair) < 2L) next
    key <- utils::URLdecode(pair[[1L]])
    value <- utils::URLdecode(paste(pair[-1L], collapse = "="))
    params[[key]] <- value
  }
  top <- params$tab %||% params$top %||% ""
  if (!nzchar(top)) return(NULL)
  nested_value <- params$subtab %||% params$nested_value %||% ""
  nested_id <- params$nested_id %||% navigation_nested_id_for_top(top)
  list(
    top = top,
    nested_id = nested_id,
    nested_value = nested_value,
    results_view = params$results_view %||% "",
    project_id = params$project %||% params$project_id %||% ""
  )
}

record_navigation_selection <- function(value, input_id, selected) {
  value <- normalize_navigation_state(value)
  selected <- as.character(selected %||% "")
  if (!length(selected) || !nzchar(selected[[1L]])) return(value)
  if (identical(input_id, "workspace_tabs")) {
    value$active_top_tab <- selected[[1L]]
  } else if (input_id %in% names(navigation_nested_fields)) {
    value[[navigation_nested_fields[[input_id]]]] <- selected[[1L]]
  }
  value
}

navigation_subtab <- function(value, top = NULL) {
  value <- normalize_navigation_state(value)
  top <- top %||% value$active_top_tab
  field <- switch(top,
    project_data = "active_project_data_subtab",
    preprocessing = "active_preprocessing_subtab",
    quality_control = "active_qc_subtab",
    essentiality = "active_essentiality_subtab",
    fitness = "active_fitness_subtab",
    go_analysis = "active_go_analysis_subtab",
    comparative_species = "active_comparative_subtab",
    ""
  )
  if (nzchar(field)) value[[field]] else NULL
}

navigation_state_client_payload <- function(value, replace = FALSE, project_id = "") {
  value <- normalize_navigation_state(value)
  top <- value$active_top_tab
  nested_id <- navigation_nested_id_for_top(top)
  nested_value <- navigation_subtab(value, top) %||% ""
  list(
    top = top,
    nested_id = nested_id,
    nested_value = nested_value,
    results_view = value$active_results_view,
    project_id = as.character(project_id %||% ""),
    replace = isTRUE(replace)
  )
}

sync_navigation_url <- function(session, value, replace = FALSE, project_id = "") {
  session$sendCustomMessage("ytab-nav-state",
                            navigation_state_client_payload(value, replace, project_id))
}

navigate_workspace <- function(session, value, top_tab, nested_id = NULL,
                               nested_value = NULL, results_view = NULL) {
  value <- record_navigation_selection(value, "workspace_tabs", top_tab)
  if (!is.null(nested_id) && !is.null(nested_value))
    value <- record_navigation_selection(value, nested_id, nested_value)
  if (!is.null(results_view))
    value <- record_navigation_selection(value, "essentiality_results_view", results_view)
  updateNavbarPage(session, "workspace_tabs", selected = value$active_top_tab)
  selected_nested <- navigation_subtab(value, value$active_top_tab)
  selected_id <- switch(value$active_top_tab,
    project_data = "project_data_tabs",
    preprocessing = "preprocessing_tabs",
    quality_control = "qc_tabs",
    essentiality = "essentiality_tabs",
    fitness = "fitness_tabs",
    go_analysis = "go_analysis_tabs",
    comparative_species = "comparative_tabs",
    NULL
  )
  session$onFlushed(function() {
    if (!is.null(selected_id) && nzchar(selected_nested %||% ""))
      bslib::nav_select(selected_id, selected = selected_nested, session = session)
    if (identical(value$active_top_tab, "essentiality") &&
        identical(value$active_essentiality_subtab, "results"))
      updateRadioButtons(session, "essentiality_results_view",
                         selected = value$active_results_view)
  }, once = TRUE)
  value
}

restore_navigation_state <- function(session, value) {
  value <- normalize_navigation_state(value)
  navigate_workspace(session, value, value$active_top_tab)
}

resolve_continue_destination <- function(project,repo_root,status,reviewed=TRUE) {
  if(!reviewed)return(list(top="preprocessing",nested="samples"))
  stage<-next_preprocessing_stage(project,repo_root,status)
  values<-c(Reference="reference",Mapping="mapping",`Hit Files`="hit_files",`Summary Tables`="summary_tables",Complete="status")
  list(top="preprocessing",nested=unname(values[[stage]]%||%"samples"))
}
