navigation_state_defaults <- function(top = "overview") {
  list(
    active_top_tab = top,
    active_preprocessing_subtab = "samples",
    active_qc_subtab = "mapping_qc",
    active_essentiality_subtab = "normalize",
    active_fitness_subtab = "fitness_design",
    active_results_view = "overview"
  )
}

navigation_nested_fields <- c(
  preprocessing_tabs = "active_preprocessing_subtab",
  qc_tabs = "active_qc_subtab",
  essentiality_tabs = "active_essentiality_subtab",
  fitness_tabs = "active_fitness_subtab",
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
    preprocessing = "active_preprocessing_subtab",
    quality_control = "active_qc_subtab",
    essentiality = "active_essentiality_subtab",
    fitness = "active_fitness_subtab",
    ""
  )
  if (nzchar(field)) value[[field]] else NULL
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
    preprocessing = "preprocessing_tabs",
    quality_control = "qc_tabs",
    essentiality = "essentiality_tabs",
    fitness = "fitness_tabs",
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
