# Workspace modules: top-level tab value -> nested navset id, state field and default subtab.
# Modules without subtabs have no nested id.
navigation_modules <- list(
  project = list(nested_id = "project_tabs", field = "active_project_subtab", default = "samples"),
  run = list(nested_id = "run_tabs", field = "active_run_subtab", default = "samples"),
  quality_control = list(nested_id = "qc_tabs", field = "active_qc_subtab", default = "overview"),
  normalization = list(nested_id = "normalization_tabs", field = "active_normalization_subtab", default = "site"),
  fitness = list(nested_id = "fitness_tabs", field = "active_fitness_subtab", default = "fitness_design"),
  essentiality = list(nested_id = "essentiality_tabs", field = "active_essentiality_subtab", default = "run"),
  gene_explorer = list(nested_id = "gene_explorer_tabs", field = "active_gene_explorer_subtab", default = "browser"),
  # Screen Planning (Methods & References became the planner's Methods section).
  screen_coverage = list(nested_id = "screen_coverage_tabs", field = "active_screen_coverage_subtab", default = "planner")
)

# Subtabs that carry their own inner tabs (third level): inner navset id -> state field, default, and the module
# subtab that holds it.
navigation_inner_navsets <- list(
  go_analysis_tabs = list(field = "active_go_analysis_subtab", default = "gene_set", top = "gene_explorer", parent = "go"),
  comparative_tabs = list(field = "active_comparative_subtab", default = "single_species", top = "gene_explorer", parent = "comparative")
)

# Tabs retired by the module reorganization, kept so old links and bookmarks still open the right place.
navigation_legacy_targets <- list(
  overview = c(top = "project", nested = "samples"),
  project_data = c(top = "project", nested = "samples"),
  pipeline = c(top = "run", nested = "profiles"),
  preprocessing = c(top = "run", nested = ""),
  results_exports = c(top = "project", nested = "samples"),
  logs_provenance = c(top = "project", nested = "provenance"),
  help = c(top = "project", nested = "help"),
  genome_browser = c(top = "gene_explorer", nested = "browser"),
  gene_domain_explorer = c(top = "gene_explorer", nested = "gene_domain"),
  go_analysis = c(top = "gene_explorer", nested = "go"),
  comparative_species = c(top = "gene_explorer", nested = "comparative")
)
# Old top tabs whose own subtabs became inner tabs.
navigation_legacy_inner <- c(go_analysis = "go_analysis_tabs", comparative_species = "comparative_tabs")
navigation_legacy_subtabs <- list(
  project_data = c(project = "samples", samples = "samples", reference = "reference", files = "samples")
)

essentiality_hybrid_input_stages <- c("combine", "combined_summary", "classifier")
navigation_merged_run_subtabs <- c("hit_files", "summary_tables", "status")
navigation_qc_legacy <- c(mapping_qc = "reads_mapping", summary_qc = "insertions", library_diagnostics = "saturation_bias", library_correlation = "agreement")

navigation_translate_legacy <- function(top, nested_value = "") {
  top <- as.character(top %||% "")[1L]
  nested_value <- as.character(nested_value %||% "")[1L]
  # MidLC moved from Essentiality ("Normalize & Choose Target") to the Normalization module, and from its own subtab
  # into Site Normalization as the "MidLC common saturation" method.
  if (!is.na(top) && identical(top, "essentiality") && nested_value %in% c("normalize", "evaluate"))
    return(list(top = "normalization", nested_value = "site"))
  if (!is.na(top) && identical(top, "normalization") && identical(nested_value, "midlc"))
    return(list(top = "normalization", nested_value = "site"))
  # Quality Control was regrouped around the QC table (Overview) into themed subtabs.
  if (!is.na(top) && identical(top, "quality_control") && nested_value %in% names(navigation_qc_legacy))
    return(list(top = "quality_control", nested_value = navigation_qc_legacy[[nested_value]]))
  # Run > Hit Files, Summary Tables and Status merged into Run > Mapping.
  if (!is.na(top) && identical(top, "run") && nested_value %in% navigation_merged_run_subtabs)
    return(list(top = "run", nested_value = "mapping"))
  # Project > Overview was retired; its stage summary lives in Run > Mapping.
  if (!is.na(top) && identical(top, "project") && identical(nested_value, "overview"))
    return(list(top = "project", nested_value = "samples"))
  # Combine Parents, Combined Summary and Classifier became the hybrid model's inputs, built automatically by Run.
  if (!is.na(top) && identical(top, "essentiality") && nested_value %in% c(essentiality_hybrid_input_stages, "hybrid_inputs"))
    return(list(top = "essentiality", nested_value = "run"))
  if (is.na(top) || !top %in% names(navigation_legacy_targets))
    return(list(top = top, nested_value = nested_value))
  target <- navigation_legacy_targets[[top]]
  if (top %in% names(navigation_legacy_inner))
    return(list(top = target[["top"]], nested_value = target[["nested"]], inner_id = navigation_legacy_inner[[top]],
                inner_value = if (nzchar(nested_value)) nested_value else ""))
  nested <- if (nzchar(target[["nested"]])) target[["nested"]] else nested_value
  mapped <- navigation_legacy_subtabs[[top]]
  if (!is.null(mapped) && nzchar(nested_value %||% "") && nested_value %in% names(mapped))
    nested <- mapped[[nested_value]]
  new_top <- if (identical(top, "project_data") && identical(nested, "reference")) "run" else target[["top"]]
  if (identical(new_top, "run") && nested %in% navigation_merged_run_subtabs) nested <- "mapping"
  list(top = new_top, nested_value = nested)
}

navigation_state_defaults <- function(top = "project") {
  top <- navigation_translate_legacy(top)$top
  state <- list(active_top_tab = top)
  for (module in c(navigation_modules, navigation_inner_navsets))
    if (length(module)) state[[module$field]] <- module$default
  state$active_results_view <- "overview"
  state
}

navigation_top_values <- names(navigation_modules)

navigation_top_nested_ids <- unlist(lapply(navigation_modules, function(module) module$nested_id))

navigation_nested_id_for_top <- function(top) {
  top <- as.character(top %||% "")
  if (!length(top) || !nzchar(top[[1L]])) return("")
  if (!top[[1L]] %in% names(navigation_top_nested_ids)) return("")
  unname(navigation_top_nested_ids[[top[[1L]]]])
}

navigation_nested_fields <- c(
  setNames(unlist(lapply(navigation_modules, function(module) module$field)), navigation_top_nested_ids),
  vapply(navigation_inner_navsets, function(inner) inner$field, ""),
  essentiality_results_view = "active_results_view"
)

normalize_navigation_state <- function(value = NULL, top = NULL) {
  defaults <- navigation_state_defaults(top %||% "project")
  if (!is.list(value)) return(defaults)
  for (name in intersect(names(defaults), names(value))) {
    candidate <- as.character(value[[name]] %||% "")
    if (length(candidate) && nzchar(candidate[[1L]])) defaults[[name]] <- candidate[[1L]]
  }
  defaults$active_top_tab <- navigation_translate_legacy(defaults$active_top_tab)$top
  if (defaults$active_essentiality_subtab %in% c("normalize", "evaluate"))
    defaults$active_essentiality_subtab <- "run"
  if (defaults$active_essentiality_subtab %in% c(essentiality_hybrid_input_stages, "hybrid_inputs"))
    defaults$active_essentiality_subtab <- "run"
  if (!defaults$active_screen_coverage_subtab %in% c("planner", "mutants"))
    defaults$active_screen_coverage_subtab <- "planner"
  if (defaults$active_qc_subtab %in% names(navigation_qc_legacy))
    defaults$active_qc_subtab <- navigation_qc_legacy[[defaults$active_qc_subtab]]
  if (!is.null(top) && nzchar(top)) defaults$active_top_tab <- navigation_translate_legacy(top)$top
  defaults
}

navigation_state_from_url <- function(value = NULL) {
  state <- normalize_navigation_state()
  if (!is.list(value)) return(state)
  target <- navigation_translate_legacy(value$top, value$nested_value)
  top <- target$top
  nested_value <- target$nested_value
  nested_id <- navigation_nested_id_for_top(top)
  results_view <- as.character(value$results_view %||% "")
  if (length(top) && top[[1L]] %in% navigation_top_values)
    state <- record_navigation_selection(state, "workspace_tabs", top[[1L]])
  if (nzchar(nested_id) && length(nested_value) && nzchar(nested_value[[1L]]))
    state <- record_navigation_selection(state, nested_id, nested_value[[1L]])
  if (nzchar(target$inner_value %||% ""))
    state <- record_navigation_selection(state, target$inner_id, target$inner_value)
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
  target <- navigation_translate_legacy(top, params$subtab %||% params$nested_value %||% "")
  list(
    top = target$top,
    nested_id = navigation_nested_id_for_top(target$top),
    nested_value = target$nested_value,
    inner_id = target$inner_id %||% "",
    inner_value = target$inner_value %||% "",
    results_view = params$results_view %||% "",
    project_id = params$project %||% params$project_id %||% ""
  )
}

record_navigation_selection <- function(value, input_id, selected) {
  value <- normalize_navigation_state(value)
  selected <- as.character(selected %||% "")
  if (!length(selected) || !nzchar(selected[[1L]])) return(value)
  if (identical(input_id, "workspace_tabs")) {
    value$active_top_tab <- navigation_translate_legacy(selected[[1L]])$top
  } else if (input_id %in% names(navigation_nested_fields)) {
    if (identical(input_id, "screen_coverage_tabs") && !selected[[1L]] %in% c("planner", "mutants"))
      selected <- "planner"
    if (identical(input_id, "qc_tabs") && selected[[1L]] %in% names(navigation_qc_legacy)) selected <- navigation_qc_legacy[[selected[[1L]]]]
    if (identical(input_id, "essentiality_tabs") && selected[[1L]] %in% c("normalize", "evaluate")) selected <- "run"
    if (identical(input_id, "essentiality_tabs") && selected[[1L]] %in% c(essentiality_hybrid_input_stages, "hybrid_inputs")) selected <- "run"
    value[[navigation_nested_fields[[input_id]]]] <- selected[[1L]]
  }
  value
}

navigation_subtab <- function(value, top = NULL) {
  value <- normalize_navigation_state(value)
  top <- top %||% value$active_top_tab
  module <- navigation_modules[[top]]
  if (length(module)) value[[module$field]] else NULL
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
  target <- navigation_translate_legacy(top_tab, nested_value)
  if (!identical(target$top, top_tab)) {
    top_tab <- target$top
    nested_value <- if (nzchar(target$nested_value)) target$nested_value else NULL
    nested_id <- if (is.null(nested_value)) NULL else navigation_nested_id_for_top(top_tab)
  }
  if (nzchar(target$inner_value %||% ""))
    value <- record_navigation_selection(value, target$inner_id, target$inner_value)
  if (!is.null(nested_id) && nested_id %in% names(navigation_inner_navsets) && !is.null(nested_value)) {
    inner <- navigation_inner_navsets[[nested_id]]
    value <- record_navigation_selection(value, nested_id, nested_value)
    top_tab <- inner$top; nested_id <- navigation_nested_id_for_top(inner$top); nested_value <- inner$parent
  }
  value <- record_navigation_selection(value, "workspace_tabs", top_tab)
  if (!is.null(nested_id) && !is.null(nested_value))
    value <- record_navigation_selection(value, nested_id, nested_value)
  if (!is.null(results_view))
    value <- record_navigation_selection(value, "essentiality_results_view", results_view)
  updateNavbarPage(session, "workspace_tabs", selected = value$active_top_tab)
  selected_nested <- navigation_subtab(value, value$active_top_tab)
  selected_id <- navigation_nested_id_for_top(value$active_top_tab)
  session$onFlushed(function() {
    if (nzchar(selected_id) && nzchar(selected_nested %||% ""))
      bslib::nav_select(selected_id, selected = selected_nested, session = session)
    for (inner_id in names(navigation_inner_navsets)) {
      inner <- navigation_inner_navsets[[inner_id]]
      if (identical(value$active_top_tab, inner$top) && identical(selected_nested, inner$parent))
        bslib::nav_select(inner_id, selected = value[[inner$field]], session = session)
    }
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
  if(!reviewed)return(list(top="run",nested="samples"))
  stage<-next_preprocessing_stage(project,repo_root,status)
  values<-c(Reference="reference",Mapping="mapping",`Hit Files`="mapping",`Summary Tables`="mapping",Complete="mapping")
  list(top="run",nested=unname(values[[stage]]%||%"samples"))
}
