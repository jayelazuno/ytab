# Which analyses the open project's experimental design supports, from `ytab design-check`, so that every module
# can say why it has nothing to show instead of leaving the user to find out by running something that cannot work.
#
# The report is read live rather than from a saved file because it is derived from the sample sheet, which the user
# edits inside the app: a stale copy would invite exactly the run the gate exists to prevent.

ytab_design_support <- function(python_bin, repo_root, project_config) {
  empty <- data.frame(analysis = character(), status = character(), why = character(),
                      command = character(), samples = character(), stringsAsFactors = FALSE)
  if (!nzchar(project_config %||% "") || !file.exists(project_config)) return(empty)
  script <- file.path(repo_root, "scripts", "local", "ytab_design_check.py")
  if (!file.exists(script)) return(empty)
  result <- tryCatch(run_process_sync(python_bin(), c(script, "--project-config", project_config, "--json"),
                                      wd = repo_root), error = function(e) list(status = 1L, stdout = ""))
  parsed <- tryCatch(jsonlite::fromJSON(result$stdout, simplifyVector = FALSE), error = function(e) NULL)
  if (!length(parsed)) return(empty)
  data.frame(
    analysis = vapply(parsed, function(x) as.character(x$analysis %||% ""), ""),
    status   = vapply(parsed, function(x) as.character(x$status %||% ""), ""),
    why      = vapply(parsed, function(x) as.character(x$why %||% ""), ""),
    command  = vapply(parsed, function(x) as.character(x$command %||% ""), ""),
    samples  = vapply(parsed, function(x) paste(unlist(x$samples %||% list()), collapse = " "), ""),
    stringsAsFactors = FALSE)
}

# The rows of the report for one analysis, matched on the prefix design-check uses ("Essentiality (HMM), strain X").
ytab_design_rows <- function(support, prefix) {
  if (!nrow(support)) return(support)
  support[startsWith(support$analysis, prefix), , drop = FALSE]
}

# NULL when the design supports the analysis (so the module shows its own empty state: nothing has been run yet),
# otherwise the explanation, which is the one thing a user cannot work out from an empty tab.
ytab_design_gate <- function(support, prefix, title) {
  rows <- ytab_design_rows(support, prefix)
  if (!nrow(support) || any(rows$status == "ready")) return(NULL)
  if (!nrow(rows)) return(NULL)
  reasons <- unique(sprintf("%s: %s", sub("^[^,]*,?\\s*", "", rows$analysis), rows$why))
  reasons <- ifelse(startsWith(reasons, ": "), sub("^: ", "", reasons), reasons)
  ytab_empty_state(title, paste0(
    "This project's experimental design does not support it. ",
    paste(reasons, collapse = " "),
    " The Samples tab shows the roles, strains, pools and conditions this is read from; the design guide",
    " (docs/design_guide.md) gives each analysis's minimum."))
}

# Sequencing runs of one library counted as replicates: design-check reports this as a "needed" step, and every
# analysis downstream of it is wrong until it is done, so it is worth saying wherever a module is used.
ytab_design_merge_warning <- function(support) {
  rows <- ytab_design_rows(support, "Merge sequencing runs")
  rows <- rows[rows$status == "needed", , drop = FALSE]
  if (!nrow(rows)) return(NULL)
  ytab_warning_banner("Some libraries were sequenced more than once",
    paste0(nrow(rows), " group(s) of runs belong to one library. Merge them (", tags$code("ytab merge-runs"),
           ") or they are counted as replicates and every test below overstates its power. ",
           rows$why[[1]], "."))
}
