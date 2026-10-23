# Quality Control -> Library Agreement: runs `ytab library-correlation` on this project's libraries, or with other
# projects' libraries, and summarizes the result (cards and pair table); qc_module.R draws the plots.

# The Library Agreement UI is in qc_module.R (qc_agreement_run_ui and the plot menu).

libcor_project_dir <- function(config) basename(dirname(dirname(config)))

# Display names for library-correlation names ("<project folder>::<sample>"): each project's sample-sheet display
# name, read when shown, so names edited later in Project > Samples carry over to earlier results. With several
# projects, a short project name (its first word) comes first.
# Roles (parent / control / treated ...) of library-correlation libraries, from each project's sample sheet
# (library_role, else condition), else the condition recorded in the run's samples.tsv.
libcor_roles <- function(names, projects_root, samples = NULL) {
  names <- as.character(names); project <- sub("::.*$", "", names); sample <- sub("^.*::", "", names)
  out <- if (is.data.frame(samples) && "condition" %in% names(samples)) as.character(samples$condition[match(names, samples$name)]) else rep(NA_character_, length(names))
  for (p in unique(project[grepl("::", names)])) {
    sheet <- tryCatch(utils::read.csv(file.path(projects_root, p, "config", "sample_sheet.csv"), stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (!is.data.frame(sheet) || !"sample" %in% names(sheet)) next
    col <- intersect(c("library_role", "condition", "guessed_condition"), names(sheet))
    if (!length(col)) next
    i <- which(project == p); k <- match(sample[i], sheet$sample); v <- as.character(sheet[[col[[1]]]])[k]
    out[i] <- ifelse(!is.na(v) & nzchar(v), v, out[i])
  }
  out <- tolower(trimws(out)); out[is.na(out) | !nzchar(out)] <- "library"
  out
}

libcor_display_names <- function(names, projects_root) {
  names <- as.character(names); project <- sub("::.*$", "", names); sample <- sub("^.*::", "", names); out <- sample
  for (p in unique(project[grepl("::", names)])) {
    i <- which(project == p); out[i] <- ytab_sample_labels(sample[i], file.path(projects_root, p))
  }
  if (length(unique(project)) > 1L) out <- paste(sub("[_ -].*$", "", project), out)
  make.unique(out, sep = " ")
}

libcor_matrix <- function(run_dir) {
  m <- ytab_read_analysis_tsv(file.path(run_dir, "correlation_matrix.tsv"))
  if (!nrow(m)) return(NULL)
  values <- as.matrix(m[, -1, drop = FALSE]); rownames(values) <- m$sample
  values
}

# Library labels for plots: the sample within one project; with several projects, a short project name (its first word),
# the sample and, when samples.tsv gives them, the pool and condition (e.g. "Gale SRR11561508 (Cg1 parent)").
libcor_short_names <- function(names, samples = NULL) {
  names <- as.character(names); sources <- sub("::.*$", "", names); sample <- sub("^.*::", "", names)
  if (length(unique(sources)) == 1L) return(sample)
  out <- paste(sub("[_ -].*$", "", sources), sample)
  if (is.data.frame(samples) && nrow(samples) && "name" %in% names(samples)) {
    i <- match(names, samples$name); extra <- trimws(paste(samples$pool[i] %||% "", samples$condition[i] %||% ""))
    out <- ifelse(!is.na(i) & nzchar(extra), sprintf("%s (%s)", out, extra), out)
  }
  out
}

libcor_heatmap_ggplot <- function(values, pairs) {
  labels <- libcor_short_names(rownames(values))
  long <- data.frame(a = factor(rep(labels, times = ncol(values)), levels = rev(labels)),
                     b = factor(rep(labels, each = nrow(values)), levels = labels), rho = as.vector(values))
  plot <- ggplot2::ggplot(long, ggplot2::aes(b, a, fill = rho)) + ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_viridis_c("Spearman rho", limits = c(min(0, min(long$rho, na.rm = TRUE)), 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = "Hermes libraries: correlation of gene insertion density") +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::coord_fixed() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
  if (nrow(values) <= 20L) plot <- plot + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 3, colour = "white")
  plot
}

qc_library_correlation_server <- function(input, output, session, repo_root, active, active_project_path, projects, jobs, python_bin) {
  runs_dir <- reactive(file.path(active()$project_root, "library_correlation"))
  runs <- reactivePoll(3000, session,
    checkFunc = function() { dir <- tryCatch(runs_dir(), error = function(e) ""); x <- ytab_analysis_runs(dir); paste(x$path, x$modified, collapse = "|") },
    valueFunc = function() ytab_analysis_runs(tryCatch(runs_dir(), error = function(e) "")))

  output$libcor_project_selector <- renderUI({
    found <- projects(); current <- libcor_project_dir(active()$project_config)
    dirs <- vapply(found$project_config, libcor_project_dir, ""); keep <- dirs != current
    choices <- if (any(keep)) stats::setNames(dirs[keep], sprintf("%s (%s)", found$project_id[keep], found$species[keep])) else character()
    selectizeInput("libcor_projects", "Compare with these projects", choices, multiple = TRUE, options = list(placeholder = "Choose projects…"))
  })
  output$libcor_ortholog_selector <- renderUI({
    files <- list.files(file.path(repo_root, "resources", "species", "orthologs"), pattern = "\\.tsv$", full.names = TRUE)
    selectInput("libcor_orthologs", "1:1 ortholog table (needed for several species)", c("None" = "", stats::setNames(files, basename(files))))
  })
  output$libcor_result_selector <- renderUI({
    x <- runs()
    if (!nrow(x)) return(tags$p(class = "text-muted", "No library correlation has been run for this project yet."))
    selectInput("libcor_result", "Result", ytab_run_choices(x))
  })

  observeEvent(input$libcor_run, {
    if (jobs$job_is_running()) { showNotification("Wait for the current job to finish.", type = "warning"); return() }
    current <- libcor_project_dir(active()$project_config)
    selected <- unique(c(current, if (identical(input$libcor_scope, "across")) input$libcor_projects %||% character()))
    if (identical(input$libcor_scope, "across") && length(selected) < 2L) { showNotification("Choose at least one other project to compare with.", type = "error"); return() }
    found <- projects(); species <- unique(found$species[vapply(found$project_config, libcor_project_dir, "") %in% selected])
    samples <- sum(found$included_count[vapply(found$project_config, libcor_project_dir, "") %in% selected])
    if (!length(selected) || samples < 2L) { showNotification("Choose projects with at least two included libraries.", type = "error"); return() }
    if (length(species) > 1L && !nzchar(input$libcor_orthologs %||% "")) {
      showNotification("These projects are from several species. Choose a 1:1 ortholog table.", type = "error"); return()
    }
    output_dir <- ytab_new_run_dir(runs_dir(), "libcor")
    args <- c(file.path(repo_root, "scripts", "local", "ytab_library_correlation.py"), as.vector(rbind("--project", selected)),
              "--metric", input$libcor_metric %||% "density", "--min-reads", as.character(max(1L, as.integer(input$libcor_min_reads %||% 1L))),
              if (nzchar(input$libcor_orthologs %||% "")) c("--orthologs", input$libcor_orthologs), "--output", output_dir)
    tryCatch({
      jobs$start_job(python_bin(), args, active_project_path(), "library_correlation", wd = repo_root, selected_items = selected)
      showNotification("Library correlation started. The result appears here when it finishes.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error"))
  }, ignoreInit = TRUE)

  selected_run <- reactive({ req(nrow(runs()) > 0); path <- input$libcor_result %||% runs()$path[[1L]]; req(path %in% runs()$path); path })
  matrix_values <- reactive(libcor_matrix(selected_run()))
  pairs <- reactive(ytab_read_analysis_tsv(file.path(selected_run(), "pairs.tsv")))
  manifest <- reactive(ytab_read_manifest(selected_run()))

  output$libcor_summary <- renderUI({
    if (!nrow(runs())) return(ytab_empty_state("No library correlation yet", "Run library correlation to compare this project's libraries with each other, or with libraries from other projects."))
    m <- manifest(); req(m); flags <- m$flagged_pairs %||% list()
    if (length(flags)) ytab_warning_banner(sprintf("%d library pairs disagree with their pool labels", length(flags)),
      "Their read correlation at shared sites does not match the pools given in the sample sheet. See the flag column in the pair table.")
  })
  output$libcor_pairs_table <- DT::renderDT({
    p <- pairs(); req(nrow(p) > 0)
    keep <- intersect(c("sample_a", "sample_b", "pair_type", "genes", "rho", "ci_low", "ci_high", "site_jaccard", "shared_site_read_rho", "clone_evidence", "flag"), names(p))
    p <- p[, keep, drop = FALSE]
    root <- dirname(active()$project_root); labels <- libcor_display_names(unique(c(p$sample_a, p$sample_b)), root)
    names_all <- unique(c(p$sample_a, p$sample_b)); p$sample_a <- labels[match(p$sample_a, names_all)]; p$sample_b <- labels[match(p$sample_b, names_all)]
    for (column in intersect(c("rho", "ci_low", "ci_high", "site_jaccard", "shared_site_read_rho"), names(p))) p[[column]] <- round(p[[column]], 3)
    names(p) <- c(sample_a = "Library A", sample_b = "Library B", pair_type = "Pair type", genes = "Genes", rho = "Spearman rho", ci_low = "CI low",
                  ci_high = "CI high", site_jaccard = "Shared sites (Jaccard)", shared_site_read_rho = "Read rho at shared sites",
                  clone_evidence = "Clone evidence", flag = "Flag")[names(p)]
    ytab_datatable(p, page_length = 10L)
  })

  download_name <- function(ext) paste0(active()$project_id, "_library_correlation_", basename(selected_run()), ext)
  output$libcor_download_matrix <- downloadHandler(filename = function() download_name("_matrix.tsv"),
    content = function(file) file.copy(file.path(selected_run(), "correlation_matrix.tsv"), file))
  output$libcor_download_pairs <- downloadHandler(filename = function() download_name("_pairs.tsv"),
    content = function(file) file.copy(file.path(selected_run(), "pairs.tsv"), file))
  invisible(list(runs = runs, selected_run = selected_run, matrix_values = matrix_values, pairs = pairs))
}
