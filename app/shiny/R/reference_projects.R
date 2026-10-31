# Reference projects: an opt-in comparison set. A user picks other projects once and every tab that can overlay them
# uses the same choice, so the comparison follows them from Quality Control to the Gene Explorer.
#
# Nothing is computed here. Each overlay reads the tables the command line already wrote for that project (its newest
# `ytab qc` run, its mapping statistics), exactly as the open project's own tab reads its own. A project with no such
# table simply contributes nothing.
#
# Comparability is the user's to judge, so every overlay carries the reference's species, transposon and library count
# rather than presenting a bare number.

YTAB_REFERENCE_COLOURS <- c("#1A73E8", "#8E24AA", "#00897B", "#E8710A", "#5D4037", "#3949AB")

# Past this many bars on one axis the libraries stop being readable, so the plot splits into a panel per project.
YTAB_REFERENCE_MAX_BARS <- 40L

ytab_reference_display_choices <- c("Automatic" = "auto", "Bars beside mine" = "bars",
                                    "One panel per project" = "panels", "Median and spread lines" = "lines")

# Automatic: bars while they still fit on one axis, panels once they do not.
ytab_reference_mode <- function(mode, own_n, reference_n) {
  mode <- mode %||% "auto"
  if (!identical(mode, "auto")) return(mode)
  if (own_n + reference_n <= YTAB_REFERENCE_MAX_BARS) "bars" else "panels"
}

# The fill scale covers the open project's QC flags and one colour per reference project.
ytab_reference_fill_values <- function(projects) c(YTAB_FLAG_COLOURS,
  stats::setNames(rep(YTAB_REFERENCE_COLOURS, length.out = max(1L, length(projects))), projects)[projects])

# Shared across tabs: the project configs the user picked, in the order they picked them.
ytab_reference_state <- function() reactiveVal(character())

# Every discoverable project except the one open, same species first so the usual comparison is at the top.
ytab_reference_choices <- function(projects, active_project) {
  if (!is.data.frame(projects) || !nrow(projects)) return(character())
  keep <- projects$project_config != as.character(active_project$project_config %||% "")
  projects <- projects[keep, , drop = FALSE]
  if (!nrow(projects)) return(character())
  species <- as.character(active_project$species %||% "")
  same <- projects$species == species & nzchar(species)
  projects <- projects[order(!same, tolower(projects$project_id)), , drop = FALSE]
  labels <- sprintf("%s (%s, %d libraries)", projects$project_id, projects$species, projects$included_count)
  group <- factor(ifelse(projects$species == species & nzchar(species), "Same species", "Other species"),
                  levels = c("Same species", "Other species"))
  Filter(length, split(stats::setNames(projects$project_config, labels), group))
}

# What a reference project is, for the label and the comparability note. Cheap enough to call per render.
ytab_reference_meta <- function(config, repo_root) {
  info <- tryCatch(read_project_summary(config, repo_root), error = function(e) NULL)
  if (is.null(info)) return(NULL)
  list(project_config = config, project_id = as.character(info$project_id %||% ""),
       label = project_display_label(info), species = as.character(info$species %||% ""),
       transposon = as.character(info$transposon %||% ""), project_root = info$project_root,
       samples = info$samples)
}

ytab_reference_metas <- function(configs, repo_root) {
  out <- lapply(as.character(configs %||% character()), ytab_reference_meta, repo_root = repo_root)
  Filter(Negate(is.null), out)
}

# The project's newest `ytab qc` run, or "" when it has none.
ytab_reference_qc_run <- function(project_root) {
  runs <- ytab_analysis_runs(file.path(project_root, "qc"))
  runs <- runs[file.exists(file.path(runs$path, "qc_table.tsv")), , drop = FALSE]
  if (!nrow(runs)) "" else runs$path[[1L]]
}

# One row per reference library for a QC-table metric: value, role and which project it came from.
# `roles` restricts the reference to the roles shown for the open project; a reference holding none of them
# contributes all of its libraries instead, and says so.
ytab_reference_metric_values <- function(metas, repo_root, metric, roles = character()) {
  rows <- lapply(metas, function(m) {
    run <- ytab_reference_qc_run(m$project_root)
    if (!nzchar(run)) return(NULL)
    t <- ytab_read_analysis_tsv(file.path(run, "qc_table.tsv"))
    if (!nrow(t) || !metric %in% names(t)) return(NULL)
    value <- suppressWarnings(as.numeric(t[[metric]]))
    role <- if ("role" %in% names(t)) tolower(as.character(t$role)) else rep("library", nrow(t))
    keep <- is.finite(value)
    # Roles were asked for but this project has none of them: fall back to all its libraries, and record that it did.
    fell_back <- length(roles) > 0L && !any(keep & role %in% roles)
    if (length(roles) && !fell_back) keep <- keep & role %in% roles
    if (!any(keep)) return(NULL)
    data.frame(.project = m$label, .species = m$species, .transposon = m$transposon,
               sample = as.character(t$sample)[keep],
               display_name = if ("display_name" %in% names(t)) as.character(t$display_name)[keep] else as.character(t$sample)[keep],
               role = role[keep], value = value[keep], fell_back = fell_back, stringsAsFactors = FALSE)
  })
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# Mapping statistics of each reference library, for the Reads & Mapping read-count plot.
ytab_reference_mapping_values <- function(metas, repo_root, column) {
  rows <- lapply(metas, function(m) {
    d <- tryCatch(read_qc_csvs(file.path(m$project_root, "mapfastq"), "mapping_stats.*\\.csv$"), error = function(e) data.frame())
    if (!is.data.frame(d) || !nrow(d) || !column %in% names(d)) return(NULL)
    value <- suppressWarnings(as.numeric(d[[column]]))
    if (!any(is.finite(value))) return(NULL)
    sample <- if ("sample" %in% names(d)) as.character(d$sample) else basename(dirname(d$.source))
    data.frame(.project = m$label, .species = m$species, .transposon = m$transposon,
               sample = sample[is.finite(value)], display_name = sample[is.finite(value)],
               role = "library", value = value[is.finite(value)], fell_back = FALSE, stringsAsFactors = FALSE)
  })
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# One band per reference project: the median its libraries sit at, and the spread around it. The interquartile range
# needs four libraries to mean anything; below that the band is the full range and the label says so.
ytab_reference_summary <- function(values) {
  if (!is.data.frame(values) || !nrow(values)) return(data.frame())
  parts <- split(values, values$.project)
  out <- lapply(names(parts), function(name) {
    v <- parts[[name]]$value; n <- length(v)
    spread <- if (n >= 4L) stats::quantile(v, c(.25, .75), names = FALSE) else range(v)
    data.frame(.project = name, species = parts[[name]]$.species[[1L]], transposon = parts[[name]]$.transposon[[1L]],
               libraries = n, median = stats::median(v), lower = spread[[1L]], upper = spread[[2L]],
               spread_label = if (n >= 4L) "interquartile range" else "full range",
               fell_back = isTRUE(parts[[name]]$fell_back[[1L]]), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# The legend names the project and how many libraries it stands for. The medians and spreads go in the caption
# instead: a legend long enough to carry them is cut off once two or three projects are chosen.
ytab_reference_legend_labels <- function(summary) sprintf("%s (n = %d)", summary$.project, summary$libraries)

# Draw the reference projects behind the open project's bars: a shaded spread and a median line each.
# `note` names the quantity when the plot itself shows more than one (the read-count plot). `horizontal` says which
# axis carries the value: upright bars put it on y, bars drawn along y put it on x.
ytab_reference_overlay <- function(p, summary, note = "", horizontal = FALSE) {
  if (!is.data.frame(summary) || !nrow(summary)) return(p)
  summary$.label <- ytab_reference_legend_labels(summary)
  summary$.label <- factor(summary$.label, levels = summary$.label)
  colours <- stats::setNames(rep(YTAB_REFERENCE_COLOURS, length.out = nrow(summary)), levels(summary$.label))
  # The spread is drawn as its two bounds rather than a shaded rectangle: a rectangle spanning the whole category axis
  # needs infinite limits, which do not survive the conversion to an interactive plot, and the bounds carry the
  # project's own colour instead of an anonymous grey.
  bounds <- data.frame(.label = rep(summary$.label, 2L), value = c(summary$lower, summary$upper))
  rule <- function(data, mapping, ...) if (horizontal) ggplot2::geom_vline(data = data, inherit.aes = FALSE, mapping, ...)
                                       else ggplot2::geom_hline(data = data, inherit.aes = FALSE, mapping, ...)
  intercept <- function(column) if (horizontal) ggplot2::aes(xintercept = .data[[column]], colour = .label)
                                else ggplot2::aes(yintercept = .data[[column]], colour = .label)
  p +
    rule(bounds, intercept("value"), linetype = "13", linewidth = .4, alpha = .75) +
    rule(summary, intercept("median"), linetype = "22", linewidth = .7) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::labs(colour = if (nzchar(note)) paste0("Reference projects (", note, ")") else "Reference projects")
}

# Rows of the same kind from each reference project, for the plots that draw a curve or a group of bars per library
# rather than one value per library. `reader` is the function the tab already uses for its own project, called with
# the reference project in place of it, so the two sides are always read the same way. Library names are prefixed
# with the project, since names repeat between projects and the plots key their series on them.
ytab_reference_curves <- function(metas, reader) {
  rows <- lapply(metas, function(m) {
    d <- tryCatch(reader(m), error = function(e) NULL)
    if (!is.data.frame(d) || !nrow(d)) return(NULL)
    # Readers that summarize a project (genes hit per role, signal per chromosome) have no sample column; only the
    # per-library ones can be labelled by library.
    if ("sample" %in% names(d)) {
      d <- qc_annotate(d, qc_library_info(m))
      d$.library <- paste0(m$label, " · ", d$.library)
    }
    d$.project <- m$label
    d
  })
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, rows)
  if (!length(rows)) return(data.frame())
  shared <- Reduce(intersect, lapply(rows, names))
  do.call(rbind, lapply(rows, function(x) x[, shared, drop = FALSE]))
}

# One reference project's count-distribution curves, from its own newest QC run.
ytab_reference_distribution <- function(meta, kind) {
  run <- ytab_reference_qc_run(meta$project_root)
  if (!nzchar(run)) return(NULL)
  ytab_read_analysis_tsv(file.path(run, "plots", paste0(kind, ".tsv.gz")))
}

# The open project's rows with the reference projects' appended, over the columns they share. Plots that key their
# series on `.library` need nothing else: the reference libraries simply become more series. `.project` is carried on
# both sides, so a plot that groups by role or pool can keep the projects apart instead of averaging across them.
ytab_reference_bind <- function(d, reference, project_label = "This project") {
  if (!is.data.frame(d) || !nrow(d)) return(d)
  if (!".project" %in% names(d)) d$.project <- project_label
  if (!is.data.frame(reference) || !nrow(reference)) return(d)
  shared <- intersect(names(d), names(reference))
  if (!length(shared)) return(d)
  rbind(d[, shared, drop = FALSE], reference[, shared, drop = FALSE])
}

# Each reference project's newest essentiality run, as its genes' HMM scores and calls. The genes are not matched
# between projects: the point is the shape of another library's calls, not a gene-by-gene comparison.
ytab_reference_essentiality <- function(metas) {
  rows <- lapply(metas, function(m) {
    runs <- tryCatch(ytab_analysis_runs(file.path(m$project_root, "essentiality")), error = function(e) data.frame())
    runs <- runs[file.exists(file.path(runs$path, "models", "genes.tsv")), , drop = FALSE]
    if (!nrow(runs)) return(NULL)
    d <- ytab_read_analysis_tsv(file.path(runs$path[[1L]], "models", "genes.tsv"))
    if (!nrow(d) || !all(c("hmm_score", "hmm_call") %in% names(d))) return(NULL)
    data.frame(.project = m$label, .species = m$species,
               hmm_score = suppressWarnings(as.numeric(d$hmm_score)),
               hmm_call = as.character(d$hmm_call), stringsAsFactors = FALSE)
  })
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# Each reference project's newest normalization run, as the same tables this tab reads for its own project. Sample
# names are prefixed with the project, so they become samples of their own in every view and in the legend; what the
# methods do to another library's counts is the point, not the values themselves.
ytab_reference_normalization <- function(metas, reader) {
  tables <- lapply(metas, function(m) {
    runs <- tryCatch(ytab_analysis_runs(file.path(m$project_root, "normalization")), error = function(e) data.frame())
    runs <- runs[dir.exists(file.path(runs$path, "experiment", "methods")), , drop = FALSE]
    if (!nrow(runs)) return(NULL)
    d <- tryCatch(reader(runs$path[[1L]]), error = function(e) NULL)
    if (!is.list(d)) return(NULL)
    lapply(d, function(x) {
      if (!is.data.frame(x) || !nrow(x) || !"display_name" %in% names(x)) return(NULL)
      x$display_name <- paste0(m$label, " · ", x$display_name)
      x$.project <- m$label
      x
    })
  })
  tables <- Filter(Negate(is.null), tables)
  if (!length(tables)) return(list())
  names_all <- unique(unlist(lapply(tables, names)))
  stats::setNames(lapply(names_all, function(key) {
    parts <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, lapply(tables, `[[`, key))
    if (!length(parts)) return(data.frame())
    shared <- Reduce(intersect, lapply(parts, names))
    do.call(rbind, lapply(parts, function(x) x[, shared, drop = FALSE]))
  }), names_all)
}

# The open project's normalization tables with the reference projects' rows appended, table by table.
ytab_reference_normalization_bind <- function(data, reference) {
  if (!length(reference)) return(data)
  for (key in intersect(names(data), names(reference))) {
    own <- data[[key]]; other <- reference[[key]]
    if (!is.data.frame(own) || !is.data.frame(other) || !nrow(other)) next
    shared <- intersect(names(own), names(other))
    if (!length(shared) || !"display_name" %in% shared) next
    data[[key]] <- rbind(own[, shared, drop = FALSE], other[, shared, drop = FALSE])
  }
  data
}

# Fills for a plot coloured by project rather than by QC flag: the open project keeps the app's own blue, each
# reference its colour from the shared list, so a project is the same colour wherever it is drawn.
ytab_reference_project_fills <- function(projects, own) {
  others <- setdiff(projects, own)
  c(stats::setNames("#1A73E8", own),
    stats::setNames(rep(YTAB_REFERENCE_COLOURS, length.out = max(1L, length(others))), others)[others])
}

# A grouping label that stays distinct across projects. Within one project nothing changes.
ytab_reference_group <- function(group, project) {
  if (length(unique(project)) < 2L) return(as.character(group))
  paste0(project, " · ", group)
}

# The plotted data with every reference library appended, so a download holds all the plot draws on and not just the
# bands' summary. Columns the two sides do not share are left empty.
ytab_reference_download <- function(data, reference, project_label) {
  values <- reference$values
  if (!is.data.frame(data) || !is.data.frame(values) || !nrow(values)) return(data)
  data$project <- as.character(project_label)
  values$project <- values$.project
  values <- values[, !startsWith(names(values), "."), drop = FALSE]
  fill <- function(d, cols) { for (col in setdiff(cols, names(d))) d[[col]] <- NA; d[, cols, drop = FALSE] }
  cols <- union(names(data), names(values))
  rbind(fill(data, cols), fill(values, cols))
}

# What the comparison is for, said plainly under every plot that draws one. It shows what other libraries look like,
# whether Hermes, Ac/Ds, piggyBac or another element; it is not a statistical comparison, and the numbers of two
# projects are not on a common footing — transposon, density, depth, strain and protocol all differ.
YTAB_REFERENCE_DISCLAIMER <- paste("Shown for a sense of what other libraries look like, not as a numerical",
  "comparison: transposon, insertion density, sequencing depth, strain and protocol all differ between projects.")

ytab_reference_caption <- function(summary, active_species) {
  if (!is.data.frame(summary) || !nrow(summary)) return("")
  parts <- sprintf("%s: %s libraries, median %s (%s %s–%s)%s", summary$.project, summary$libraries,
                   qc_format_value(summary$median), summary$spread_label,
                   qc_format_value(summary$lower), qc_format_value(summary$upper),
                   ifelse(summary$fell_back, "; all its libraries, it has none of the roles shown here", ""))
  paste0("Reference projects. ", paste(parts, collapse = "; "), ". ", YTAB_REFERENCE_DISCLAIMER)
}

# The same note for the plots that draw reference libraries without a median and spread to report (curves, grouped
# bars, gene tracks): which projects are drawn, and what the comparison is for.
ytab_reference_note <- function(metas) {
  if (!length(metas)) return("")
  who <- vapply(metas, function(m) sprintf("%s (%s%s)", m$label, m$species,
    if (nzchar(m$transposon %||% "")) paste0(", ", m$transposon) else ""), "")
  paste0("Reference projects: ", paste(who, collapse = "; "), ". ", YTAB_REFERENCE_DISCLAIMER)
}

# The selector itself. Opt-in but in plain sight: a control hidden behind a disclosure at the foot of a long options
# panel is one nobody finds. A chooser rather than a list of tick boxes, so the projects are out of the way until
# they are wanted; clicking a chosen one in the list removes it again.
# `display_modes` offers the choice between bars, panels and lines. A view that has only one way of drawing a
# reference (the gene track) leaves it out.
ytab_reference_projects_ui <- function(id, title = "Compare with other projects", display_modes = TRUE) {
  tags$div(class = "ytab-reference-projects",
    tags$hr(),
    tags$label(class = "control-label", title),
    uiOutput(paste0(id, "_reference_selector")),
    if (isTRUE(display_modes)) uiOutput(paste0(id, "_reference_display")),
    helpText("Optional. Choose a project to draw its libraries beside this project's; click it again in the list,",
             "or its ×, to remove them. The choice stays as you move between tabs.",
             "It is there to show what other libraries look like — Hermes, Ac/Ds, piggyBac or another element —",
             "not to compare numbers between projects."))
}

# Selectize hides what is already chosen and treats a second click as a no-op, so a chosen project can only be
# dropped from the chip. Keeping the chosen rows in the list and turning a click on one into a removal makes the same
# click do both jobs, which is what a user expects of a list they are ticking through.
YTAB_REFERENCE_TOGGLE_JS <- "function() {
  var self = this;
  self.$dropdown.on('mousedown', '[data-value]', function(event) {
    var value = this.getAttribute('data-value');
    if (self.items.indexOf(value) === -1) return;
    event.preventDefault(); event.stopPropagation();
    self.removeItem(value); self.refreshOptions(false);
  });
}"

# Renders the selector and keeps `state` in step with it. Returns nothing; read the choice from `state()`.
ytab_reference_server <- function(id, input, output, session, repo_root, active, projects, state, display_modes = TRUE) {
  control <- paste0(id, "_reference_projects")
  output[[paste0(id, "_reference_selector")]] <- renderUI({
    choices <- ytab_reference_choices(tryCatch(projects(), error = function(e) data.frame()), active())
    if (!length(unlist(choices))) return(tags$p(class = "text-muted", "No other project is available to compare with."))
    selectizeInput(control, NULL, choices, selected = isolate(state()), multiple = TRUE,
      options = list(placeholder = "Choose projects to compare with…", hideSelected = FALSE,
                     plugins = list("remove_button"), onInitialize = I(YTAB_REFERENCE_TOGGLE_JS)))
  })
  output[[paste0(id, "_reference_display")]] <- renderUI({
    if (!isTRUE(display_modes) || !length(state())) return(NULL)
    selectInput(paste0(id, "_reference_mode"), "Show reference libraries as", ytab_reference_display_choices,
                selected = isolate(input[[paste0(id, "_reference_mode")]]) %||% "auto")
  })
  # Both outputs can be scrolled out of sight inside a long control panel; without this Shiny suspends them and the
  # controls never appear.
  for (name in paste0(id, c("_reference_selector", "_reference_display")))
    outputOptions(output, name, suspendWhenHidden = FALSE)
  observeEvent(input[[control]], state(as.character(input[[control]] %||% character())), ignoreNULL = FALSE)
}
