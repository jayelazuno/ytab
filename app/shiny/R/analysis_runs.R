# Shared helpers for analyses run through `ytab` commands that write one run folder with a manifest.json
# (library correlation, normalization, HMM, multi-condition, genetic interactions).

# Run folders directly under `dir` that contain `marker`, newest first.
ytab_analysis_runs <- function(dir, marker = "manifest.json") {
  empty <- data.frame(run_id = character(), path = character(), modified = as.POSIXct(character()), stringsAsFactors = FALSE)
  if (!nzchar(dir %||% "") || !dir.exists(dir)) return(empty)
  manifests <- file.path(list.dirs(dir, recursive = FALSE), marker)
  manifests <- manifests[file.exists(manifests)]
  if (!length(manifests)) return(empty)
  paths <- dirname(manifests)
  runs <- data.frame(run_id = basename(paths), path = paths, modified = file.mtime(manifests), stringsAsFactors = FALSE)
  runs[order(runs$modified, decreasing = TRUE), , drop = FALSE]
}

# Readable names for run folders: "<kind>_<YYYYMMDD>_<HHMMSS>" becomes "QC table, 21 Oct 2026, 13:13"; other folder
# names are shown with underscores as spaces and the folder's date.
ytab_run_kinds <- c(qc = "QC table", libcor = "Library correlation", norm = "Normalization", essentiality = "Essentiality",
                    multicondition = "Multi-condition", gi = "Genetic interactions")

ytab_run_label <- function(run_id, modified = NULL) {
  m <- regmatches(run_id, regexec("^(.*)_([0-9]{8})_([0-9]{6})$", run_id))[[1]]
  when <- if (length(m) == 4L) as.POSIXct(paste(m[3], m[4]), format = "%Y%m%d %H%M%S") else modified
  kind <- if (length(m) == 4L) unname(ytab_run_kinds[m[2]] %||% NA) else NA
  if (length(m) == 4L && is.na(kind)) kind <- gsub("_", " ", m[2])
  if (length(m) != 4L) kind <- gsub("_", " ", run_id)
  stamp <- if (length(when) && !is.na(when)) format(when, "%d %b %Y, %H:%M") else ""
  trimws(paste0(kind, if (nzchar(stamp)) paste0(", ", sub("^0", "", stamp))))
}

ytab_run_choices <- function(runs) {
  if (!nrow(runs)) return(character())
  stats::setNames(runs$path, vapply(seq_len(nrow(runs)), function(i) ytab_run_label(runs$run_id[[i]], runs$modified[[i]]), ""))
}

ytab_new_run_dir <- function(dir, prefix = "run") file.path(dir, paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))

ytab_read_analysis_tsv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA", "None"))
}

ytab_read_manifest <- function(run_dir) read_json_safely(file.path(run_dir, "manifest.json"))

# House plot style: bars with black borders, black bold axis titles and tick labels, x and y axis lines drawn (not on
# heatmaps). Applied to every app plot, on screen (ytab_ggplotly) and in downloads (ytab_save_ggplot).
YTAB_BAR_BORDER_MM <- .8
YTAB_AXIS_LINE_MM <- .8
YTAB_TICK_MM <- .6

ytab_style_ggplot <- function(plot) {
  if (!inherits(plot, "ggplot")) return(plot)
  tiles <- any(vapply(plot$layers, function(l) inherits(l$geom, c("GeomTile", "GeomRaster")), TRUE))
  for (layer in plot$layers) {
    if (inherits(layer$geom, c("GeomCol", "GeomBar", "GeomRect")) && !inherits(layer$geom, "GeomTile") && is.null(layer$mapping$colour)) {
      layer$aes_params$colour <- "black"; layer$aes_params$linewidth <- YTAB_BAR_BORDER_MM
    }
  }
  plot + ggplot2::theme(axis.title = ggplot2::element_text(face = "bold", colour = "black"),
                        axis.text = ggplot2::element_text(face = "bold", colour = "black"),
                        axis.line = if (tiles) ggplot2::element_blank() else ggplot2::element_line(colour = "black", linewidth = YTAB_AXIS_LINE_MM),
                        axis.ticks = if (tiles) ggplot2::element_blank() else ggplot2::element_line(colour = "black", linewidth = YTAB_TICK_MM),
                        panel.border = if (tiles) plot$theme$panel.border %||% ggplot2::element_blank() else ggplot2::element_blank())
}

# Bold axis titles and tick labels in a plotly figure (this plotly.js has no font weight, so the text is wrapped in <b>)
# and black axis lines and ticks.
# native = TRUE for figures built directly with plot_ly: axis lines are switched on and automatic tick labels use a bold
# font face (they have no tick text to wrap).
ytab_plotly_style <- function(widget, axis_lines = TRUE, native = FALSE) {
  built <- plotly::plotly_build(widget)
  bold <- function(x) ifelse(grepl("^<b>", x), x, paste0("<b>", x, "</b>"))
  for (name in grep("^[xy]axis[0-9]*$", names(built$x$layout), value = TRUE)) {
    axis <- built$x$layout[[name]]
    title <- if (is.list(axis$title)) axis$title$text else axis$title
    if (length(title) && nzchar(title %||% "")) { if (is.list(axis$title)) axis$title$text <- bold(title) else axis$title <- bold(title) }
    if (length(axis$ticktext)) axis$ticktext <- bold(as.character(axis$ticktext))
    axis$tickfont <- utils::modifyList(axis$tickfont %||% list(), c(list(color = "black"), if (native && !length(axis$ticktext)) list(family = "Arial Black, Arial, sans-serif")))
    if (native) axis <- utils::modifyList(axis, list(showline = TRUE, ticks = "outside"))
    # ggplotly keeps axis.line's colour but not its width; set it in pixels (1 mm = 96 / 25.4 px).
    if (isTRUE(axis_lines) && !identical(axis$showline, FALSE)) axis <- utils::modifyList(axis, list(linecolor = "black", tickcolor = "black",
      linewidth = YTAB_AXIS_LINE_MM * 96 / 25.4, tickwidth = YTAB_TICK_MM * 96 / 25.4))
    built$x$layout[[name]] <- axis
  }
  # Heatmap cells: ggplotly drops tile borders, so cells are separated by gaps over a dark plot background
  # (not for heatmaps with more than 80 rows or columns, whose cells would vanish behind the gaps).
  for (i in seq_along(built$x$data)) {
    trace <- built$x$data[[i]]
    if (identical(trace$type, "heatmap") && length(dim(trace$z)) == 2L && max(dim(trace$z)) <= 80L) {
      built$x$data[[i]]$xgap <- 1.5; built$x$data[[i]]$ygap <- 1.5; built$x$layout$plot_bgcolor <- "#595959"
    }
  }
  built
}

ytab_ggplotly <- function(plot, tooltip = "text") ytab_plotly_style(plotly::ggplotly(ytab_style_ggplot(plot), tooltip = tooltip))

# Writes a ggplot to PNG or PDF according to the file extension, in the house plot style.
ytab_save_ggplot <- function(file, plot, width = 9, height = 7) {
  device <- if (grepl("\\.pdf$", file, ignore.case = TRUE)) "pdf" else "png"
  ggplot2::ggsave(file, ytab_style_ggplot(plot), device = device, width = width, height = height, dpi = 200, bg = "white")
}

ytab_plot_downloads_ui <- function(id, tables = c(table = "Download table")) {
  tags$div(class = "ytab-actions",
    downloadButton(paste0(id, "_download_png"), "Download plot (PNG)"),
    downloadButton(paste0(id, "_download_pdf"), "Download plot (PDF)"),
    lapply(names(tables), function(name) downloadButton(paste0(id, "_download_", name), tables[[name]])))
}

# Methodology reference for a method adapted from another tool; the only place such tools are named in the app.
ytab_reference <- function(...) tags$p(class = "ytab-reference", tags$b("Reference. "), ...)

ytab_methods_drawer <- function(...) tags$details(class = "ytab-technical-details", tags$summary("Methods"), ...)
