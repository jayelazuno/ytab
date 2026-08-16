essentiality_generated_plot_inventory <- function(project) {
  root <- file.path(project$project_root, "classifier")
  files <- if (dir.exists(root)) list.files(root, pattern = "\\.png$", recursive = TRUE, full.names = TRUE) else character()
  files <- files[file.exists(files)]
  if (!length(files)) return(data.frame())
  project_root <- normalizePath(project$project_root, winslash = "/", mustWork = FALSE)
  data.frame(
    file = files,
    filename = basename(files),
    target = basename(dirname(dirname(files))),
    title = tools::file_path_sans_ext(gsub("_", " ", basename(files))),
    served_url = vapply(files, function(path) {
      path <- normalizePath(path, winslash = "/", mustWork = FALSE)
      rel <- substring(path, nchar(project_root) + 2L)
      paste0("ytab-project-output/", diagnostic_encode_relative_path(rel))
    }, ""),
    stringsAsFactors = FALSE
  )
}

essentiality_generated_plot_cards <- function(project) {
  plots <- essentiality_generated_plot_inventory(project)
  if (!nrow(plots)) return(tags$p(class = "text-muted", "No classifier-generated plot files are available."))
  plots <- plots[order(plots$target, plots$filename), , drop = FALSE]
  tags$div(class = "ytab-plot-grid", lapply(seq_len(nrow(plots)), function(i) {
    tags$article(
      class = "ytab-plot-card ytab-release-card",
      title = plots$filename[[i]],
      tags$div(class = "ytab-plot-card-header",
               tags$h4(paste("Generated:", plots$title[[i]])),
               tags$span(class = "ytab-status-badge", "Static generated image")),
      tags$p(class = "text-muted", paste("Target:", plots$target[[i]])),
      tags$div(class = "ytab-diagnostic-preview",
        tags$img(src = plots$served_url[[i]], alt = plots$filename[[i]],
                 loading = "lazy", style = "width:100%;max-height:520px;object-fit:contain")
      ),
      tags$details(tags$summary("Show filename"), tags$code(plots$filename[[i]]))
    )
  }))
}
