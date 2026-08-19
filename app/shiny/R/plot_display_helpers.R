ytab_plot_output <- function(output_id, height = "medium") {
  plotOutput(output_id, width = "100%", height = ytab_plot_height_px(height))
}

ytab_plot_width_class <- function(choice) {
  choice <- choice %||% "standard"
  if (identical(choice, "wide")) return("ytab-plot-width-wide")
  if (identical(choice, "full")) return("ytab-plot-width-full")
  "ytab-plot-width-standard"
}

ytab_plot_frame <- function(plot_ui, width = "standard", type = "app-rendered") {
  tags$div(
    class = paste("ytab-plot-frame", paste0("ytab-plot-frame-", type),
                  ytab_plot_width_class(width)),
    plot_ui
  )
}

ytab_static_image_ui <- function(src, alt = "", max_height = "760px", show_filename = FALSE) {
  tagList(
    tags$div(
      class = "ytab-static-image",
      tags$img(src = src, alt = alt, loading = "lazy",
               style = paste0("max-height:", max_height, ";"))
    ),
    if (isTRUE(show_filename)) tags$small(class = "text-muted", alt)
  )
}

ytab_static_image_output_card <- function(title, output_id, description = NULL,
                                          height = "large", downloads = NULL,
                                          details = NULL) {
  tags$article(
    class = "ytab-static-image-card ytab-release-card",
    tags$div(
      class = "ytab-plot-card-header",
      tags$div(
        tags$h4(title),
        if (!is.null(description)) tags$p(class = "text-muted", description)
      ),
      if (!is.null(downloads)) tags$div(class = "ytab-actions", downloads)
    ),
    tags$div(
      class = "ytab-static-image ytab-static-image-output",
      imageOutput(output_id, width = "100%", height = ytab_plot_height_px(height))
    ),
    details
  )
}

ytab_generated_file_gallery <- function(files, title_col = "title", url_col = "served_url",
                                        filename_col = "filename") {
  if (!is.data.frame(files) || !nrow(files)) {
    return(ytab_empty_state("No generated files are available."))
  }
  tags$div(
    class = "ytab-generated-file-gallery",
    lapply(seq_len(nrow(files)), function(i) {
      title <- as.character(files[[title_col]][[i]] %||% files[[filename_col]][[i]])
      url <- as.character(files[[url_col]][[i]] %||% "")
      filename <- as.character(files[[filename_col]][[i]] %||% title)
      tags$article(
        class = "ytab-static-image-card ytab-release-card",
        tags$div(class = "ytab-plot-card-header",
                 tags$h4(title)),
        if (nzchar(url)) ytab_static_image_ui(url, filename, "520px"),
        tags$details(tags$summary("Show filename"), tags$code(filename))
      )
    })
  )
}
