ytab_page_header <- function(title, subtitle = NULL, actions = NULL) {
  tags$header(
    class = "ytab-release-header",
    tags$div(
      tags$h2(title),
      if (!is.null(subtitle) && nzchar(as.character(subtitle)))
        tags$p(class = "ytab-release-subtitle", subtitle)
    ),
    if (!is.null(actions)) tags$div(class = "ytab-release-actions", actions)
  )
}

ytab_section_header <- function(title, description = NULL) {
  tags$div(
    class = "ytab-section-header",
    tags$h3(title),
    if (!is.null(description) && nzchar(as.character(description)))
      tags$p(description)
  )
}

ytab_status_badge <- function(label, status = "neutral") {
  tags$span(class = paste("ytab-status-badge", paste0("ytab-status-", status)), label)
}

ytab_result_card <- function(title, ..., status = NULL) {
  tags$section(
    class = "ytab-result-card ytab-release-card",
    tags$div(
      class = "ytab-result-heading",
      tags$h4(title),
      if (!is.null(status)) ytab_status_badge(status)
    ),
    ...
  )
}

ytab_metric_card <- function(label, value, caption = NULL, status = "neutral") {
  tags$article(
    class = paste("ytab-metric-card", paste0("ytab-metric-", status)),
    tags$span(class = "ytab-metric-label", label),
    tags$strong(class = "ytab-metric-value", value),
    if (!is.null(caption) && nzchar(as.character(caption)))
      tags$span(class = "ytab-metric-caption", caption)
  )
}

ytab_warning_banner <- function(title, message = NULL) {
  tags$div(
    class = "ytab-warning-banner",
    tags$strong(title),
    if (!is.null(message)) tags$p(message)
  )
}

ytab_control_panel <- function(title, ...) {
  tags$aside(
    class = "ytab-control-panel",
    tags$h4(title),
    ...
  )
}

ytab_two_column_layout <- function(controls, main, controls_width = "320px") {
  tags$div(
    class = "ytab-two-column-layout",
    style = paste0("--ytab-controls-width:", controls_width, ";"),
    tags$div(class = "ytab-two-column-controls", controls),
    tags$div(class = "ytab-two-column-main", main)
  )
}

ytab_plot_card <- function(title, plot_ui, description = NULL, controls = NULL, downloads = NULL) {
  tags$article(
    class = "ytab-plot-card ytab-release-card",
    tags$div(
      class = "ytab-plot-card-header",
      tags$div(
        tags$h4(title),
        if (!is.null(description)) tags$p(class = "text-muted", description)
      ),
      if (!is.null(controls)) controls
    ),
    tags$div(class = "ytab-plot-card-body", plot_ui),
    if (!is.null(downloads)) tags$div(class = "ytab-actions", downloads)
  )
}

ytab_static_image_card <- function(title, src, alt = title, caption = NULL,
                                   max_height = "760px", download = NULL) {
  tags$article(
    class = "ytab-static-image-card ytab-release-card",
    tags$div(
      class = "ytab-plot-card-header",
      tags$div(
        tags$h4(title),
        if (!is.null(caption)) tags$p(class = "text-muted", caption)
      ),
      download
    ),
    tags$div(
      class = "ytab-static-image",
      tags$img(src = src, alt = alt, loading = "lazy",
               style = paste0("max-height:", max_height, ";"))
    )
  )
}

ytab_table_card <- function(title, table_ui, description = NULL, downloads = NULL) {
  tags$section(
    class = "ytab-table-card ytab-release-card",
    tags$div(
      class = "ytab-table-card-header",
      tags$div(
        tags$h4(title),
        if (!is.null(description)) tags$p(class = "text-muted", description)
      ),
      if (!is.null(downloads)) tags$div(class = "ytab-actions", downloads)
    ),
    table_ui
  )
}

ytab_download_card <- function(title, ...) {
  tags$section(
    class = "ytab-download-card ytab-release-card",
    tags$h4(title),
    tags$div(class = "ytab-actions", ...)
  )
}

ytab_technical_details <- function(summary = "Technical details", ...) {
  tags$details(
    class = "ytab-technical-details ytab-release-details",
    tags$summary(summary),
    ...
  )
}

ytab_empty_state <- function(title, message = NULL, action = NULL) {
  tags$div(
    class = "ytab-empty-state",
    tags$strong(title),
    if (!is.null(message)) tags$p(message),
    action
  )
}
