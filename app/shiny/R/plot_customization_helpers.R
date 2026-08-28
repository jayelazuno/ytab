ytab_plot_style_choices <- function() {
  c("Clean" = "clean", "Compact" = "compact", "Presentation" = "presentation")
}

ytab_plot_height_choices <- function() {
  c("Small" = "small", "Medium" = "medium", "Large" = "large", "Extra large" = "xlarge")
}

ytab_text_size_choices <- function() {
  c("Small" = "small", "Medium" = "medium", "Large" = "large")
}

ytab_plot_width_choices <- function() {
  c("Standard" = "standard", "Wide" = "wide", "Full panel" = "full")
}

ytab_label_mode_choices <- function() {
  c("Full labels" = "full", "Compact labels" = "compact", "Hide labels" = "hide")
}

ytab_grid_choices <- function() {
  c("Show grid" = "show", "Hide grid" = "hide")
}

ytab_label_angle_choices <- function() {
  c("0°" = "0", "30°" = "30", "45°" = "45", "90°" = "90")
}

ytab_static_image_mode_choices <- function() {
  c("Fit to panel" = "fit", "Natural size" = "natural")
}

ytab_plot_height_px <- function(choice, custom = NULL, default = "medium") {
  choice <- choice %||% default
  if (identical(choice, "custom") && !is.null(custom)) {
    value <- suppressWarnings(as.integer(custom))
    if (!is.na(value)) return(paste0(max(250L, min(value, 1400L)), "px"))
  }
  heights <- c(small = "380px", medium = "520px", large = "700px", xlarge = "860px")
  unname(heights[[choice]] %||% heights[[default]] %||% heights[["medium"]])
}

ytab_text_size_px <- function(choice) {
  sizes <- c(small = 12L, medium = 15L, large = 19L)
  unname(sizes[[choice %||% "medium"]] %||% 15L)
}

ytab_plot_customization_controls <- function(
    id,
    include_points = FALSE,
    include_bars = FALSE,
    include_heatmap = FALSE,
    include_labels = TRUE,
    include_grid = TRUE,
    include_value_labels = include_bars,
    default_height = "medium",
    default_width = "standard",
    default_bar_orientation = "vertical",
    default_label_angle = "30",
    default_show_value_labels = TRUE) {
  ns_id <- function(suffix) paste0(id, "_", suffix)
  controls <- list(
    selectInput(ns_id("plot_style"), "Plot style", ytab_plot_style_choices(), selected = "clean"),
    selectInput(ns_id("plot_height"), "Plot height", ytab_plot_height_choices(), selected = default_height),
    selectInput(ns_id("plot_width"), "Plot width", ytab_plot_width_choices(), selected = default_width),
    selectInput(ns_id("text_size"), "Text size", ytab_text_size_choices(), selected = "medium")
  )
  if (include_labels) {
    controls <- c(controls, list(
      selectInput(ns_id("label_mode"), "Label mode", ytab_label_mode_choices(), selected = "full"),
      selectInput(ns_id("label_angle"), "Label angle", ytab_label_angle_choices(), selected = default_label_angle)
    ))
  }
  if (include_grid)
    controls <- c(controls, list(selectInput(ns_id("grid"), "Grid", ytab_grid_choices(), selected = "show")))
  if (include_points) {
    controls <- c(controls, list(
      sliderInput(ns_id("point_size"), "Point size", min = 0.4, max = 4, value = 1.5, step = 0.1),
      sliderInput(ns_id("point_opacity"), "Point opacity", min = 0.1, max = 1, value = 0.55, step = 0.05)
    ))
  }
  if (include_bars)
    controls <- c(controls, list(selectInput(ns_id("bar_orientation"), "Bar orientation",
                                             c("Horizontal when labels are long" = "auto",
                                               "Horizontal" = "horizontal",
                                               "Vertical" = "vertical"),
                                             selected = default_bar_orientation)))
  if (include_value_labels)
    controls <- c(controls, list(checkboxInput(ns_id("show_value_labels"),
                                               "Show value labels", default_show_value_labels)))
  if (include_heatmap) {
    controls <- c(controls, list(
      selectInput(ns_id("heatmap_height"), "Heatmap height", ytab_plot_height_choices(), selected = "large"),
      sliderInput(ns_id("row_label_size"), "Row label size", min = 6, max = 18, value = 10, step = 1),
      selectInput(ns_id("column_label_angle"), "Column label angle", ytab_label_angle_choices(), selected = "45"),
      checkboxInput(ns_id("show_row_labels"), "Show row labels", TRUE),
      selectInput(ns_id("heatmap_palette"), "Palette",
                  c("Blue-white-red" = "blue_red", "Viridis" = "viridis",
                    "Magma" = "magma", "Gray" = "gray"),
                  selected = "blue_red")
    ))
  }
  tags$details(
    class = "ytab-plot-controls",
    open = "open",
    tags$summary("Plot display options"),
    tags$div(class = "ytab-control-grid", controls),
    tags$p(class = "ytab-display-only-note",
           "Display-only controls; underlying result tables and statistics are unchanged.")
  )
}

ytab_plot_customization_values <- function(input, id) {
  value <- function(suffix, fallback = NULL) input[[paste0(id, "_", suffix)]] %||% fallback
  list(
    style = value("plot_style", "clean"),
    height = ytab_plot_height_px(value("plot_height", "medium")),
    width = value("plot_width", "standard"),
    text_size = ytab_text_size_px(value("text_size", "medium")),
    label_mode = value("label_mode", "full"),
    label_angle = suppressWarnings(as.integer(value("label_angle", 30))),
    grid = value("grid", "show"),
    point_size = suppressWarnings(as.numeric(value("point_size", 1.5))),
    point_opacity = suppressWarnings(as.numeric(value("point_opacity", 0.55))),
    bar_orientation = value("bar_orientation", "vertical"),
    show_value_labels = isTRUE(value("show_value_labels", TRUE)),
    heatmap_height = ytab_plot_height_px(value("heatmap_height", "large"), default = "large"),
    row_label_size = suppressWarnings(as.numeric(value("row_label_size", 10))),
    column_label_angle = suppressWarnings(as.integer(value("column_label_angle", 45))),
    heatmap_palette = value("heatmap_palette", "blue_red"),
    show_row_labels = isTRUE(value("show_row_labels", TRUE))
  )
}

ytab_with_plot_display_options <- function(input, id, expr) {
  get <- function(suffix, fallback = NULL) input[[paste0(id, "_", suffix)]] %||% fallback
  old <- options(
    ytab.plot.style = get("plot_style", "clean"),
    ytab.plot.text_size = get("text_size", "medium"),
    ytab.plot.label_angle = get("label_angle", "30"),
    ytab.plot.label_mode = get("label_mode", "full"),
    ytab.plot.grid = get("grid", "show"),
    ytab.plot.bar_orientation = get("bar_orientation", "vertical"),
    ytab.plot.show_value_labels = isTRUE(get("show_value_labels", TRUE))
  )
  on.exit(options(old), add = TRUE)
  force(expr)
}
