genome_browser_ui <- function() {
  ytab_two_column_layout(
    controls = ytab_control_panel(
      "Genome browser controls",
      uiOutput("genome_browser_status"),
      textInput("genome_browser_locus", "Locus", "",
                placeholder = "CP048230.1:1-50000 or chromosome:start-end"),
      actionButton("genome_browser_go", "Go to locus", class = "btn-primary"),
      tags$div(
        class = "ytab-actions",
        actionButton("genome_browser_pan_left", "\u2190 Pan left"),
        actionButton("genome_browser_pan_right", "Pan right \u2192")
      ),
      tags$hr(),
      selectInput("genome_browser_track_preset", "Track preset",
                  choices = c("All tracks" = "all", "Custom" = "custom"),
                  selected = "all"),
      tags$p(class = "text-muted ytab-control-help",
             "Choose which insertion tracks are loaded into IGV."),
      uiOutput("genome_browser_preset_message"),
      selectizeInput("genome_browser_tracks", "Tracks to display",
                     choices = character(), multiple = TRUE,
                     options = list(placeholder = "Tracks discovered from CreateHitFile outputs")),
      tags$div(class = "ytab-actions",
               actionButton("genome_browser_reload", "Reload browser"),
               actionButton("genome_browser_clear_tracks", "Clear tracks")),
      tags$details(
        class = "ytab-technical-details",
        tags$summary("Technical details"),
        verbatimTextOutput("genome_browser_technical")
      )
    ),
    main = panel_card(
      "Genome Browser",
      tags$p(class = "ytab-stage-purpose",
             "Embedded IGV view of the current project reference and insertion tracks."),
      tags$div(id = "ytab_igv_warning", class = "ytab-warning-banner", style = "display:none;"),
      uiOutput("genome_browser_lane_key"),
      tags$div(id = "ytab_igv_browser", class = "ytab-igv-browser", tabindex = "0"),
      tags$p(class = "text-muted ytab-control-help",
             "Tip: click the browser, then use ← and → to pan across nearby genes."),
      tags$details(
        class = "ytab-more-options",
        tags$summary("Tracks loaded"),
        DT::DTOutput("genome_browser_track_table")
      )
    )
  )
}
