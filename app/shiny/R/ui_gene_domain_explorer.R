gene_domain_explorer_ui <- function() {
  ytab_two_column_layout(
    controls = ytab_control_panel(
        "Insertion Explorer controls",
        uiOutput("gene_domain_project_indicator"),
        textInput("gene_domain_query", "Gene name or ID", "",
                  placeholder = "systematic name, standard name, locus tag, alias"),
        actionButton("gene_domain_search", "Search gene", class = "btn-primary"),
        tags$hr(),
        selectInput(
          "gene_domain_track_source", "Track source",
          choices = c("Raw CreateHitFile tracks" = "raw",
                      "Combined parent hit file" = "combined_parent",
                      "All available tracks" = "all"),
          selected = "raw"
        ),
        selectInput("gene_domain_track_preset", "Track preset",
                    choices = c("All tracks" = "all",
                                "Parents only" = "parents",
                                "Treated only" = "treated",
                                "Matched pairs" = "matched_pairs",
                                "Custom" = "custom"),
                    selected = "all"),
        tags$p(class = "text-muted ytab-control-help",
               "Choose which insertion tracks are displayed for the selected gene."),
        uiOutput("gene_domain_preset_message"),
        selectizeInput("gene_domain_samples", "Tracks to display",
                       choices = character(), multiple = TRUE,
                       options = list(placeholder = "All tracks by default")),
        numericInput("gene_domain_flank_bp", "Flank bp", value = 1000,
                     min = 0, max = 100000, step = 100),
        checkboxInput("gene_domain_show_domains", "Show domains when available", TRUE),
        checkboxInput("gene_domain_show_direction", "Show gene direction", TRUE),
        numericInput("gene_domain_width_px", "Image width", value = 1800,
                     min = 900, max = 4000, step = 100),
        selectInput("gene_domain_display_height", "Preview height",
                    choices = ytab_plot_height_choices(),
                    selected = "large"),
        selectInput("gene_domain_label_mode", "Label mode",
                    choices = c("Full sample names" = "full",
                                "Compact display labels" = "compact"),
                    selected = "full"),
        checkboxInput("gene_domain_show_site_counts", "Show site counts in labels", TRUE),
        actionButton("gene_domain_generate", "Generate figure", class = "btn-primary")
    ),
    main = panel_card(
        "Gene & Domain Insertion Explorer",
        tags$p(class = "ytab-stage-purpose",
               "Search a gene, choose existing insertion tracks, and plot insertions across the gene region."),
        uiOutput("gene_domain_status"),
        uiOutput("gene_domain_gene_summary"),
        DT::DTOutput("gene_domain_candidates"),
        uiOutput("gene_domain_domain_status"),
        uiOutput("gene_domain_figure_ui"),
        tags$details(
          class = "ytab-more-options",
          tags$summary("Insertion table"),
          DT::DTOutput("gene_domain_insertions")
        ),
        tags$details(
          class = "ytab-technical-details",
          tags$summary("Technical details"),
          verbatimTextOutput("gene_domain_technical")
        )
    )
  )
}
