# Gene & Domain Insertions: everything that answers "what happened at this gene", for one species and across them.
#
# The four inner tabs were previously split between this module and a separate "Across Species" subtab, whose own
# Single Species View showed the same gene from the same projects as this explorer but as tables only. They are one
# question asked at four widths, so they are one navset now, and the explorer itself is the Single Species View.

gene_domain_single_species_ui <- function() {
  ytab_two_column_layout(
    controls = ytab_control_panel(
        "Single species controls",
        uiOutput("gene_domain_project_indicator"),
        textInput("gene_domain_query", "Gene name or ID", "",
                  placeholder = "systematic name, standard name, locus tag, alias"),
        actionButton("gene_domain_search", "Search gene", class = "btn-primary"),
        helpText("The gene and its insertion figure appear as soon as the gene resolves."),
        tags$hr(),
        selectInput("gene_domain_track_style", "Figure",
                    choices = c("Insertion sites per track" = "sites",
                                "Reads per insertion site, along the coding sequence" = "reads"),
                    selected = "sites"),
        uiOutput("gene_domain_style_note"),
        tags$hr(),
        tags$p(tags$b("Compare with other datasets")),
        ytab_reference_projects_ui("gene_domain", display_modes = FALSE),
        tags$details(class = "ytab-more-options", tags$summary("Tracks"),
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
          uiOutput("gene_domain_preset_message"),
          selectizeInput("gene_domain_samples", "Tracks to display",
                         choices = character(), multiple = TRUE,
                         options = list(placeholder = "All tracks by default"))),
        tags$details(class = "ytab-more-options", tags$summary("Display"),
          numericInput("gene_domain_flank_bp", "Flank bp", value = 1000,
                       min = 0, max = 100000, step = 100),
          checkboxInput("gene_domain_show_domains", "Show domains when available", TRUE),
          checkboxInput("gene_domain_show_direction", "Show gene direction", TRUE),
          numericInput("gene_domain_width_px", "Image width", value = 1800,
                       min = 900, max = 4000, step = 100),
          selectInput("gene_domain_display_height", "Preview height",
                      choices = c("Fit to the window" = "fit", ytab_plot_height_choices()),
                      selected = "fit"),
          selectInput("gene_domain_label_mode", "Label mode",
                      choices = c("Full sample names" = "full",
                                  "Compact display labels" = "compact"),
                      selected = "full"),
          checkboxInput("gene_domain_show_site_counts", "Show site counts in labels", TRUE))
    ),
    main = panel_card(
        "Single species view",
        tags$p(class = "ytab-stage-purpose",
               "Search a gene to see where insertions fall across it in this dataset, and optionally in other datasets."),
        uiOutput("gene_domain_status"),
        uiOutput("gene_domain_gene_summary"),
        uiOutput("gene_domain_domain_status"),
        # The figure's own element is built once, here, and only shown or hidden. Rebuilding an imageOutput inside
        # a renderUI races its own value: the server pushes the new image at the element being replaced, so the
        # fresh one stays empty and the figure appears to need a second search to change.
        uiOutput("gene_domain_figure_message"),
        conditionalPanel(
          condition = "output.gene_domain_has_figure === true",
          tags$article(
            class = "ytab-static-image-card ytab-release-card",
            tags$div(class = "ytab-plot-card-header",
              tags$div(tags$h4("Gene/domain insertion figure"),
                       uiOutput("gene_domain_figure_caption")),
              tags$div(class = "ytab-actions",
                       downloadButton("gene_domain_download_png", "Download PNG"),
                       downloadButton("gene_domain_download_table", "Download table"))),
            tags$div(class = "ytab-static-image ytab-static-image-output",
                     imageOutput("gene_domain_figure", width = "100%", height = "auto")))),
        tags$details(
          class = "ytab-more-options",
          tags$summary("Matching genes"),
          DT::DTOutput("gene_domain_candidates")
        ),
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

gene_domain_explorer_ui <- function() {
  navset_tab(
    id = "comparative_tabs",
    nav_panel("Single Species View", value = "single_species", gene_domain_single_species_ui()),
    nav_panel("Comparative View", value = "comparative", comparative_comparison_ui()),
    nav_panel("Gene Group Analysis", value = "gene_group", comparative_gene_group_ui()),
    nav_panel("Essentiality across species", value = "essentiality_species", cross_species_ui())
  )
}
