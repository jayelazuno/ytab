
# Species names are set in italics wherever they are shown, as a species name is written. Two spellings, because
# the figures are drawn once for the browser and once for download: plotly renders the HTML tags in a label,
# ggsave does not and needs plotmath instead.
ytab_species_italic <- function(species, html = TRUE) {
  species <- as.character(species); species[is.na(species)] <- ""
  if (html) paste0("<i>", htmltools::htmlEscape(species), "</i>")
  else parse(text = sprintf('italic("%s")', gsub('"', "", species, fixed = TRUE)))
}

# A column read back from a result table, as text with the blanks blank. `ytab_read_analysis_tsv` turns an empty
# field into NA, and `nzchar(NA)` is TRUE, so a missing call would otherwise pass a "has a value" test and then
# turn an all()/any() over it into NA - which is the "missing value where TRUE/FALSE needed" an if() then raises.
ytab_text_column <- function(data, name) {
  value <- if (name %in% names(data)) as.character(data[[name]]) else rep("", nrow(data))
  value[is.na(value)] <- ""
  value
}

comparative_controls_card <- function(title, ...) {
  tags$aside(class = "ytab-comparative-controls ytab-control-panel", panel_card(title, ...))
}

comparative_result_card <- function(title, ...) {
  tags$section(class = "ytab-comparative-results ytab-release-card", panel_card(title, ...))
}

comparative_layout <- function(controls, results) {
  tags$div(
    class = "row ytab-comparative-layout",
    tags$div(class = "col-sm-4 col-md-3", controls),
    tags$div(class = "col-sm-8 col-md-9", results)
  )
}

# The former Single Species View lived here as three tables of the same gene in one project. The Gene & Domain
# Insertions explorer answers that question with the insertion figure as well, so it is the Single Species View now
# (ui_gene_domain_explorer.R) and this one is gone rather than kept as a second way in.

comparative_species_checkbox_ui <- function() {
  tagList(
    checkboxGroupInput("comparative_species_set", "Species",
                       choices = character(), selected = character())
  )
}

comparative_comparison_ui <- function() {
  comparative_layout(
    comparative_controls_card(
      "Comparative controls",
      textInput("comparative_gene", "Gene name or ID", "",
                placeholder = "any species' name or systematic id"),
      actionButton("comparative_compare", "Compare across species", class = "btn-primary"),
      helpText("Each species contributes one library set: the densest that can support an essentiality call. ",
               "A species whose libraries are too sparse is named as such rather than given a call."),
      # The same two readings the single species view offers, under the same names. Switching between them only
      # redraws: the profile table already holds the reads at every site, so neither reading needs another search.
      radioButtons("comparative_track_style", "Show",
                   choices = c("Reads per insertion site" = "reads", "Insertion sites" = "sites"),
                   selected = "reads"),
      tags$details(class = "ytab-more-options", tags$summary("Settings"),
        numericInput("comparative_min_density", "Insertion sites per coding kb a species must reach",
                     value = 18, min = 1, max = 200, step = 1),
        helpText("About 36 for the label-free HMM, about 18 where curated labels are available ",
                 "(docs/methods/00_recommendations.md).")),
      ytab_plot_downloads_ui("comparative", c(profiles = "Download profiles (TSV)",
                                              representatives = "Download representatives (TSV)"))
    ),
    comparative_result_card(
      "Comparative View",
      tags$p(class = "ytab-stage-purpose",
             "One gene across species: where insertions fall in it, and what each species' libraries call it."),
      uiOutput("comparative_state"),
      ytab_plot_card("Insertions across species",
                     plotly::plotlyOutput("comparative_plot", height = "620px"),
                     description = uiOutput("comparative_plot_caption")),
      tags$details(class = "ytab-more-options", tags$summary("Which libraries represent each species"),
                   DT::DTOutput("comparative_representatives")),
      tags$details(class = "ytab-more-options", tags$summary("Methods"),
                   tags$p(tags$b("Representative. "),
                          "Among a species' essentiality runs, the one whose libraries reach the highest pooled ",
                          "insertion density (sites per coding kb, median gene), provided it reaches the density ",
                          "a call needs. Density is what decides whether a call means anything, and the deepest ",
                          "set is usually the unselected wild-type pool, so the rule needs no list of strain names."),
                   tags$p(tags$b("Finding the gene. "),
                          "The query is matched against every represented species' own gene identifiers, names ",
                          "and recorded aliases at once, rather than resolved in one species and carried to the ",
                          "others through the ortholog tables. A gene that is present in a species and named as ",
                          "asked for is then found there even when no table connects it to S. cerevisiae, which ",
                          "is the usual situation for a family that expanded in that lineage. Each panel says ",
                          "which of the three routes reached it: a curated ortholog table (strongest, and the ",
                          "best coverage for the distant species), a shared OrthoFinder gene family (relatedness ",
                          "rather than correspondence, and the only route able to report a family that a ",
                          "one-to-one table cannot express), or the gene's name alone (weakest). The routes are ",
                          "unioned, so a gap in any one source costs an edge rather than a wrong answer."),
                   tags$p(tags$b("Orthologs. "),
                          "Through the S. cerevisiae hub, from each project's own ortholog table (reference: ",
                          "orthology_file), as whole families rather than one-to-one pairs. Where a family has ",
                          "expanded in a species, every member gets its own panel, numbered within the family: the ",
                          "expansion is usually the point of the comparison, and keeping only unambiguous pairs ",
                          "would report such a species as having no counterpart when it has several. A species ",
                          "with no counterpart at all keeps an empty panel saying so, because a species dropped ",
                          "from the figure and a species missing the gene look identical."),
                   tags$p(tags$b("Profile. "),
                          "Reads at each insertion site over the coding sequence and flanks of 25% of its length, ",
                          "summed over the representative's libraries. One panel per species, named for the species ",
                          "and its call, each scaled within itself, so what compares between species of different ",
                          "depth is where the insertions fall. Reads per insertion site puts the read count on a log ",
                          "axis; insertion sites draws every site as an equal tick with the reads at it in its width, ",
                          "as the single species view does, and reports n per species. n counts the sites inside the ",
                          "coding sequence only: the flank sites are drawn, because that contrast is what the figure ",
                          "is for, but a gene's density is a property of its coding sequence."),
                   tags$p("Command: ", tags$code("ytab gene-profiles"), "."))
    )
  )
}

comparative_gene_group_ui <- function() {
  comparative_layout(
    comparative_controls_card(
      "Gene group controls",
      fileInput("comparative_group_file", "Upload CSV gene group",
                accept = c(".csv", "text/csv")),
      textAreaInput("comparative_group_text", "Paste gene list", "",
                    rows = 6, placeholder = "One gene ID or name per line"),
      actionButton("comparative_group_analyze", "Analyze group",
                   class = "btn-primary"),
      helpText("Every species contributes the densest library set that can support an essentiality call, ",
               "the same representative Comparative View uses. Each gene is matched in every species at once."),
      radioButtons("comparative_group_track_style", "Show",
                   choices = c("Reads per insertion site" = "reads", "Insertion sites" = "sites"),
                   selected = "reads"),
      tags$details(class = "ytab-more-options", tags$summary("Settings"),
        numericInput("comparative_group_min_density", "Insertion sites per coding kb a species must reach",
                     value = 18, min = 1, max = 200, step = 1)),
      ytab_plot_downloads_ui("comparative_group", c(calls = "Download calls (TSV)",
                                                    profiles = "Download profiles (TSV)"))
    ),
    comparative_result_card(
      "Gene Group Analysis",
      tags$p(class = "ytab-stage-purpose",
             "A list of genes, and what each species' representative libraries call each of them."),
      uiOutput("comparative_group_state"),
      uiOutput("comparative_group_summary"),
      ytab_plot_card("Insertions across species, gene by gene",
                     plotly::plotlyOutput("comparative_group_profiles", height = "700px"),
                     description = uiOutput("comparative_group_profiles_caption")),
      tags$details(class = "ytab-more-options", tags$summary("Summary of the calls"),
                   ytab_plot_card("Every gene in the group",
                                  plotly::plotlyOutput("comparative_group_plot", height = "560px"),
                                  description = uiOutput("comparative_group_plot_caption"))),
      tags$details(class = "ytab-more-options", tags$summary("Genes that matched nothing"),
                   DT::DTOutput("comparative_group_unmapped")),
      tags$details(class = "ytab-more-options", tags$summary("Which libraries represent each species"),
                   DT::DTOutput("comparative_group_availability")),
      tags$details(class = "ytab-more-options", tags$summary("Result rows"),
                   DT::DTOutput("comparative_group_results"))
    )
  )
}

# Rows whose names or IDs contain the query (any case), or whose ID is one of `ids` (the gene's orthologs, matched
# exactly; a GFF prefix such as "gene:" is ignored).



comparative_availability_display <- function(availability) {
  if (!is.data.frame(availability) || !nrow(availability)) return(data.frame())
  data.frame(
    Species = as.character(availability$label %||% availability$species),
    `YTAB project outputs` = as.character(availability$status %||% ""),
    Projects = as.integer(availability$project_count %||% 0L),
    `Essentiality results` = as.integer(availability$classifier_projects %||% 0L),
    `Fitness results` = as.integer(availability$fitness_projects %||% 0L),
    `Insertion summaries` = as.integer(availability$summary_projects %||% 0L),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Essentiality call -> colour, shared with the Essentiality module so a call means one thing everywhere in YTAB.

# The first column present out of `names`, as character; "" when none is.

# Rows reduced to one score and one call per species (and per gene, for a group), which is what both plots draw.
# Essentiality rows carry hmm_call/hmm_score; fitness rows a call and a log2 fold change; insertion summaries
# neither, and those are shown as presence only.

# Stack result tables that do not share a column set. Species differ in which models their run produced - a project
# with curated labels carries the hybrid's score and its agreement with the HMM, one without carries neither - so
# rbind() refuses them. The union is kept, with NA where a species has no such column.

comparative_server <- function(input, output, session, repo_root, active = function() list(),
                               python_bin = locate_python_executable) {
  keep_or <- function(value, fallback) if (length(value) && nzchar(value[[1]])) value else fallback
  species_manifest <- reactive(load_comparative_species_manifest(repo_root))
  # Both tabs now resolve genes through `ytab gene-profiles`, which reads the ortholog tables, the orthogroups
  # and each species' own annotation itself. The orthogroup and id-alias readers that used to be wired in here
  # went with it, and with them a C. glabrata-specific CAGL-to-GWK60 lookup: translating one species' gene ids
  # by name was the kind of per-species path that cannot be carried to a species nobody has added yet.

  observe({
    species <- species_manifest()
    enabled <- species[species$enabled, , drop = FALSE]
    choices <- as.list(setNames(enabled$species, enabled$label))
    # default: the open project's species, else the first enabled species
    current <- tolower(as.character(tryCatch(active()$species, error = function(e) "") %||% ""))
    first <- if (current %in% enabled$species) current else if (length(choices)) unname(unlist(choices, use.names = FALSE))[[1]] else ""
    updateSelectInput(session, "comparative_group_species",
                      choices = choices,
                      selected = keep_or(isolate(input$comparative_group_species), first))
  })

  # rendered with its choices: an empty group rendered after the observer's update would otherwise replace them
  output$comparative_species_checks <- renderUI({
    species <- species_manifest(); enabled <- species[species$enabled, , drop = FALSE]
    checkboxGroupInput("comparative_species_set", "Species", choices = stats::setNames(enabled$species, enabled$label),
                       selected = isolate(input$comparative_species_set) %||% enabled$species)
  })

  # One gene across species, computed by `ytab gene-profiles`: it picks each species' representative library set,
  # maps the gene through the ortholog families and reads the insertion profile. The app only displays the
  # result, as it does for every other analysis.
  comparative_run <- eventReactive(input$comparative_compare, {
    query <- trimws(input$comparative_gene %||% "")
    if (!nzchar(query)) return(list(status = "no_query"))
    out <- file.path(repo_root, "output", "gene_profiles", gsub("[^A-Za-z0-9_.-]", "_", query))
    args <- c(file.path(repo_root, "scripts", "local", "ytab_gene_profiles.py"),
              "--gene", query, "--output", out,
              "--min-density", as.character(input$comparative_min_density %||% 18), "--json")
    result <- tryCatch(run_process_sync(python_bin(), args, wd = repo_root),
                       error = function(e) list(status = 1L, stdout = "", stderr = conditionMessage(e)))
    if (!identical(as.integer(result$status), 0L))
      return(list(status = "failure", error = paste(c(result$stderr, result$stdout), collapse = "\n")))
    # The manifest states the flank the CLI actually drew, so the axis is labelled and sized from the run rather
    # than from a constant repeated here that could drift away from it.
    manifest <- tryCatch(jsonlite::fromJSON(file.path(out, "manifest.json")), error = function(e) list())
    list(status = "success", query = query, dir = out,
         flank = manifest$flank_fraction %||% NA_real_,
         representatives = ytab_read_analysis_tsv(file.path(out, "representatives.tsv")),
         profiles = ytab_read_analysis_tsv(file.path(out, "profiles.tsv")))
  }, ignoreInit = TRUE)

  output$comparative_state <- renderUI({
    if ((input$comparative_compare %||% 0L) == 0L)
      return(ytab_empty_state("Search a gene to compare it across species",
        "Each species contributes the densest library set that can support an essentiality call; species whose libraries cannot are named rather than left blank."))
    run <- comparative_run()
    if (identical(run$status, "no_query")) return(tags$p(class = "ytab-warning", "Enter a gene name or ID first."))
    if (identical(run$status, "failure"))
      return(ytab_result_card("Could not compare this gene across species", status = "warning",
                              tags$p(run$error %||% "")))
    NULL
  })

  # The panels a comparison should show: one per counterpart the gene has in each represented species, and one
  # for a represented species that has no counterpart at all. A species that cannot support a call in the first
  # place is not a panel - it is reported in the banner and the table - so the figure stays about the gene.
  comparative_panels <- function(reps, html = TRUE) {
    kept <- reps[as.character(reps$represented) %in% c("1", "TRUE", "True"), , drop = FALSE]
    if (!nrow(kept)) return(data.frame())
    n <- suppressWarnings(as.integer(kept$n_orthologs)); n[is.na(n)] <- 0L
    index <- suppressWarnings(as.integer(kept$ortholog_index)); index[is.na(index)] <- 0L
    gene <- ytab_text_column(kept, "gene")
    gene_name <- ytab_text_column(kept, "gene_name")
    shown_gene <- ifelse(nzchar(gene_name) & gene_name != gene, sprintf("%s (%s)", gene_name, gene), gene)
    call <- ytab_text_column(kept, "call")
    # A family that expanded is numbered in its own panels, so nine C. glabrata panels for one S. cerevisiae gene
    # read as one family rather than nine unrelated results.
    family <- ifelse(n > 1L, sprintf(" (%d/%d)", index, n), "")
    # How the gene was reached is part of what the panel claims. A gene found through an ortholog table is
    # asserted to correspond to the query; one found only because it carries the same name is not, and saying so
    # on the panel is the difference between reporting evidence and implying it.
    matched <- ytab_text_column(kept, "matched_by")
    route <- ifelse(matched == "orthogroup", " \u00b7 same gene family",
             ifelse(matched == "name", " \u00b7 name match only", ""))
    # The call goes on its own line. In one line a systematic id and a family number together overran the width
    # of a column panel, and it was the call - the part being compared - that was cut off.
    # Italic only in the browser's copy. plotly renders the tag in a strip annotation; ggsave would print it.
    shown_species <- if (html) paste0("<i>", kept$label, "</i>") else as.character(kept$label)
    label <- ifelse(!nzchar(gene),
                    sprintf("%s\n%s", shown_species, ifelse(nzchar(kept$status), kept$status, "no ortholog")),
                    sprintf("%s \u00b7 %s%s\n%s%s", shown_species, shown_gene, family,
                            ifelse(nzchar(call), call, "no call"), route))
    out <- data.frame(species = as.character(kept$species), gene = gene, label = as.character(kept$label),
                      status = as.character(kept$status), panel_label = label, matched_by = matched,
                      has_gene = nzchar(gene), index = index, stringsAsFactors = FALSE)
    out <- out[order(out$label, out$index), , drop = FALSE]
    out$key <- paste(out$species, out$gene, sep = "\u0001")
    out$panel <- factor(out$panel_label, levels = unique(out$panel_label))
    out
  }

  comparative_figure <- function(html = TRUE) {
    run <- comparative_run(); req(identical(run$status, "success"))
    profiles <- run$profiles; reps <- run$representatives
    panels <- comparative_panels(reps, html)
    validate(need(nrow(panels) > 0, "No species has libraries dense enough to place this gene."))
    style <- input$comparative_track_style %||% "reads"
    levels_ <- levels(panels$panel)

    # Every represented species keeps its panel, including one with no counterpart: an absent species and an
    # absent gene look the same once the species is dropped, and they are not the same answer.
    if (nrow(profiles)) {
      profiles$key <- paste(profiles$species, profiles$gene, sep = "\u0001")
      profiles <- profiles[profiles$key %in% panels$key, , drop = FALSE]
      profiles$panel <- factor(panels$panel_label[match(profiles$key, panels$key)], levels = levels_)
    }
    drawn <- if (nrow(profiles)) profiles else profiles[0, , drop = FALSE]
    empty <- panels[!panels$key %in% unique(drawn$key), , drop = FALSE]

    # The panels are laid out in columns once there are more than a few, so that an expanded family fits the
    # window instead of running off the bottom of it.
    ncol <- if (nrow(panels) > 15L) 4L else if (nrow(panels) > 8L) 3L else if (nrow(panels) > 4L) 2L else 1L
    # The flank the CLI drew, reported in its manifest, so the empty panels span the same axis as the drawn ones
    # rather than a width guessed here.
    flank <- suppressWarnings(as.numeric(run$flank %||% NA)); if (is.na(flank)) flank <- .25
    x_range <- range(c(-flank, 1 + flank, if (nrow(drawn)) drawn$cds_fraction))
    top <- if (identical(style, "reads"))
      stats::setNames(vapply(split(drawn$reads, drawn$panel, drop = FALSE),
                             function(x) if (length(x)) max(x) else 10, numeric(1)), levels_) else
      stats::setNames(rep(1.08, length(levels_)), levels_)
    band <- data.frame(panel = factor(levels_, levels = levels_),
                       ymin = if (identical(style, "reads")) 1 else 0,
                       ymax = as.numeric(top[levels_]), stringsAsFactors = FALSE)

    plot <- ggplot2::ggplot() +
      ggplot2::geom_rect(data = band, ggplot2::aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax),
                         fill = "grey92", inherit.aes = FALSE)
    if (nrow(drawn)) {
      drawn$call <- panels$panel_label[match(drawn$key, panels$key)]
      drawn$hover <- sprintf("%s\n%s\n%s reads at %.0f%% along the coding sequence",
                             drawn$label, drawn$gene, drawn$reads, 100 * drawn$cds_fraction)
      if (identical(style, "reads")) {
        # Height is the read count on a log axis, so the axis states what the tick means. The stem starts at one
        # read, which is the floor of the measurement: a site exists because at least one read landed there.
        plot <- plot +
          ggplot2::geom_segment(data = drawn,
                                ggplot2::aes(x = cds_fraction, xend = cds_fraction, y = 1, yend = reads,
                                             colour = label, text = hover), linewidth = .3) +
          ggplot2::geom_point(data = drawn, ggplot2::aes(x = cds_fraction, y = reads, colour = label, text = hover),
                              size = .5)
      } else {
        # Every site an equal tick, as the single species view draws them, with the reads at it in the line width.
        # The width is binned rather than continuous: plotly gives each distinct line width its own trace, so a
        # width per read count turned a three-species figure into nearly two hundred traces.
        drawn$width_bin <- cut(drawn$reads, breaks = c(-Inf, 10, 100, Inf), labels = c(".35", ".7", "1.2"))
        plot <- plot +
          ggplot2::geom_segment(data = drawn,
                                ggplot2::aes(x = cds_fraction, xend = cds_fraction, y = 0, yend = 1,
                                             colour = label, text = hover,
                                             linewidth = as.numeric(as.character(width_bin))))
      }
    }
    if (identical(style, "sites")) {
      # How many sites the panel is reporting, counted the way the single species view counts them: inside the
      # coding sequence only. The flank sites are drawn, because that is the contrast the figure is for, but they
      # are not part of n. `in_cds` arrives from the TSV as the strings True/False, so it is compared as text -
      # read as a logical it would be NA, and every site would fall into the count.
      counts <- if (nrow(drawn))
        do.call(rbind, lapply(split(drawn, drawn$panel, drop = TRUE), function(d) data.frame(
          panel = d$panel[[1]], n = sum(tolower(as.character(d$in_cds)) %in% c("true", "1", "yes")),
          stringsAsFactors = FALSE))) else NULL
      if (!is.null(counts) && nrow(counts))
        plot <- plot + ggplot2::geom_text(data = counts,
                                          ggplot2::aes(x = x_range[2], y = 1.045, label = sprintf("n = %d in CDS", n)),
                                          hjust = 1, vjust = 0, size = 3.4, fontface = "bold",
                                          colour = "grey20", inherit.aes = FALSE)
      plot <- plot + ggplot2::scale_linewidth_identity() +
        ggplot2::scale_y_continuous(breaks = NULL, limits = c(0, 1.14)) +
        ggplot2::labs(y = "Insertion sites")
    } else {
      plot <- plot + ggplot2::scale_y_log10() + ggplot2::labs(y = "Reads per insertion site (log scale)")
    }
    if (nrow(empty))
      # Said in the panel rather than left blank: an empty panel on its own is read as a gene with no insertions,
      # which is the opposite of what a missing ortholog means.
      plot <- plot + ggplot2::geom_text(
        data = data.frame(panel = factor(empty$panel_label, levels = levels_),
                          note = ifelse(empty$has_gene, empty$status, "no ortholog"), stringsAsFactors = FALSE),
        ggplot2::aes(x = mean(x_range), y = if (identical(style, "sites")) .5 else 3, label = note),
        size = 3.6, fontface = "italic", colour = "grey45", inherit.aes = FALSE)
    plot +
      ggplot2::facet_wrap(~panel, ncol = ncol, scales = "free_y", strip.position = "top", drop = FALSE) +
      ggplot2::scale_colour_manual(values = ytab_species_colour(unique(panels$label)),
                                   breaks = unique(panels$label), name = NULL) +
      ggplot2::scale_x_continuous(breaks = c(0, .5, 1), labels = c("start", "50%", "stop"),
                                  limits = x_range) +
      ggplot2::labs(x = sprintf("Position along the coding sequence (shaded; flanks are %g%% of its length)",
                                100 * flank)) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(legend.position = "none",
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face = "bold", hjust = 0, size = 10, lineheight = 1.15),
                     # The log axis runs to four digits, which left the title sitting on the tick labels.
                     axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
                     axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)))
  }

  output$comparative_plot <- plotly::renderPlotly(ytab_ggplotly(comparative_figure(TRUE), tooltip = "text"))

  output$comparative_plot_caption <- renderUI({
    run <- comparative_run(); req(identical(run$status, "success"))
    reps <- run$representatives
    represented <- as.character(reps$represented) %in% c("1", "TRUE", "True")
    profiled <- reps[reps$status == "profiled", , drop = FALSE]
    # Species, not rows: an expanded family contributes several rows for one species, and counting rows would
    # report nine C. glabrata panels as nine species.
    species_profiled <- length(unique(profiled$species))
    families <- profiled[suppressWarnings(as.integer(profiled$n_orthologs)) > 1L, , drop = FALSE]
    absent <- reps[represented & !nzchar(as.character(reps$gene %||% "")), , drop = FALSE]
    unavailable <- reps[!represented, , drop = FALSE]
    per_species <- function(d) {
      parts <- lapply(split(d, d$species), function(x) sprintf(
        "%s: %s (%s, %s sites/coding kb)", ytab_species_italic(x$label[[1]]),
        paste(unique(ifelse(nzchar(x$call), x$call, "no call")), collapse = "/"),
        x$project[[1]], x$sites_per_coding_kb[[1]]))
      paste(unlist(parts), collapse = " \u00b7 ")
    }
    tagList(
      tags$p(class = "text-muted", HTML(sprintf(
        "%s in %d species%s. %s Each panel is named for its species, its counterpart of the gene and that gene's call, and is scaled within itself, so what compares between panels is where the insertions fall, not how many.",
        htmltools::htmlEscape(run$query), species_profiled,
        if (nrow(families)) sprintf(" (%s)", paste(unique(sprintf("%s has %s counterparts",
          ytab_species_italic(families$label), families$n_orthologs)), collapse = "; ")) else "",
        if (identical(input$comparative_track_style %||% "reads", "reads"))
          "The axis is reads at each insertion site, on a log scale." else
          "Every insertion site is one tick, with the reads at it in its width."))),
      if (nrow(profiled)) tags$p(class = "text-muted", HTML(per_species(profiled))),
      if (nrow(absent))
        tags$p(class = "text-muted", HTML(sprintf("No counterpart of this gene in %s; shown as an empty panel.",
                                             paste(unique(ytab_species_italic(absent$label)), collapse = ", ")))),
      {
        route <- as.character(profiled$matched_by %||% "")
        notes <- c(
          if (any(route == "orthogroup")) sprintf(
            "Some panels (%s) are in the same gene family as the query rather than named as its counterpart by an ortholog table: relatedness, not correspondence.",
            paste(unique(ytab_species_italic(profiled$label[route == "orthogroup"])), collapse = ", ")),
          if (any(route == "name")) sprintf(
            "Some panels (%s) carry this gene's name with nothing connecting them to the rest, so they are shown on the name alone; a lineage-specific expansion often keeps the name without having a counterpart.",
            paste(unique(ytab_species_italic(profiled$label[route == "name"])), collapse = ", ")))
        if (length(notes)) tags$p(class = "text-muted", HTML(paste(notes, collapse = " ")))
      },
      if (nrow(unavailable))
        ytab_warning_banner("Not every species could be compared",
          HTML(paste(sprintf("%s: %s.", ytab_species_italic(unavailable$label),
                             htmltools::htmlEscape(unavailable$status)), collapse = " "))))
  })

  # Species columns in the tables are italicised too, which means the column is HTML and the table must be told
  # not to escape it; every other column stays escaped.
  ytab_italicise_species_column <- function(data, column = "label") {
    if (is.data.frame(data) && column %in% names(data)) data[[column]] <- ytab_species_italic(data[[column]])
    data
  }

  output$comparative_representatives <- DT::renderDT({
    run <- comparative_run(); req(identical(run$status, "success"))
    reps <- run$representatives
    if (!nrow(reps)) return(NULL)
    shown <- data.frame(Species = ytab_species_italic(reps$label),
                        Gene = reps$gene, Name = reps$gene_name,
                        `Matched by` = reps$matched_by, Call = reps$call,
                        `P(ES)` = reps$P_ES, `Sites per coding kb` = reps$sites_per_coding_kb,
                        Project = reps$project, Libraries = reps$n_libraries, Status = reps$status,
                        check.names = FALSE, stringsAsFactors = FALSE)
    DT::datatable(shown, rownames = FALSE, selection = "none", escape = FALSE,
                  options = list(dom = "t", scrollX = TRUE))
  })

  comparative_download_name <- function(suffix)
    paste0("gene_profiles_", gsub("[^A-Za-z0-9_.-]", "_", comparative_run()$query %||% "gene"), suffix)
  for (format in c("png", "pdf")) local({
    ext <- format
    output[[paste0("comparative_download_", ext)]] <- downloadHandler(
      filename = function() comparative_download_name(paste0(".", ext)),
      content = function(file) ytab_save_ggplot(file, comparative_figure(FALSE), width = 10, height = 6))
  })
  for (table in c("profiles", "representatives")) local({
    key <- table
    output[[paste0("comparative_download_", key)]] <- downloadHandler(
      filename = function() comparative_download_name(paste0("_", key, ".tsv")),
      content = function(file) file.copy(file.path(comparative_run()$dir, paste0(key, ".tsv")), file))
  })

  group_genes <- eventReactive(input$comparative_group_analyze, {
    pasted <- unlist(strsplit(input$comparative_group_text %||% "", "\\s|,|;"))
    pasted <- pasted[nzchar(trimws(pasted))]
    uploaded <- character()
    file <- input$comparative_group_file
    if (!is.null(file$datapath) && file.exists(file$datapath)) {
      table <- tryCatch(read.csv(file$datapath, stringsAsFactors = FALSE,
                                 check.names = FALSE), error = function(e) data.frame())
      if (nrow(table) && ncol(table)) uploaded <- as.character(table[[1]])
    }
    unique(trimws(c(pasted, uploaded)))
  }, ignoreInit = TRUE)

  # The same command and the same representative rule as Comparative View, asked for many genes at once. It
  # reports calls rather than insertion profiles, which is what the heatmap needs and costs no hit-file reading,
  # so a group of a hundred genes answers in about the time one profiled gene takes.
  group_run <- reactive({
    genes <- group_genes(); req(length(genes) > 0)
    dir <- file.path(repo_root, "output", "gene_profiles", "_groups",
                     substr(digest::digest(sort(genes)), 1, 16))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    list_file <- file.path(dir, "genes.txt")
    writeLines(genes, list_file)
    args <- c(file.path(repo_root, "scripts", "local", "ytab_gene_profiles.py"),
              "--genes", list_file, "--output", dir,
              "--min-density", as.character(input$comparative_group_min_density %||% 18), "--json")
    result <- tryCatch(run_process_sync(python_bin(), args, wd = repo_root),
                       error = function(e) list(status = 1L, stdout = "", stderr = conditionMessage(e)))
    if (!identical(as.integer(result$status), 0L))
      return(list(status = "failure", error = paste(c(result$stderr, result$stdout), collapse = "\n")))
    manifest <- tryCatch(jsonlite::fromJSON(file.path(dir, "manifest.json")), error = function(e) list())
    profiles_file <- file.path(dir, "profiles.tsv")
    list(status = "success", dir = dir, genes = genes, manifest = manifest,
         flank = manifest$flank_fraction %||% NA_real_,
         calls = ytab_read_analysis_tsv(file.path(dir, "group_calls.tsv")),
         profiles = if (file.exists(profiles_file)) ytab_read_analysis_tsv(profiles_file) else data.frame())
  })

  output$comparative_group_state <- renderUI({
    if ((input$comparative_group_analyze %||% 0L) == 0L)
      return(ytab_empty_state("Paste or upload a gene list",
        "Each gene is looked for in every species at once and reported with the call its representative libraries give it."))
    run <- group_run()
    if (identical(run$status, "failure"))
      return(ytab_result_card("Could not analyse this gene group", status = "warning", tags$p(run$error %||% "")))
    NULL
  })

  # One tile per gene per species. A species can hold several counterparts of one gene, so the tile states what
  # they agree on: a family where some members are essential and some are not is "mixed", which is a result in
  # its own right and must not be flattened to either call.
  group_frame <- reactive({
    run <- group_run(); req(identical(run$status, "success"))
    calls <- run$calls
    validate(need(nrow(calls) > 0, "No gene in this list matched anything."))
    shown <- calls[as.character(calls$represented) %in% c("1", "TRUE", "True"), , drop = FALSE]
    validate(need(nrow(shown) > 0, "No species has libraries dense enough to call these genes."))
    parts <- lapply(split(shown, list(shown$query, shown$label), drop = TRUE), function(d) {
      call <- ytab_text_column(d, "call")
      gene <- ytab_text_column(d, "gene")
      route <- ytab_text_column(d, "matched_by")
      made <- call[nzchar(call)]
      state <- if (!length(made)) {
        if (any(nzchar(gene))) "no call" else "no counterpart"
      } else if (all(made == "ES")) "ES" else if (all(made != "ES")) "not ES" else "mixed"
      data.frame(query = d$query[[1]], label = d$label[[1]], state = state,
                 n = sum(nzchar(gene)),
                 routes = paste(sort(unique(route[nzchar(route)])), collapse = ", "),
                 genes = paste(gene[nzchar(gene)], collapse = ", "),
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, parts)
  })

  output$comparative_group_summary <- renderUI({
    run <- group_run(); req(identical(run$status, "success"))
    frame <- group_frame()
    resolved <- length(unique(frame$query[frame$n > 0]))
    tags$div(class = "ytab-stat-grid",
      tags$div(tags$b(length(run$genes)), "Genes in the list"),
      tags$div(tags$b(resolved), "Found in at least one species"),
      tags$div(tags$b(length(run$genes) - resolved), "Matched nothing"))
  })

  # The same figure Comparative View draws, once per gene in the list: one panel per counterpart, the species
  # and the gene naming it, the heights scaled within each panel. Built here rather than shared with
  # comparative_figure() because the strip has to carry the queried gene as well, and conflating the two
  # signatures made both harder to read than the handful of lines they have in common.
  group_profile_figure <- function(html = TRUE) {
    run <- group_run(); req(identical(run$status, "success"))
    calls <- run$calls; profiles <- run$profiles
    shown <- calls[as.character(calls$represented) %in% c("1", "TRUE", "True") &
                     nzchar(ytab_text_column(calls, "gene")), , drop = FALSE]
    validate(need(nrow(shown) > 0, "No gene in this list matched anything in a species with dense enough libraries."))
    style <- input$comparative_group_track_style %||% "reads"
    gene <- ytab_text_column(shown, "gene")
    # The species' own name for the gene beside its systematic id. Without it a family reads as the same gene
    # repeated: YGR088W next to YDR256C is only obviously CTT1 beside CTA1 once both are named.
    gene_name <- ytab_text_column(shown, "gene_name")
    shown_gene <- ifelse(nzchar(gene_name) & gene_name != gene, sprintf("%s (%s)", gene_name, gene), gene)
    call <- ytab_text_column(shown, "call")
    route <- ytab_text_column(shown, "matched_by")
    species <- if (html) paste0("<i>", shown$label, "</i>") else as.character(shown$label)
    mark <- ifelse(route == "orthogroup", " \u00b7 same gene family",
            ifelse(route == "name", " \u00b7 name match only", ""))
    panels <- data.frame(
      key = paste(shown$query, shown$species, gene, sep = "\u0001"),
      panel_label = sprintf("%s \u2014 %s\n%s \u00b7 %s%s", shown$query, species, shown_gene,
                            ifelse(nzchar(call), call, "no call"), mark),
      query = shown$query, label = as.character(shown$label), stringsAsFactors = FALSE)
    panels <- panels[order(panels$query, panels$label, panels$key), , drop = FALSE]
    panels <- panels[!duplicated(panels$key), , drop = FALSE]
    levels_ <- unique(panels$panel_label)

    drawn <- profiles
    if (nrow(drawn)) {
      drawn$key <- paste(drawn$query, drawn$species, drawn$gene, sep = "\u0001")
      drawn <- drawn[drawn$key %in% panels$key, , drop = FALSE]
      drawn$panel <- factor(panels$panel_label[match(drawn$key, panels$key)], levels = levels_)
    }
    empty <- panels[!panels$key %in% unique(drawn$key), , drop = FALSE]
    flank <- suppressWarnings(as.numeric(run$flank %||% NA)); if (is.na(flank)) flank <- .25
    x_range <- range(c(-flank, 1 + flank, if (nrow(drawn)) drawn$cds_fraction))
    ncol <- if (length(levels_) > 15L) 4L else if (length(levels_) > 8L) 3L else if (length(levels_) > 4L) 2L else 1L
    top <- if (identical(style, "reads"))
      stats::setNames(vapply(split(drawn$reads, drawn$panel, drop = FALSE),
                             function(x) if (length(x)) max(x) else 10, numeric(1)), levels_) else
      stats::setNames(rep(1.08, length(levels_)), levels_)
    band <- data.frame(panel = factor(levels_, levels = levels_),
                       ymin = if (identical(style, "reads")) 1 else 0,
                       ymax = as.numeric(top[levels_]), stringsAsFactors = FALSE)
    plot <- ggplot2::ggplot() +
      ggplot2::geom_rect(data = band, ggplot2::aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax),
                         fill = "grey92", inherit.aes = FALSE)
    if (nrow(drawn)) {
      drawn$hover <- sprintf("%s in %s\n%s\n%s reads at %.0f%% along the coding sequence",
                             drawn$query, drawn$label, drawn$gene, drawn$reads, 100 * drawn$cds_fraction)
      if (identical(style, "reads")) {
        plot <- plot +
          ggplot2::geom_segment(data = drawn, ggplot2::aes(x = cds_fraction, xend = cds_fraction, y = 1,
                                                           yend = reads, colour = label, text = hover),
                                linewidth = .3) +
          ggplot2::geom_point(data = drawn, ggplot2::aes(x = cds_fraction, y = reads, colour = label,
                                                         text = hover), size = .5)
      } else {
        drawn$width_bin <- cut(drawn$reads, breaks = c(-Inf, 10, 100, Inf), labels = c(".35", ".7", "1.2"))
        plot <- plot +
          ggplot2::geom_segment(data = drawn, ggplot2::aes(x = cds_fraction, xend = cds_fraction, y = 0, yend = 1,
                                                           colour = label, text = hover,
                                                           linewidth = as.numeric(as.character(width_bin))))
      }
    }
    if (identical(style, "sites")) {
      counts <- if (nrow(drawn))
        do.call(rbind, lapply(split(drawn, drawn$panel, drop = TRUE), function(d) data.frame(
          panel = d$panel[[1]], n = sum(tolower(as.character(d$in_cds)) %in% c("true", "1", "yes")),
          stringsAsFactors = FALSE))) else NULL
      if (!is.null(counts) && nrow(counts))
        plot <- plot + ggplot2::geom_text(data = counts, ggplot2::aes(x = x_range[2], y = 1.045,
                                          label = sprintf("n = %d in CDS", n)), hjust = 1, vjust = 0,
                                          size = 3.2, fontface = "bold", colour = "grey20", inherit.aes = FALSE)
      plot <- plot + ggplot2::scale_linewidth_identity() +
        ggplot2::scale_y_continuous(breaks = NULL, limits = c(0, 1.14)) +
        ggplot2::labs(y = "Insertion sites")
    } else {
      plot <- plot + ggplot2::scale_y_log10() + ggplot2::labs(y = "Reads per insertion site (log scale)")
    }
    if (nrow(empty))
      plot <- plot + ggplot2::geom_text(
        data = data.frame(panel = factor(empty$panel_label, levels = levels_),
                          note = "no insertion in the gene or its flanks", stringsAsFactors = FALSE),
        ggplot2::aes(x = mean(x_range), y = if (identical(style, "sites")) .5 else 3, label = note),
        size = 3.2, fontface = "italic", colour = "grey45", inherit.aes = FALSE)
    plot +
      ggplot2::facet_wrap(~panel, ncol = ncol, scales = "free_y", strip.position = "top", drop = FALSE) +
      ggplot2::scale_colour_manual(values = ytab_species_colour(unique(panels$label)),
                                   breaks = unique(panels$label), name = NULL) +
      ggplot2::scale_x_continuous(breaks = c(0, .5, 1), labels = c("start", "50%", "stop"), limits = x_range) +
      ggplot2::labs(x = sprintf("Position along the coding sequence (shaded; flanks are %g%% of its length)",
                                100 * flank)) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(legend.position = "none", strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(face = "bold", hjust = 0, size = 9, lineheight = 1.15),
                     axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
                     axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)))
  }

  output$comparative_group_profiles <- plotly::renderPlotly(
    ytab_ggplotly(group_profile_figure(TRUE), tooltip = "text"))

  output$comparative_group_profiles_caption <- renderUI({
    run <- group_run(); req(identical(run$status, "success"))
    reps <- run$manifest$representatives
    if (is.null(reps) || !length(reps)) return(NULL)
    # The manifest keys representatives by species key; the display label lives in the calls table.
    label_of <- stats::setNames(as.character(run$calls$label), as.character(run$calls$species))
    tags$p(class = "text-muted", HTML(sprintf(
      "%s Each panel is one gene in one species, scaled within itself. Drawn from: %s.",
      if (identical(input$comparative_group_track_style %||% "reads", "reads"))
        "The axis is reads at each insertion site, on a log scale." else
        "Every insertion site is one tick, with the reads at it in its width.",
      paste(vapply(names(reps), function(s) sprintf("%s from %s (%s libraries, %s sites/coding kb)",
        ytab_species_italic(label_of[[s]] %||% s), htmltools::htmlEscape(reps[[s]]$project),
        length(reps[[s]]$libraries), reps[[s]]$sites_per_coding_kb), character(1)), collapse = "; "))))
  })

  output$comparative_group_download_profiles <- downloadHandler(
    filename = function() "gene_group_profiles.tsv",
    content = function(file) file.copy(file.path(group_run()$dir, "profiles.tsv"), file))

  output$comparative_group_plot <- plotly::renderPlotly({
    frame <- group_frame()
    validate(need(nrow(frame) > 0, "No gene in this list has a result yet."))
    ytab_ggplotly(group_heatmap(frame, html = TRUE), tooltip = "text")
  })

  # One builder for both renderings: the browser gets the HTML spelling of the italics, a download the plotmath one.
  group_heatmap <- function(frame, html = TRUE) {
    # Genes with the most essential calls first, so the interesting block is at the top.
    score <- stats::aggregate(list(k = as.integer(frame$state == "ES")), list(query = frame$query), sum)
    frame$query <- factor(frame$query, score$query[order(score$k, decreasing = FALSE)])
    frame$state <- factor(frame$state, c("ES", "mixed", "not ES", "no call", "no counterpart"))
    species <- levels(factor(frame$label))
    ggplot2::ggplot(frame, ggplot2::aes(label, query, fill = state,
      text = sprintf("%s in %s\n%s\n%s%s", query, label, state,
                     ifelse(n > 0, sprintf("%d counterpart(s): %s", n, genes), "no counterpart"),
                     ifelse(nzchar(routes), sprintf("\nmatched by: %s", routes), "")))) +
      ggplot2::geom_tile(colour = "white", linewidth = .4) +
      ggplot2::scale_fill_manual(values = c(ES = "#b2182b", mixed = "#ef8a62", `not ES` = "#67a9cf",
                                            `no call` = "grey80", `no counterpart` = "grey92"),
                                 drop = FALSE, name = NULL) +
      ggplot2::scale_x_discrete(labels = ytab_species_italic(species, html)) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold"),
                     axis.text.y = ggplot2::element_text(size = 8))
  }

  output$comparative_group_plot_caption <- renderUI({
    frame <- group_frame(); if (!nrow(frame)) return(NULL)
    tagList(
      tags$p(class = "text-muted", sprintf(
        "%d genes across %d species. Red essential, amber mixed within a species' counterparts, blue not essential, grey no counterpart or no call.",
        length(unique(frame$query)), length(unique(frame$label)))),
      if (any(frame$routes != "" & frame$routes != "ortholog"))
        tags$p(class = "text-muted",
               "Some tiles rest on a shared gene family or on the gene's name alone rather than on a curated ortholog; the tooltip names the route for each."))
  })

  output$comparative_group_unmapped <- DT::renderDT({
    frame <- group_frame()
    missing <- setdiff(group_run()$genes, unique(frame$query[frame$n > 0]))
    if (!length(missing)) return(NULL)
    DT::datatable(data.frame(gene = missing, stringsAsFactors = FALSE), rownames = FALSE, selection = "none")
  })

  output$comparative_group_availability <- DT::renderDT({
    run <- group_run(); req(identical(run$status, "success"))
    calls <- run$calls
    cols <- intersect(c("label", "project", "sites_per_coding_kb", "represented", "status"), names(calls))
    data <- unique(calls[, cols, drop = FALSE])
    names(data)[names(data) == "label"] <- "Species"
    data$Species <- ytab_species_italic(data$Species)
    DT::datatable(data, rownames = FALSE, selection = "none", escape = FALSE,
                  options = list(dom = "t", scrollX = TRUE, pageLength = 25))
  })

  output$comparative_group_results <- DT::renderDT({
    run <- group_run(); req(identical(run$status, "success"))
    shown <- run$calls
    if ("label" %in% names(shown)) shown$label <- ytab_species_italic(shown$label)
    DT::datatable(shown, rownames = FALSE, selection = "none", escape = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  group_ggplot <- reactive({
    frame <- group_frame(); req(nrow(frame) > 0)
    group_heatmap(frame, html = FALSE)
  })

  output$comparative_group_download_png <- downloadHandler(
    filename = function() "gene_group_calls.png",
    content = function(file) ytab_save_ggplot(file, group_ggplot(), width = 9, height = 7))
  output$comparative_group_download_pdf <- downloadHandler(
    filename = function() "gene_group_calls.pdf",
    content = function(file) ytab_save_ggplot(file, group_ggplot(), width = 9, height = 7))
  output$comparative_group_download_calls <- downloadHandler(
    filename = function() "gene_group_calls.tsv",
    content = function(file) file.copy(file.path(group_run()$dir, "group_calls.tsv"), file))
}
