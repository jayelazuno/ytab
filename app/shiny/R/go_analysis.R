go_analysis_ontology_choices <- function() {
  c("Biological Process (BP)" = "BP",
    "Molecular Function (MF)" = "MF",
    "Cellular Component (CC)" = "CC")
}

go_analysis_output_root <- function(project) {
  file.path(project$project_root, "go_analysis")
}

go_analysis_discover_results <- function(project) {
  root <- go_analysis_output_root(project)
  if (!dir.exists(root)) return(list())
  manifests <- list.files(root, pattern = "^manifest\\.json$", recursive = TRUE,
                          full.names = TRUE)
  lapply(manifests, function(path) {
    meta <- tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE),
                     error = function(e) list())
    out <- as.character(meta$output_dir %||% dirname(path))
    list(
      analysis_id = as.character(meta$analysis_id %||% basename(dirname(path))),
      status = as.character(meta$status %||% ""),
      source_label = as.character(meta$source_label %||% ""),
      completed_at = as.character(meta$completed_at %||% ""),
      output_dir = out,
      manifest = path,
      topgo_results = as.character(meta$topgo_results %||%
                                     file.path(out, "tables", "topgo_results.csv")),
      significant_results = as.character(meta$significant_results %||%
                                           file.path(out, "tables", "topgo_significant.csv")),
      semantic_heatmap = as.character(meta$semantic_heatmap %||%
                                        file.path(out, "tables", "semantic_reduced_heatmap.csv")),
      target_file = as.character(meta$target_file %||% file.path(out, "inputs", "target_genes.csv")),
      target_gene_count = as.integer(meta$target_gene_count %||% NA_integer_),
      mapped_target_rows = as.integer(meta$mapped_target_rows %||% NA_integer_),
      universe_gene_rows = as.integer(meta$universe_gene_rows %||% NA_integer_)
    )
  })
}

go_analysis_latest_result <- function(results) {
  if (!length(results)) return(NULL)
  times <- vapply(results, function(x) {
    info <- file.info(x$manifest)
    if (!is.na(info$mtime)) as.numeric(info$mtime) else 0
  }, numeric(1))
  results[[which.max(times)]]
}

go_analysis_read_csv <- function(path) {
  if (!nzchar(path %||% "") || !file.exists(path)) return(data.frame())
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) data.frame())
}

go_analysis_gene_set_table <- function(selected_hits = data.frame(),
                                       pasted = "",
                                       uploaded_path = NULL,
                                       source = "fitness") {
  if (identical(source, "fitness")) {
    x <- selected_hits
    if (!is.data.frame(x) || !nrow(x)) {
      return(data.frame(input_gene = character(), source = character(),
                        stringsAsFactors = FALSE))
    }
    id <- ytab_first_nonempty(
      if ("feature_id" %in% names(x)) x$feature_id else "",
      if ("Original feature ID" %in% names(x)) x[["Original feature ID"]] else "",
      if ("cagl_display_id" %in% names(x)) x$cagl_display_id else "",
      if ("CAGL ID" %in% names(x)) x[["CAGL ID"]] else "",
      if ("gene_display_name" %in% names(x)) x$gene_display_name else "",
      if ("Gene name" %in% names(x)) x[["Gene name"]] else ""
    )
    out <- data.frame(input_gene = id, source = "Current Fitness selected top hits",
                      stringsAsFactors = FALSE)
    if ("selected_rank" %in% names(x)) out$Rank <- x$selected_rank
    if ("hit_direction" %in% names(x)) out$`Hit direction` <- x$hit_direction
    if ("cagl_display_id" %in% names(x)) out$`CAGL ID` <- x$cagl_display_id
    if ("gene_display_name" %in% names(x)) out$`Gene name` <- x$gene_display_name
    out <- out[nzchar(trimws(out$input_gene)), , drop = FALSE]
    return(out[!duplicated(out$input_gene), , drop = FALSE])
  }
  if (identical(source, "upload") && !is.null(uploaded_path) &&
      nzchar(uploaded_path) && file.exists(uploaded_path)) {
    text <- paste(readLines(uploaded_path, warn = FALSE), collapse = "\n")
  } else {
    text <- pasted
  }
  values <- unlist(strsplit(as.character(text %||% ""), "[,;[:space:]]+"), use.names = FALSE)
  values <- unique(trimws(values[nzchar(trimws(values))]))
  data.frame(input_gene = values,
             source = if (identical(source, "upload")) "Uploaded gene list" else "Pasted gene list",
             stringsAsFactors = FALSE)
}

go_analysis_fitness_selected_hits <- function(result, mode = "combined", pair = "",
                                              direction = "both", n = 10L,
                                              min_support = 0L, repo_root = NULL) {
  if (is.null(result)) return(data.frame())
  mode <- mode %||% "combined"
  pair <- pair %||% ""
  direction <- direction %||% "both"
  n <- suppressWarnings(as.integer(n %||% 10L))
  if (is.na(n) || n < 0L) n <- 10L
  min_support <- if (identical(mode, "combined")) suppressWarnings(as.integer(min_support %||% 0L)) else 0L
  if (is.na(min_support) || min_support < 0L) min_support <- 0L
  x <- fitness_ma_top_hits_table(result, mode, pair, direction, n, "top", "", min_support)
  if (nrow(x) && "feature_id" %in% names(x)) x <- ytab_join_glabrata_display(x, repo_root, "feature_id")
  x
}

go_analysis_selected_hit_display_table <- function(data) {
  if (!is.data.frame(data) || !nrow(data)) return(data.frame())
  keep <- intersect(c("selected_rank", "cagl_display_id", "gene_display_name",
                      "cg_to_sc_relationship_display", "hit_direction", "pair",
                      "log2fc", "mean_abundance", "rank_z_strength",
                      "supporting_pool_ids"), names(data))
  out <- data[, keep, drop = FALSE]
  for (nm in intersect(c("log2fc", "mean_abundance", "rank_z_strength"), names(out))) {
    out[[nm]] <- round(as.numeric(out[[nm]]), 3)
  }
  names(out)[names(out) == "selected_rank"] <- "Rank"
  names(out)[names(out) == "cagl_display_id"] <- "CAGL ID"
  names(out)[names(out) == "gene_display_name"] <- "Gene name"
  names(out)[names(out) == "cg_to_sc_relationship_display"] <- "Cg-to-Sc relationship"
  names(out)[names(out) == "hit_direction"] <- "Hit direction"
  names(out)[names(out) == "pair"] <- "Pair"
  names(out)[names(out) == "log2fc"] <- "Directional log2FC"
  names(out)[names(out) == "mean_abundance"] <- "CPM/read support"
  names(out)[names(out) == "rank_z_strength"] <- "Local z-score support"
  names(out)[names(out) == "supporting_pool_ids"] <- "Supporting pool IDs"
  out
}

go_analysis_format_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "",
         ifelse(x < 0.001, formatC(x, format = "e", digits = 2),
                formatC(x, format = "f", digits = 4)))
}

go_analysis_format_number <- function(x, digits = 2L) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

go_analysis_ontology_label <- function(code, fallback = "") {
  code <- as.character(code %||% "")
  labels <- c(BP = "Biological Process", MF = "Molecular Function", CC = "Cellular Component")
  out <- unname(labels[code])
  out[is.na(out) | !nzchar(out)] <- as.character(fallback %||% "")[seq_along(out)][is.na(out) | !nzchar(out)]
  out
}

go_analysis_wrap_text <- function(x, width = 56L) {
  x <- as.character(x %||% "")
  if (requireNamespace("stringr", quietly = TRUE)) {
    x <- stringr::str_wrap(x, width = width)
  }
  x
}

go_analysis_go_display_table <- function(data, ontology = "BP", p_cutoff = 0.05,
                                         top_n = 15L) {
  if (!is.data.frame(data) || !nrow(data)) return(data.frame())
  required <- c("GO.ID", "Term", "ontology_code", "weight01_numeric")
  if (!all(required %in% names(data))) return(data.frame())
  ontology <- as.character(ontology %||% "BP")
  p_cutoff <- suppressWarnings(as.numeric(p_cutoff %||% 0.05))
  if (!is.finite(p_cutoff)) p_cutoff <- 0.05
  top_n <- suppressWarnings(as.integer(top_n %||% 15L))
  if (is.na(top_n) || top_n < 1L) top_n <- 15L
  x <- data[data$ontology_code == ontology, , drop = FALSE]
  if (!nrow(x)) return(data.frame())
  x$weight01_numeric <- suppressWarnings(as.numeric(x$weight01_numeric))
  x$classic_numeric <- suppressWarnings(as.numeric(x$classic_numeric))
  x <- x[is.finite(x$weight01_numeric), , drop = FALSE]
  x <- x[order(x$weight01_numeric), , drop = FALSE]
  if (nrow(x)) {
    x$go_plot_rank <- seq_len(nrow(x))
    x$go_displayed_in_plot <- x$weight01_numeric < p_cutoff & x$go_plot_rank <= top_n
  }
  score <- -log10(pmax(x$weight01_numeric, 1e-300))
  out <- data.frame(
    `Plot term` = ifelse(isTRUE(length(x$go_displayed_in_plot) > 0L) &
                           !is.na(x$go_displayed_in_plot) &
                           x$go_displayed_in_plot, "yes", ""),
    Ontology = go_analysis_ontology_label(x$ontology_code, x$ontology_label),
    `GO ID` = x$GO.ID,
    Term = go_analysis_wrap_text(x$Term, width = 64L),
    `Annotated genes` = x$Annotated,
    `Significant genes` = x$Significant,
    Expected = go_analysis_format_number(x$Expected, 2L),
    `Classic Fisher P` = go_analysis_format_p(x$classic_numeric),
    `weight01 Fisher P` = go_analysis_format_p(x$weight01_numeric),
    `-log10(weight01 P)` = go_analysis_format_number(score, 2L),
    `Target genes with GO` = x$target_genes_with_GO,
    `Universe genes with GO` = x$universe_genes_with_GO,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  gene_cols <- intersect(c("Genes", "genes", "target_genes", "target_gene_names",
                           "significant_genes", "Significant genes list"), names(x))
  if (length(gene_cols)) {
    out$`Genes contributing to term` <- go_analysis_wrap_text(x[[gene_cols[[1L]]]], width = 80L)
  }
  out
}

go_analysis_semantic_display_table <- function(data, ontology = "all") {
  if (!is.data.frame(data) || !nrow(data)) return(data.frame())
  required <- c("go", "parent", "score", "size", "term", "parentTerm",
                "ontology_code")
  if (!all(required %in% names(data))) return(data.frame())
  ontology <- as.character(ontology %||% "all")
  x <- data
  if (!identical(ontology, "all")) {
    x <- x[x$ontology_code == ontology, , drop = FALSE]
  }
  if (!nrow(x)) return(data.frame())
  x$score <- suppressWarnings(as.numeric(x$score))
  x$termUniqueness <- suppressWarnings(as.numeric(x$termUniqueness))
  x$termUniquenessWithinCluster <- suppressWarnings(as.numeric(x$termUniquenessWithinCluster))
  x$termDispensability <- suppressWarnings(as.numeric(x$termDispensability))
  x <- x[order(x$ontology_code, -x$score, x$parentTerm, x$term), , drop = FALSE]
  out <- data.frame(
    Ontology = go_analysis_ontology_label(x$ontology_code, x$ontology_label),
    `Representative GO ID` = x$parent,
    `Representative theme` = go_analysis_wrap_text(x$parentTerm, width = 62L),
    `Member GO ID` = x$go,
    `Member GO term` = go_analysis_wrap_text(x$term, width = 62L),
    `Cluster size` = x$size,
    `Best enrichment score` = go_analysis_format_number(x$score, 2L),
    `Term uniqueness` = go_analysis_format_number(x$termUniqueness, 2L),
    `Uniqueness in cluster` = go_analysis_format_number(x$termUniquenessWithinCluster, 2L),
    Dispensability = go_analysis_format_number(x$termDispensability, 2L),
    `Semantic threshold` = go_analysis_format_number(x$semantic_similarity_threshold, 2L),
    `Semantic method` = x$semantic_similarity_method,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out
}

go_analysis_display_datatable <- function(data, page_length = 10L) {
  if (!is.data.frame(data) || !nrow(data)) return(NULL)
  DT::datatable(
    data,
    rownames = FALSE,
    filter = "top",
    class = "compact stripe hover ytab-go-result-table",
    selection = "none",
    escape = TRUE,
    options = list(
      pageLength = page_length,
      scrollX = TRUE,
      autoWidth = TRUE,
      columnDefs = list(
        list(className = "dt-left ytab-go-term-cell", targets = "_all")
      )
    )
  )
}

go_analysis_plot_top_terms <- function(data, ontology = "BP", top_n = 15L,
                                       p_cutoff = 0.05, text_size = 13L,
                                       plot_options = NULL) {
  old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE)
  if (!is.data.frame(data) || !nrow(data) || !"ontology_code" %in% names(data)) {
    plot.new(); title("No GO enrichment result is available")
    return(invisible(NULL))
  }
  data <- data[data$ontology_code == ontology, , drop = FALSE]
  if (!nrow(data) || !"weight01_numeric" %in% names(data)) {
    plot.new(); title(paste("No", ontology, "GO terms are available"))
    return(invisible(NULL))
  }
  ontology_values <- if ("ontology_label" %in% names(data)) unique(as.character(data$ontology_label)) else ontology
  ontology_values <- ontology_values[nzchar(ontology_values)]
  ontology_label <- if (length(ontology_values)) ontology_values[[1L]] else ontology
  data$weight01_numeric <- suppressWarnings(as.numeric(data$weight01_numeric))
  data <- data[is.finite(data$weight01_numeric) & data$weight01_numeric > 0 &
                 data$weight01_numeric < p_cutoff, , drop = FALSE]
  data <- data[order(data$weight01_numeric), , drop = FALSE]
  if (nrow(data) > top_n) data <- data[seq_len(top_n), , drop = FALSE]
  if (!nrow(data)) {
    plot.new(); title(paste0(ontology_label, ": no enriched terms below p < ", p_cutoff))
    return(invisible(NULL))
  }
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    top_n <- max(1L, suppressWarnings(as.integer(top_n %||% 15L)))
    text_size <- suppressWarnings(as.numeric(text_size %||% 13L))
    if (!is.finite(text_size)) text_size <- 13
    style <- as.character(plot_options$style %||% "clean")
    show_grid <- identical(plot_options$grid %||% "hide", "show")
    label_mode <- as.character(plot_options$label_mode %||% "full")
    wrap_width <- if (identical(label_mode, "compact")) 42L else 62L
    terms <- as.character(data$Term)
    if (requireNamespace("stringr", quietly = TRUE)) {
      terms <- stringr::str_wrap(terms, width = wrap_width)
    }
    plot_df <- data.frame(
      Term_wrapped = terms,
      score = -log10(pmax(data$weight01_numeric, 1e-300)),
      stringsAsFactors = FALSE
    )
    plot_df$Term_wrapped <- factor(plot_df$Term_wrapped, levels = rev(plot_df$Term_wrapped))
    grid_color <- if (identical(style, "presentation")) "#D8E0EA" else "#E5EAF0"
    bar_fill <- if (identical(style, "compact")) "#4477AA" else "#2B79B9"
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = score, y = Term_wrapped)) +
      ggplot2::geom_col(width = 0.74, fill = bar_fill, color = "black", linewidth = 0.95) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
      ggplot2::labs(
        title = paste0("Top enriched GO terms — ", ontology_label),
        x = "-log10(topGO weight01 Fisher p-value)",
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = text_size) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = if (show_grid) ggplot2::element_line(color = grid_color, linewidth = 0.35) else ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = text_size + 3, color = "black", margin = ggplot2::margin(b = 14)),
        axis.title.x = ggplot2::element_text(face = "bold", color = "black", size = text_size, margin = ggplot2::margin(t = 12)),
        axis.text.x = ggplot2::element_text(face = "bold", color = "black", size = max(9, text_size - 2)),
        axis.text.y = ggplot2::element_text(face = "bold", color = "black", size = max(9, text_size - 2), margin = ggplot2::margin(r = 10)),
        plot.margin = ggplot2::margin(t = 14, r = 34, b = 18, l = 46)
      )
    print(p)
    return(invisible(p))
  }
  terms <- as.character(data$Term)
  if (requireNamespace("stringr", quietly = TRUE)) terms <- stringr::str_wrap(terms, width = 58)
  scores <- -log10(pmax(data$weight01_numeric, 1e-300))
  par(mar = c(5.5, 13, 4, 2))
  barplot(rev(scores), horiz = TRUE, names.arg = rev(terms), las = 1,
          col = "#2B79B9", border = NA, cex.names = text_size / 16,
          xlab = "-log10(topGO weight01 Fisher p-value)",
          main = paste0("Top enriched GO terms — ", ontology_label))
  grid(nx = NA, ny = NULL, col = "#E5EAF0")
}

go_analysis_plot_semantic_heatmap <- function(data, text_size = 13L, top_n = 35L,
                                              plot_options = NULL) {
  old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE)
  required <- c("representative_term", "ontology_code", "score")
  if (!is.data.frame(data) || !nrow(data) || !all(required %in% names(data))) {
    plot.new(); title("No semantic-reduction heatmap is available")
    return(invisible(NULL))
  }
  data$score <- suppressWarnings(as.numeric(data$score))
  data <- data[is.finite(data$score) & data$score > 0, , drop = FALSE]
  if (!nrow(data)) {
    plot.new(); title("No semantic-reduced GO terms are available")
    return(invisible(NULL))
  }
  data <- data[order(-data$score), , drop = FALSE]
  keep_terms <- unique(data$representative_term)
  top_n <- max(5L, suppressWarnings(as.integer(top_n %||% 35L)))
  if (length(keep_terms) > top_n) keep_terms <- keep_terms[seq_len(top_n)]
  data <- data[data$representative_term %in% keep_terms, , drop = FALSE]
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    text_size <- suppressWarnings(as.numeric(text_size %||% 13L))
    if (!is.finite(text_size)) text_size <- 13
    style <- as.character(plot_options$style %||% "clean")
    plot_max <- 6
    data$score_capped <- pmin(data$score, plot_max)
    term_scores <- stats::aggregate(score ~ representative_term, data, max)
    term_scores <- term_scores[order(term_scores$score, decreasing = TRUE), , drop = FALSE]
    term_levels <- rev(term_scores$representative_term)
    ontology_levels <- intersect(c("BP", "MF", "CC"), unique(as.character(data$ontology_code)))
    if (!length(ontology_levels)) ontology_levels <- unique(as.character(data$ontology_code))
    if (requireNamespace("stringr", quietly = TRUE)) {
      labels <- setNames(stringr::str_wrap(term_levels, width = 58), term_levels)
    } else {
      labels <- setNames(term_levels, term_levels)
    }
    data$representative_term <- factor(data$representative_term, levels = term_levels)
    data$ontology_code <- factor(data$ontology_code, levels = ontology_levels)
    panel_border <- if (identical(style, "presentation")) "black" else "grey45"
    p <- ggplot2::ggplot(data, ggplot2::aes(x = ontology_code, y = representative_term, fill = score_capped)) +
      ggplot2::geom_tile(color = "grey45", linewidth = 0.55) +
      ggplot2::scale_y_discrete(labels = labels, drop = FALSE) +
      ggplot2::scale_fill_gradientn(
        colours = c("grey90", "#FDB3B3", "#D7191C", "#7F0000"),
        values = c(0, 0.35, 0.72, 1),
        limits = c(0, plot_max),
        breaks = c(0, 3, 6),
        labels = c("0", "3", "\u22656"),
        name = "-log10(topGO\nweight01 P)"
      ) +
      ggplot2::labs(
        title = "Semantic-reduced GO themes",
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = text_size) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = panel_border, fill = NA, linewidth = 0.8),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = text_size + 3, color = "black", margin = ggplot2::margin(b = 14)),
        axis.text.x = ggplot2::element_text(face = "bold", color = "black", size = max(9, text_size - 1)),
        axis.text.y = ggplot2::element_text(face = "bold", color = "black", size = max(8, text_size - 3), margin = ggplot2::margin(r = 8)),
        legend.title = ggplot2::element_text(face = "bold", color = "black", size = max(9, text_size - 2)),
        legend.text = ggplot2::element_text(face = "bold", color = "black", size = max(8, text_size - 3)),
        legend.position = "right",
        plot.margin = ggplot2::margin(t = 14, r = 34, b = 18, l = 52)
      )
    print(p)
    return(invisible(p))
  }
  ont <- intersect(c("BP", "MF", "CC"), unique(data$ontology_code))
  mat <- matrix(0, nrow = length(keep_terms), ncol = length(ont),
                dimnames = list(keep_terms, ont))
  for (i in seq_len(nrow(data))) {
    mat[data$representative_term[[i]], data$ontology_code[[i]]] <-
      max(mat[data$representative_term[[i]], data$ontology_code[[i]]], data$score[[i]], na.rm = TRUE)
  }
  labels <- rownames(mat)
  if (requireNamespace("stringr", quietly = TRUE)) labels <- stringr::str_wrap(labels, width = 50)
  par(mar = c(5, 13, 4, 5))
  image(t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = hcl.colors(80, "Reds 3"),
        xlab = "", ylab = "", main = "Semantic-reduced GO themes")
  axis(1, at = seq(0, 1, length.out = ncol(mat)), labels = colnames(mat), las = 1,
       cex.axis = text_size / 15)
  axis(2, at = seq(0, 1, length.out = nrow(mat)), labels = rev(labels), las = 1,
       cex.axis = text_size / 17)
  box()
  mtext("-log10(topGO weight01 Fisher p-value)", side = 4, line = 2.5, cex = text_size / 16)
}

go_analysis_plot_height <- function(input, id) {
  height <- ytab_plot_customization_values(input, id)$height
  height <- sub("px$", "", as.character(height))
  value <- suppressWarnings(as.integer(height))
  if (is.na(value)) 700L else value
}

go_analysis_ui <- function() {
  navset_tab(id = "go_analysis_tabs",
    nav_panel("Gene set", value = "gene_set",
      ytab_two_column_layout(
        controls = ytab_control_panel("GO analysis setup",
          radioButtons("go_gene_source", "Gene source",
            choices = c("Current Fitness selected top hits" = "fitness",
                        "Paste gene list" = "paste",
                        "Upload gene list" = "upload"),
            selected = "fitness"),
          conditionalPanel("input.go_gene_source == 'fitness'",
            tags$div(class = "ytab-control-subsection",
              uiOutput("go_fitness_result_selector"),
              selectInput("go_fitness_mode", "Comparison view",
                          choices = c("Combined across pools" = "combined",
                                      "Individual treated-control pair" = "individual"),
                          selected = "combined"),
              uiOutput("go_fitness_support_control"),
              conditionalPanel("input.go_fitness_mode == 'individual'",
                               uiOutput("go_fitness_pair_control")),
              selectInput("go_fitness_hit_direction", "Hit direction",
                          choices = c("Depleted in treatment" = "depleted",
                                      "Enriched in treatment" = "enriched",
                                      "Both depleted and enriched" = "both"),
                          selected = "both"),
              numericInput("go_fitness_top_n", "Number of genes", value = 10,
                           min = 1, max = 500, step = 1),
              uiOutput("go_fitness_selection_status")
            )
          ),
          conditionalPanel("input.go_gene_source == 'paste'",
            textAreaInput("go_pasted_genes", "Pasted genes", rows = 5,
                          placeholder = "CAGL IDs, GWK60 IDs, S. cerevisiae ORFs, or gene names")
          ),
          conditionalPanel("input.go_gene_source == 'upload'",
            fileInput("go_gene_file", "Upload gene list CSV/TXT", accept = c(".csv", ".tsv", ".txt"))
          ),
          checkboxGroupInput("go_ontologies", "Ontologies",
                             choices = go_analysis_ontology_choices(),
                             selected = c("BP", "MF", "CC")),
          numericInput("go_node_size", "topGO node size", value = 10, min = 1, step = 1),
          numericInput("go_p_cutoff", "Significance cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
          checkboxInput("go_semantic_reduction", "Run semantic reduction with rrvgo", TRUE),
          numericInput("go_rrvgo_threshold", "rrvgo similarity threshold", value = 0.70, min = 0.1, max = 0.99, step = 0.05),
          textInput("go_analysis_id", "Analysis ID", value = "current_selected_hits"),
          actionButton("run_go_analysis", "Run GO analysis", class = "btn-primary"),
          tags$details(class = "ytab-technical-details",
                       tags$summary("Technical details"),
                       verbatimTextOutput("go_analysis_status"))
        ),
        main = tagList(
          ytab_section_header("Selected gene set",
                              "GO enrichment uses S. cerevisiae ORF IDs resolved through the existing C. glabrata annotation lookup."),
          uiOutput("go_gene_set_summary"),
          DT::DTOutput("go_gene_set_table"),
          tags$div(class = "ytab-actions",
                   downloadButton("download_go_targets_current", "Download selected GO genes CSV"))
        )
      )
    ),
    nav_panel("GO enrichment", value = "enrichment",
      ytab_two_column_layout(
        controls = ytab_control_panel("GO plot controls",
          selectInput("go_result_ontology", "Ontology", choices = go_analysis_ontology_choices(), selected = "BP"),
          numericInput("go_top_terms", "Top terms to display", value = 15, min = 1, max = 60, step = 1),
          ytab_plot_customization_controls("go_enrichment", include_bars = FALSE,
                                           include_labels = TRUE, default_height = "large",
                                           default_bar_orientation = "horizontal")
        ),
        main = tagList(
          ytab_section_header("GO enrichment"),
          ytab_plot_frame(plotOutput("go_enrichment_plot", height = "700px"),
                          width = "wide", type = "app-rendered"),
          DT::DTOutput("go_enrichment_table"),
          tags$div(class = "ytab-actions",
                   downloadButton("download_go_topgo_results", "Download topGO results CSV"),
                   downloadButton("download_go_significant_results", "Download significant GO terms CSV"))
        )
      )
    ),
    nav_panel("Semantic reduction", value = "semantic_reduction",
      ytab_two_column_layout(
        controls = ytab_control_panel("Semantic reduction",
          selectInput("go_semantic_ontology", "Ontology",
                      choices = c("All ontologies" = "all", go_analysis_ontology_choices()),
                      selected = "all"),
          numericInput("go_semantic_top_terms", "Representative themes to display", value = 35, min = 5, max = 100, step = 5),
          ytab_plot_customization_controls("go_semantic", include_heatmap = FALSE,
                                           include_labels = FALSE, include_grid = FALSE,
                                           default_height = "large")
        ),
        main = tagList(
          ytab_section_header("Semantic-reduced GO themes",
                              "Representative GO themes after semantic-similarity reduction."),
          ytab_plot_frame(plotOutput("go_semantic_heatmap", height = "700px"),
                          width = "wide", type = "app-rendered"),
          DT::DTOutput("go_semantic_table"),
          tags$div(class = "ytab-actions",
                   downloadButton("download_go_semantic_results", "Download semantic reduction CSV"))
        )
      )
    )
  )
}

go_analysis_server <- function(input, output, session, repo_root, active,
                               active_project_path, jobs, python_bin,
                               selected_hits_reactive,
                               fitness_result_reactive = NULL) {
  go_fitness_result <- reactive({
    if (is.function(fitness_result_reactive)) fitness_result_reactive() else NULL
  })
  go_fitness_selected_hits <- reactive({
    result <- go_fitness_result()
    if (is.null(result)) return(data.frame())
    go_analysis_fitness_selected_hits(
      result = result,
      mode = input$go_fitness_mode %||% "combined",
      pair = input$go_fitness_pair %||% "",
      direction = input$go_fitness_hit_direction %||% "both",
      n = input$go_fitness_top_n %||% 10L,
      min_support = input$go_fitness_min_support %||% 0L,
      repo_root = repo_root
    )
  })
  go_gene_set <- reactive({
    upload <- input$go_gene_file
    uploaded_path <- if (is.null(upload)) NULL else upload$datapath
    selected_hits <- if (identical(input$go_gene_source %||% "fitness", "fitness")) {
      go_fitness_selected_hits()
    } else {
      selected_hits_reactive()
    }
    go_analysis_gene_set_table(
      selected_hits = selected_hits,
      pasted = input$go_pasted_genes %||% "",
      uploaded_path = uploaded_path,
      source = input$go_gene_source %||% "fitness"
    )
  })

  output$go_fitness_result_selector <- renderUI({
    result <- go_fitness_result()
    if (is.null(result)) return(tags$p(class = "text-muted", "No Fitness result is available for selected-hit gene sets."))
    tags$p(class = "text-muted", sprintf("Fitness result: %s", basename(result$output_dir %||% "")))
  })

  output$go_fitness_pair_control <- renderUI({
    pairs <- fitness_ma_pair_choices(go_fitness_result())
    selectInput("go_fitness_pair", "Pair", choices = if (length(pairs)) pairs else "No pair data available")
  })

  output$go_fitness_support_control <- renderUI({
    if (!identical(input$go_fitness_mode %||% "combined", "combined")) return(NULL)
    result <- go_fitness_result()
    pools <- length(fitness_ma_available_pools(result))
    if (!pools) return(tags$p(class = "text-muted", "Concordance filter is unavailable because by-pool support data are missing."))
    choices <- c("All support classes" = "0",
                 setNames(as.character(seq_len(pools)),
                          paste0(">= ", seq_len(pools), " supporting pool",
                                 ifelse(seq_len(pools) == 1L, "", "s"))))
    selected <- as.character(input$go_fitness_min_support %||% "0")
    if (!selected %in% unname(choices)) selected <- "0"
    selectInput("go_fitness_min_support", "Concordance filter",
                choices = choices, selected = selected)
  })

  output$go_fitness_selection_status <- renderUI({
    if (!identical(input$go_gene_source %||% "fitness", "fitness")) return(NULL)
    x <- go_fitness_selected_hits()
    if (!nrow(x)) return(tags$p(class = "text-muted", "No Fitness selected hits pass the current GO gene-set controls."))
    tags$p(class = "text-muted", sprintf("%d Fitness selected hits will be used as GO target genes.", nrow(x)))
  })

  output$go_gene_set_summary <- renderUI({
    x <- go_gene_set()
    if (!nrow(x)) return(ytab_empty_state("No genes selected", "Choose a gene source or enter a gene list."))
    tags$p(sprintf("%d genes are selected for GO analysis.", nrow(x)))
  })

  output$go_gene_set_table <- DT::renderDT({
    if (identical(input$go_gene_source %||% "fitness", "fitness")) {
      x <- go_analysis_selected_hit_display_table(go_fitness_selected_hits())
    } else {
      x <- go_gene_set()
    }
    if (!nrow(x)) return(NULL)
    ytab_datatable(x, page_length = 10L, filter = "top", scroll_x = TRUE)
  })

  go_results_tick <- reactiveVal(0L)

  go_results <- reactive({
    go_results_tick()
    active()
    go_analysis_discover_results(active())
  })

  latest_go_result <- reactive({
    go_analysis_latest_result(go_results())
  })

  observe({
    if (!identical(jobs$stage(), "go_analysis")) return(NULL)
    if (!jobs$job_is_running()) return(NULL)
    invalidateLater(1500, session)
    jobs$poll_job()
    if (!jobs$job_is_running()) {
      go_results_tick(go_results_tick() + 1L)
    }
    invisible(NULL)
  })

  output$go_analysis_status <- renderText({
    job <- jobs$current_job()
    latest <- latest_go_result()
    paste(c(
      paste("Current job stage:", jobs$stage()),
      paste("Command:", jobs$command()),
      paste("Status:", as.character(job$status %||% "idle")),
      if (!is.null(latest)) paste("Latest GO result:", latest$analysis_id) else "Latest GO result: none",
      jobs$read_job_stdout(),
      jobs$read_job_stderr()
    ), collapse = "\n")
  })

  observeEvent(input$run_go_analysis, {
    if (jobs$job_is_running()) {
      showNotification("A job is already running.", type = "warning")
      return()
    }
    genes <- go_gene_set()
    if (!nrow(genes)) {
      showNotification("No GO target genes are selected.", type = "error")
      return()
    }
    project <- active()
    input_dir <- file.path(project$project_root, "go_analysis", "inputs")
    dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
    target_file <- file.path(input_dir, "current_go_target_genes.csv")
    write.csv(genes, target_file, row.names = FALSE)
    ontologies <- paste(input$go_ontologies %||% c("BP", "MF", "CC"), collapse = ",")
    analysis_id <- input$go_analysis_id %||% "current_selected_hits"
    args <- c(
      file.path(repo_root, "scripts/local/ytab_run_go_analysis.py"),
      "--project-config", active_project_path(),
      "--analysis-id", analysis_id,
      "--gene-file", target_file,
      "--ontologies", ontologies,
      "--node-size", as.character(as.integer(input$go_node_size %||% 10)),
      "--p-cutoff", as.character(as.numeric(input$go_p_cutoff %||% 0.05)),
      "--top-terms", as.character(as.integer(input$go_top_terms %||% 15)),
      "--source-label", input$go_gene_source %||% "fitness",
      "--force"
    )
    if (isTRUE(input$go_semantic_reduction)) {
      args <- c(args, "--semantic-reduction",
                "--rrvgo-threshold", as.character(as.numeric(input$go_rrvgo_threshold %||% 0.70)))
    }
    jobs$start_job(python_bin(), args, active_project_path(), "go_analysis",
                   wd = repo_root, selected_items = genes$input_gene)
    showNotification("GO analysis started.", type = "message")
  }, ignoreInit = TRUE)

  topgo_table <- reactive({
    result <- latest_go_result()
    if (is.null(result)) return(data.frame())
    go_analysis_read_csv(result$topgo_results)
  })

  significant_table <- reactive({
    result <- latest_go_result()
    if (is.null(result)) return(data.frame())
    go_analysis_read_csv(result$significant_results)
  })

  semantic_table <- reactive({
    result <- latest_go_result()
    if (is.null(result)) return(data.frame())
    go_analysis_read_csv(result$semantic_heatmap)
  })

  output$go_enrichment_plot <- renderPlot({
    values <- ytab_plot_customization_values(input, "go_enrichment")
    go_analysis_plot_top_terms(topgo_table(), input$go_result_ontology %||% "BP",
                               input$go_top_terms %||% 15,
                               input$go_p_cutoff %||% 0.05,
                               values$text_size,
                               values)
  }, height = function() go_analysis_plot_height(input, "go_enrichment"))

  output$go_enrichment_table <- DT::renderDT({
    x <- go_analysis_go_display_table(
      topgo_table(),
      ontology = input$go_result_ontology %||% "BP",
      p_cutoff = input$go_p_cutoff %||% 0.05,
      top_n = input$go_top_terms %||% 15L
    )
    go_analysis_display_datatable(x, page_length = 10L)
  })

  output$go_semantic_heatmap <- renderPlot({
    values <- ytab_plot_customization_values(input, "go_semantic")
    x <- semantic_table()
    ontology <- input$go_semantic_ontology %||% "all"
    if (!identical(ontology, "all") && is.data.frame(x) && "ontology_code" %in% names(x)) {
      x <- x[x$ontology_code == ontology, , drop = FALSE]
    }
    go_analysis_plot_semantic_heatmap(x, values$text_size,
                                      input$go_semantic_top_terms %||% 35L,
                                      values)
  }, height = function() go_analysis_plot_height(input, "go_semantic"))

  output$go_semantic_table <- DT::renderDT({
    result <- latest_go_result()
    if (is.null(result)) return(NULL)
    path <- file.path(result$output_dir, "tables", "rrvgo_reduced_terms.csv")
    x <- go_analysis_read_csv(path)
    x <- go_analysis_semantic_display_table(
      x,
      ontology = input$go_semantic_ontology %||% "all"
    )
    go_analysis_display_datatable(x, page_length = 10L)
  })

  output$go_result_summary <- renderUI({
    result <- latest_go_result()
    if (is.null(result)) return(ytab_empty_state("No GO result yet", "Run GO analysis from the Gene set tab."))
    tags$dl(class = "ytab-meta",
      tags$dt("Analysis ID"), tags$dd(result$analysis_id),
      tags$dt("Target genes"), tags$dd(result$target_gene_count),
      tags$dt("Mapped targets"), tags$dd(result$mapped_target_rows),
      tags$dt("Universe rows"), tags$dd(result$universe_gene_rows),
      tags$dt("Output directory"), tags$dd(result$output_dir)
    )
  })

  output$download_go_topgo_results <- downloadHandler(
    filename = function() paste0(active()$project_id, ".go.topgo_results.csv"),
    content = function(file) {
      x <- topgo_table()
      if (!nrow(x)) x <- data.frame(message = "No topGO results available.")
      write.csv(x, file, row.names = FALSE)
    }
  )
  output$download_go_significant_results <- downloadHandler(
    filename = function() paste0(active()$project_id, ".go.significant_terms.csv"),
    content = function(file) {
      x <- significant_table()
      if (!nrow(x)) x <- data.frame(message = "No significant GO terms available.")
      write.csv(x, file, row.names = FALSE)
    }
  )
  output$download_go_semantic_results <- downloadHandler(
    filename = function() paste0(active()$project_id, ".go.semantic_reduction.csv"),
    content = function(file) {
      result <- latest_go_result()
      x <- if (is.null(result)) data.frame() else
        go_analysis_read_csv(file.path(result$output_dir, "tables", "rrvgo_reduced_terms.csv"))
      if (!nrow(x)) x <- data.frame(message = "No semantic reduction results available.")
      write.csv(x, file, row.names = FALSE)
    }
  )
  output$download_go_targets_current <- downloadHandler(
    filename = function() paste0(active()$project_id, ".go.selected_genes.csv"),
    content = function(file) {
      x <- if (identical(input$go_gene_source %||% "fitness", "fitness")) {
        go_analysis_selected_hit_display_table(go_fitness_selected_hits())
      } else {
        go_gene_set()
      }
      if (!nrow(x)) x <- data.frame(message = "No GO target genes selected.")
      write.csv(x, file, row.names = FALSE)
    }
  )
}
