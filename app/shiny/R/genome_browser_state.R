genome_browser_reference_dir <- function(project) {
  reference <- project$reference %||% list()
  fasta <- ytab_resolve_path(reference$fasta %||% "", project$repo_root %||% "")
  if (nzchar(fasta)) dirname(fasta) else ""
}

# Genome browser files are written by the pipeline (scripts/local/ytab_browser_tracks.py,
# ytab.domain_explorer.browser_tracks) into <project>/genome_browser/; the app only reads them.
genome_browser_asset_dir <- function(project) {
  path <- file.path(project$project_root %||% "", "genome_browser")
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

genome_browser_manifest <- function(project) {
  path <- file.path(genome_browser_asset_dir(project), "manifest.json")
  m <- if (file.exists(path)) tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL) else NULL
  if (is.null(m) || !identical(as.integer(m$contract %||% 0L), 2L)) NULL else m
}

# Runs the browser-tracks pipeline (it rebuilds only when an input is newer than its manifest). TRUE when the files
# changed, so the caller can refresh what it shows.
genome_browser_ensure <- function(project, python, repo_root, config_path) {
  before <- file.mtime(file.path(genome_browser_asset_dir(project), "manifest.json"))
  result <- tryCatch(run_process_sync(python, c(file.path(repo_root, "scripts/local/ytab_browser_tracks.py"), "--project-config", config_path), wd = repo_root),
                     error = function(e) list(status = 1L, stderr = conditionMessage(e)))
  after <- file.mtime(file.path(genome_browser_asset_dir(project), "manifest.json"))
  list(ok = identical(as.integer(result$status %||% 1L), 0L), message = result$stderr %||% "", changed = !identical(before, after))
}

genome_browser_alias_rows <- function(project) {
  m <- genome_browser_manifest(project); path <- file.path(genome_browser_asset_dir(project), m$aliases %||% "")
  if (is.null(m) || !file.exists(path)) return(data.frame())
  rows <- strsplit(readLines(path, warn = FALSE), "\t")
  data.frame(canonical = vapply(rows, `[[`, "", 1L), stringsAsFactors = FALSE)
}

genome_browser_asset_url <- function(project, path) {
  genome_browser_file_url(path, genome_browser_asset_dir(project), "ytab-igv-assets")
}

genome_browser_contig_to_display <- function(project) {
  map <- qc_plot_glabrata_chromosome_map(project)
  if (!nrow(map)) return(character())
  stats::setNames(gsub(" ", "", map$display), map$contig)
}

genome_browser_display_reference <- function(project) {
  m <- genome_browser_manifest(project); dir <- genome_browser_asset_dir(project)
  if (is.null(m)) return(list(fasta = "", fai = ""))
  list(fasta = file.path(dir, m$reference), fai = file.path(dir, m$index))
}

genome_browser_gene_bed_file <- function(project) {
  m <- genome_browser_manifest(project)
  if (is.null(m)) "" else file.path(genome_browser_asset_dir(project), m$genes)
}

genome_browser_file_url <- function(path, root, route) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!nzchar(path) || !nzchar(root) || !startsWith(path, paste0(root, "/"))) return("")
  relative <- substring(path, nchar(root) + 2L)
  encoded <- diagnostic_encode_relative_path(relative)
  if (!nzchar(encoded)) "" else paste0(route, "/", encoded)
}

genome_browser_compact_sample_label <- function(sample) {
  label <- gsub("-", " ", sample, fixed = TRUE)
  trimws(label)
}

# One insertion track per sample, as listed by the pipeline; labels, roles and pools come from the sample sheet.
genome_browser_track_inventory <- function(project) {
  m <- genome_browser_manifest(project)
  if (is.null(m) || !length(m$tracks)) return(data.frame())
  dir <- genome_browser_asset_dir(project)
  sample <- vapply(m$tracks, function(t) as.character(t$sample), "")
  data <- data.frame(
    sample = sample,
    track_name = ytab_sample_labels(sample, project),
    display_label = vapply(ytab_sample_labels(sample, project), genome_browser_compact_sample_label, character(1)),
    source_file = file.path(dir, vapply(m$tracks, function(t) as.character(t$file), "")),
    format = "insertion_bed",
    stringsAsFactors = FALSE
  )
  data <- gene_domain_order_tracks(data, project$samples %||% data.frame())
  data$track_id <- make.names(data$sample, unique = TRUE)
  rownames(data) <- NULL
  data
}

genome_browser_preset_choices <- function(data) {
  if (!all(c("role", "pool") %in% names(data))) data <- gene_domain_order_tracks(data)
  choices <- gene_domain_preset_choices(data)
  choices <- choices[!unname(choices) %in% "matched_pairs"]
  pools <- gene_domain_matched_pools(data)
  if (length(pools)) {
    pool_labels <- paste("Pool", pools, "pair")
    names(pools) <- pool_labels
    choices <- c(
      choices[!unname(choices) %in% "custom"],
      stats::setNames(paste0("pool", pools, "_pair"), pool_labels),
      choices[unname(choices) %in% "custom"]
    )
  }
  choices
}

genome_browser_preset_rows <- function(preset, data) {
  gene_domain_preset_track_rows(preset, data)
}

genome_browser_reference_config <- function(project) {
  display_ref <- genome_browser_display_reference(project)
  fasta <- display_ref$fasta
  fai <- display_ref$fai
  if (!file.exists(fasta) || !file.exists(fai)) return(NULL)
  fasta_url <- genome_browser_asset_url(project, fasta)
  if (!nzchar(fasta_url)) return(NULL)
  alias_rows <- genome_browser_alias_rows(project)
  tracks <- list()
  gene_bed <- genome_browser_gene_bed_file(project)
  if (file.exists(gene_bed)) {
    gene_url <- genome_browser_asset_url(project, gene_bed)
    if (nzchar(gene_url)) {
      tracks <- list(list(
        id = "gene_annotations",
        name = "Gene annotations",
        type = "annotation",
        format = "bed",
        url = gene_url,
        displayMode = "COLLAPSED",
        height = 55,
        order = 1000000,
        showTrackLabel = FALSE,
        searchable = TRUE,
        visibilityWindow = 10000
      ))
    }
  }
  list(
    id = paste0("ytab_", project$project_id %||% "project"),
    name = paste(project$display_name %||% project$project_id %||% "YTAB project", "reference"),
    fastaURL = fasta_url,
    indexed = FALSE,
    chromosomeOrder = if (nrow(alias_rows)) as.list(alias_rows$canonical) else NULL,
    tracks = tracks
  )
}

genome_browser_default_locus <- function(project) {
  reference <- project$reference %||% list()
  fasta <- ytab_resolve_path(reference$fasta %||% "", project$repo_root %||% "")
  fai <- paste0(fasta, ".fai")
  if (!file.exists(fai)) return("")
  index <- tryCatch(read.delim(fai, header = FALSE, stringsAsFactors = FALSE), error = function(e) data.frame())
  if (!nrow(index) || ncol(index) < 2L) return("")
  # the first chromosome in the reference's display order (chromosome_names.tsv), else the first FASTA record
  first <- 1L
  names_file <- file.path(dirname(fasta), "chromosome_names.tsv")
  if (file.exists(names_file)) {
    order <- tryCatch(read.delim(names_file, stringsAsFactors = FALSE), error = function(e) data.frame())
    if (all(c("contig", "order") %in% names(order))) {
      rank <- suppressWarnings(as.numeric(order$order[match(as.character(index[[1]]), order$contig)]))
      if (any(!is.na(rank))) first <- which.min(ifelse(is.na(rank), Inf, rank))
    }
  }
  chrom <- as.character(index[[1]][[first]])
  display <- qc_plot_chromosome_display(chrom, project)
  display <- gsub(" ", "", display)
  length <- suppressWarnings(as.integer(index[[2]][[first]]))
  end <- if (is.na(length)) 50000L else min(length, 50000L)
  paste0(display, ":1-", end)
}

genome_browser_track_configs <- function(project, rows) {
  if (!is.data.frame(rows) || !nrow(rows)) return(list())
  root <- project$project_root %||% ""
  configs <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, , drop = FALSE]
    format <- row$format[[1]]
    source_file <- row$source_file[[1]]
    if (format %in% c("insertion_bed", "bedgraph")) {
      url <- genome_browser_asset_url(project, source_file)
    } else {
      url <- genome_browser_file_url(source_file, root, "ytab-project-output")
    }
    role <- row$role[[1]] %||% ""
    color <- if (identical(role, "treated")) "rgb(190, 95, 65)" else if (identical(role, "parent")) "rgb(70, 125, 165)" else "rgb(90, 90, 90)"
    list(
      id = row$track_id[[1]] %||% make.names(row$sample[[1]]),
      name = as.character(i),
      type = if (format %in% c("bigwig", "bedgraph", "wig")) "wig" else "annotation",
      format = if (identical(format, "insertion_bed")) "bed" else format,
      url = url,
      color = color,
      height = if (identical(format, "insertion_bed")) 42 else 55,
      displayMode = if (identical(format, "insertion_bed")) "COLLAPSED" else NULL,
      showTrackLabel = FALSE,
      searchable = if (identical(format, "insertion_bed")) FALSE else NULL,
      visibilityWindow = if (identical(format, "insertion_bed")) 1000000 else NULL,
      autoscale = TRUE,
      showDataRange = FALSE,
      min = -10000,
      max = 10000
    )
  })
  Filter(function(x) nzchar(x$url %||% ""), configs)
}

genome_browser_session_config <- function(project, rows, locus = "") {
  genome <- genome_browser_reference_config(project)
  if (is.null(genome)) return(NULL)
  locus <- trimws(as.character(locus %||% ""))
  if (!nzchar(locus)) locus <- genome_browser_default_locus(project)
  list(
    project_id = project$project_id %||% "",
    genome = genome,
    locus = locus,
    showTrackLabels = FALSE,
    tracks = genome_browser_track_configs(project, rows)
  )
}
