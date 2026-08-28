genome_browser_reference_dir <- function(project) {
  reference <- project$reference %||% list()
  fasta <- ytab_resolve_path(reference$fasta %||% "", project$repo_root %||% "")
  if (nzchar(fasta)) dirname(fasta) else ""
}

genome_browser_alias_rows <- function(project) {
  map <- qc_plot_glabrata_chromosome_map(project)
  if (!nrow(map)) return(data.frame())
  compact <- gsub(" ", "", map$display)
  data.frame(
    canonical = compact,
    alias_1 = map$display,
    alias_2 = map$contig,
    alias_3 = paste0(compact, "_C_glabrata_CBS138"),
    stringsAsFactors = FALSE
  )
}

genome_browser_alias_dir <- function(project) {
  path <- file.path(tempdir(), "ytab_igv_aliases", project$project_id %||% "project")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

genome_browser_alias_file <- function(project) {
  rows <- genome_browser_alias_rows(project)
  if (!nrow(rows)) return("")
  path <- file.path(genome_browser_alias_dir(project), "chromosome_aliases.tsv")
  write.table(rows, path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

genome_browser_asset_dir <- function(project) {
  path <- file.path(tempdir(), "ytab_igv_assets", project$project_id %||% "project")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

genome_browser_asset_url <- function(project, path) {
  genome_browser_file_url(path, genome_browser_asset_dir(project), "ytab-igv-assets")
}

genome_browser_contig_to_display <- function(project) {
  map <- qc_plot_glabrata_chromosome_map(project)
  if (!nrow(map)) return(character())
  stats::setNames(gsub(" ", "", map$display), map$contig)
}

genome_browser_write_indexed_fasta <- function(input_fasta, output_fasta, project) {
  contig_map <- genome_browser_contig_to_display(project)
  lines <- readLines(input_fasta, warn = FALSE)
  headers <- grep("^>", lines)
  if (!length(headers)) return(FALSE)
  output_fai <- paste0(output_fasta, ".fai")
  dir.create(dirname(output_fasta), recursive = TRUE, showWarnings = FALSE)
  con <- file(output_fasta, open = "wb")
  on.exit(close(con), add = TRUE)
  offset <- 0L
  fai <- list()
  write_bytes <- function(text) {
    bytes <- charToRaw(text)
    writeBin(bytes, con)
    length(bytes)
  }
  for (i in seq_along(headers)) {
    start <- headers[[i]]
    end <- if (i < length(headers)) headers[[i + 1L]] - 1L else length(lines)
    original <- sub("^>(\\S+).*$", "\\1", lines[[start]])
    chrom <- contig_map[[original]] %||% original
    offset <- offset + write_bytes(paste0(">", chrom, " ", original, "\n"))
    seq_lines <- lines[(start + 1L):end]
    seq <- paste(seq_lines[!grepl("^>", seq_lines)], collapse = "")
    seq <- gsub("\\s+", "", seq)
    seq_len <- nchar(seq, type = "bytes")
    seq_offset <- offset
    width <- 80L
    if (seq_len > 0L) {
      starts <- seq(1L, seq_len, by = width)
      for (s in starts) {
        chunk <- substr(seq, s, min(seq_len, s + width - 1L))
        offset <- offset + write_bytes(paste0(chunk, "\n"))
      }
    }
    fai[[length(fai) + 1L]] <- data.frame(
      chrom = chrom,
      length = seq_len,
      offset = seq_offset,
      line_bases = width,
      line_width = width + 1L,
      stringsAsFactors = FALSE
    )
  }
  fai <- do.call(rbind, fai)
  write.table(fai, output_fai, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  TRUE
}

genome_browser_display_reference <- function(project) {
  reference <- project$reference %||% list()
  source_fasta <- ytab_resolve_path(reference$fasta %||% "", project$repo_root %||% "")
  if (!file.exists(source_fasta)) return(list(fasta = "", fai = ""))
  ref_dir <- file.path(genome_browser_asset_dir(project), "reference")
  output_fasta <- file.path(ref_dir, "ytab_igv_reference.fna")
  output_fai <- paste0(output_fasta, ".fai")
  stale <- !file.exists(output_fasta) || !file.exists(output_fai) ||
    file.info(output_fasta)$mtime < file.info(source_fasta)$mtime
  if (isTRUE(stale)) genome_browser_write_indexed_fasta(source_fasta, output_fasta, project)
  list(fasta = output_fasta, fai = output_fai)
}

genome_browser_gff_attribute <- function(attributes, key) {
  pattern <- paste0("(^|;)", key, "=([^;]+)")
  hit <- regexec(pattern, attributes)
  value <- regmatches(attributes, hit)
  vapply(value, function(x) {
    if (length(x) < 3L) return("")
    URLdecode(x[[3]])
  }, character(1))
}

genome_browser_clean_display_label <- function(x) {
  x <- trimws(as.character(x %||% ""))
  x[is.na(x) | x %in% c("0", "NA", "N/A", "None", "none", "NULL", "null", ".", "-")] <- ""
  x
}

genome_browser_first_display_label <- function(...) {
  values <- list(...)
  max_len <- max(c(1L, vapply(values, length, integer(1))))
  out <- rep("", max_len)
  for (value in values) {
    value <- genome_browser_clean_display_label(value)
    if (length(value) == 1L && max_len > 1L) value <- rep(value, max_len)
    if (length(value) < max_len) value <- rep_len(value, max_len)
    fill <- !nzchar(out) & nzchar(value)
    out[fill] <- value[fill]
  }
  out
}

genome_browser_gene_bed_file <- function(project) {
  reference <- project$reference %||% list()
  gff <- ytab_resolve_path(reference$gff %||% "", project$repo_root %||% "")
  if (!file.exists(gff)) return("")
  out <- file.path(genome_browser_asset_dir(project), "reference", "genes.v5.bed")
  stale <- !file.exists(out) || file.info(out)$mtime < file.info(gff)$mtime
  if (!isTRUE(stale)) return(out)
  data <- tryCatch(read.delim(gff, header = FALSE, sep = "\t", quote = "",
                              comment.char = "#", stringsAsFactors = FALSE),
                   error = function(e) data.frame())
  if (!nrow(data) || ncol(data) < 9L) return("")
  data <- data[data[[3]] == "gene", , drop = FALSE]
  if (!nrow(data)) return("")
  contig_map <- genome_browser_contig_to_display(project)
  chrom <- contig_map[as.character(data[[1]])]
  chrom[is.na(chrom)] <- as.character(data[[1]][is.na(chrom)])
  name <- genome_browser_gff_attribute(as.character(data[[9]]), "Name")
  locus <- genome_browser_gff_attribute(as.character(data[[9]]), "locus_tag")
  id <- genome_browser_gff_attribute(as.character(data[[9]]), "ID")
  feature_id <- genome_browser_first_display_label(locus, id, name)
  label <- genome_browser_first_display_label(name, locus, id)
  if (exists("ytab_join_glabrata_display", mode = "function") && nzchar(project$repo_root %||% "")) {
    anno <- tryCatch(
      ytab_join_glabrata_display(
        data.frame(feature_id = feature_id, stringsAsFactors = FALSE),
        project$repo_root,
        "feature_id"
      ),
      error = function(e) data.frame()
    )
    if (is.data.frame(anno) && "scer_gene_name" %in% names(anno)) {
      mapped_name <- genome_browser_first_display_label(
        anno$scer_gene_name,
        if ("cagl_id" %in% names(anno)) anno$cagl_id else ""
      )
      label <- mapped_name
    } else if (is.data.frame(anno) && "gene_display_name" %in% names(anno)) {
      mapped_name <- genome_browser_first_display_label(
        if ("scer_gene_name" %in% names(anno)) anno$scer_gene_name else "",
        if ("cagl_id" %in% names(anno)) anno$cagl_id else ""
      )
      label <- mapped_name
    }
  }
  label <- genome_browser_clean_display_label(label)
  keep_label <- nzchar(label)
  if (!any(keep_label)) return("")
  data <- data[keep_label, , drop = FALSE]
  chrom <- chrom[keep_label]
  label <- label[keep_label]
  bed <- data.frame(
    chrom = unname(chrom),
    start = pmax(0L, as.integer(data[[4]]) - 1L),
    end = as.integer(data[[5]]),
    name = label,
    score = 0L,
    strand = as.character(data[[7]]),
    stringsAsFactors = FALSE
  )
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  write.table(bed, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  out
}

genome_browser_hit_count_column <- function(data) {
  hits <- names(data)[tolower(trimws(names(data))) %in% c("hit count", "hit_count", "count")]
  if (length(hits)) hits[[1]] else ""
}

genome_browser_hit_position_column <- function(data) {
  hits <- names(data)[tolower(trimws(names(data))) %in% c("hit position", "hit_position", "position", "pos")]
  if (length(hits)) hits[[1]] else ""
}

genome_browser_chromosome_column <- function(data) {
  hits <- names(data)[tolower(trimws(names(data))) %in% c("chromosome", "chrom", "chr", "seqid")]
  if (length(hits)) hits[[1]] else ""
}

genome_browser_display_insertion_bed <- function(project, sample, source_file) {
  if (!file.exists(source_file)) return("")
  out_dir <- file.path(genome_browser_asset_dir(project), "insertion_tracks")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_dir, paste0(make.names(sample %||% tools::file_path_sans_ext(basename(source_file))), ".insertions.v3.bed"))
  stale <- !file.exists(out) || file.info(out)$mtime < file.info(source_file)$mtime
  if (!isTRUE(stale)) return(out)
  data <- tryCatch(read.delim(source_file, header = TRUE, sep = "\t", quote = "",
                              check.names = FALSE, stringsAsFactors = FALSE),
                   error = function(e) data.frame())
  chrom_col <- genome_browser_chromosome_column(data)
  pos_col <- genome_browser_hit_position_column(data)
  count_col <- genome_browser_hit_count_column(data)
  if (!nrow(data) || !nzchar(chrom_col) || !nzchar(pos_col)) return("")
  pos <- suppressWarnings(as.integer(data[[pos_col]]))
  keep <- !is.na(pos) & pos > 0L
  if (!any(keep)) return("")
  data <- data[keep, , drop = FALSE]
  pos <- pos[keep]
  contig_map <- genome_browser_contig_to_display(project)
  chrom <- contig_map[as.character(data[[chrom_col]])]
  chrom[is.na(chrom)] <- as.character(data[[chrom_col]][is.na(chrom)])
  bed <- data.frame(
    chrom = unname(chrom),
    start = pmax(0L, pos - 1L),
    end = pos,
    stringsAsFactors = FALSE
  )
  bed <- bed[order(bed$chrom, bed$start, bed$end), , drop = FALSE]
  write.table(bed, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  out
}

genome_browser_display_bedgraph <- function(project, source_file) {
  if (!file.exists(source_file)) return("")
  out_dir <- file.path(genome_browser_asset_dir(project), "tracks")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_dir, basename(source_file))
  stale <- !file.exists(out) || file.info(out)$mtime < file.info(source_file)$mtime
  if (!isTRUE(stale)) return(out)
  data <- tryCatch(read.delim(source_file, header = FALSE, sep = "\t", quote = "",
                              stringsAsFactors = FALSE),
                   error = function(e) data.frame())
  if (!nrow(data) || ncol(data) < 4L) return("")
  contig_map <- genome_browser_contig_to_display(project)
  chrom <- contig_map[as.character(data[[1]])]
  chrom[is.na(chrom)] <- as.character(data[[1]][is.na(chrom)])
  data[[1]] <- unname(chrom)
  write.table(data[, 1:4, drop = FALSE], out, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  out
}

genome_browser_file_url <- function(path, root, route) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!nzchar(path) || !nzchar(root) || !startsWith(path, paste0(root, "/"))) return("")
  relative <- substring(path, nchar(root) + 2L)
  encoded <- diagnostic_encode_relative_path(relative)
  if (!nzchar(encoded)) "" else paste0(route, "/", encoded)
}

genome_browser_track_format <- function(path) {
  if (grepl("_hits\\.txt$", basename(path), ignore.case = TRUE)) return("insertion_bed")
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("bw", "bigwig")) return("bigwig")
  if (ext %in% c("bedgraph", "bdg")) return("bedgraph")
  if (ext %in% "wig") return("wig")
  if (ext %in% "bed") return("bed")
  ext
}

genome_browser_track_priority <- function(path) {
  format <- genome_browser_track_format(path)
  if (identical(format, "insertion_bed")) return(0L)
  if (identical(format, "bedgraph")) return(1L)
  if (identical(format, "bigwig")) return(2L)
  if (identical(format, "wig")) return(3L)
  if (identical(format, "bed")) return(4L)
  9L
}

genome_browser_compact_sample_label <- function(sample) {
  label <- gsub("-parent-", " ", sample, fixed = TRUE)
  label <- gsub("1_5mM-", "", label, fixed = TRUE)
  label <- gsub("-treated", "-treated", label, fixed = TRUE)
  label <- gsub("-", " ", label, fixed = TRUE)
  trimws(label)
}

genome_browser_track_inventory <- function(project) {
  root <- project$project_root %||% ""
  track_root <- file.path(root, "create_hit_file")
  if (!dir.exists(track_root)) return(data.frame())
  hit_files <- list.files(
    track_root,
    pattern = "_hits\\.txt$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  files <- hit_files
  if (!length(files)) files <- list.files(
    track_root,
    pattern = "\\.(bw|bigwig|bedgraph|wig|bed)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(hit_files))
    files <- files[grepl("\\.insertions\\.(bw|bigwig|bedgraph|wig)$|\\.filter_1\\.bed$", files, ignore.case = TRUE)]
  if (!length(files)) return(data.frame())
  sample <- basename(dirname(files))
  data <- data.frame(
    sample = sample,
    track_name = sample,
    display_label = vapply(sample, genome_browser_compact_sample_label, character(1)),
    source_file = normalizePath(files, winslash = "/", mustWork = FALSE),
    format = vapply(files, genome_browser_track_format, character(1)),
    priority = vapply(files, genome_browser_track_priority, integer(1)),
    stringsAsFactors = FALSE
  )
  data <- data[order(data$sample, data$priority, data$source_file), , drop = FALSE]
  data <- data[!duplicated(data$sample), , drop = FALSE]
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
        name = "",
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
  chrom <- as.character(index[[1]][[1]])
  display <- qc_plot_chromosome_display(chrom, project)
  display <- gsub(" ", "", display)
  length <- suppressWarnings(as.integer(index[[2]][[1]]))
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
    if (identical(format, "insertion_bed")) {
      source_file <- genome_browser_display_insertion_bed(project, row$sample[[1]], source_file)
      url <- genome_browser_asset_url(project, source_file)
    } else if (identical(format, "bedgraph")) {
      source_file <- genome_browser_display_bedgraph(project, source_file)
      url <- genome_browser_asset_url(project, source_file)
    } else {
      url <- genome_browser_file_url(source_file, root, "ytab-project-output")
    }
    role <- row$role[[1]] %||% ""
    color <- if (identical(role, "treated")) "rgb(190, 95, 65)" else if (identical(role, "parent")) "rgb(70, 125, 165)" else "rgb(90, 90, 90)"
    list(
      id = row$track_id[[1]] %||% make.names(row$sample[[1]]),
      name = "",
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
    tracks = genome_browser_track_configs(project, rows)
  )
}
