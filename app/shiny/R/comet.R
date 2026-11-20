# Comet plots: a gene's reads under treatment against its reads in the control, coloured by gene set.
#
# The overlay is assembled in the app so that choosing "Comet" draws a plot immediately: the two axes are columns
# the results table already carries, and the sets are a lookup over files that ship with the species. Nothing here
# computes a statistic -- a gene is matched to a named set by its identifier, its own name or its ortholog's name,
# which is presentation, the same kind of work as colouring points by their call. `ytab comet` remains for
# publication figures and writes the same assignment from the same files.
#
# Ortholog names matter more than they look. In a non-model species most genes carry no name of their own: in
# C. glabrata about 746 of 5,234 genes are named in the reference annotation while 4,812 reach an S. cerevisiae
# ortholog, so a set declared on ortholog names covers the genome and the same set declared on native names covers
# a seventh of it.

COMET_UNASSIGNED <- "other"

# Every GMT shipped for a species: a named list of files, labelled by file name.
comet_set_files <- function(repo_root, species_dir) {
  dir <- tryCatch(file.path(repo_root, species_dir, "genesets"), error = function(e) "")
  if (!length(dir) || !nzchar(dir) || !dir.exists(dir)) return(character())
  files <- list.files(dir, "\\.gmt$", full.names = TRUE)
  stats::setNames(files, tools::file_path_sans_ext(basename(files)))
}

# One GMT: set name -> member vector. Name, description, then members, tab separated; # and blank lines ignored.
comet_read_gmt <- function(path) {
  if (!length(path) || !file.exists(path)) return(list())
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !startsWith(trimws(lines), "#")]
  out <- list()
  for (line in lines) {
    parts <- trimws(strsplit(line, "\t", fixed = TRUE)[[1L]])
    if (length(parts) < 3L || !nzchar(parts[[1L]])) next
    members <- parts[-(1:2)]
    out[[parts[[1L]]]] <- toupper(members[nzchar(members)])
  }
  out
}

# Gene identifier -> the ortholog's gene name, from whichever orthology table of the species fits these ids.
# The table is chosen by matching, not by name, so a species is supported by dropping its table in the directory.
comet_orthologs <- function(repo_root, species_dir, ids) {
  dir <- tryCatch(file.path(repo_root, species_dir, "orthologs"), error = function(e) "")
  if (!length(dir) || !nzchar(dir) || !dir.exists(dir) || !length(ids)) return(NULL)
  ids <- unique(as.character(ids))
  for (path in list.files(dir, "\\.txt$|\\.tsv$", full.names = TRUE)) {
    d <- tryCatch(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                                    colClasses = "character", nrows = -1L), error = function(e) NULL)
    if (is.null(d) || ncol(d) < 2L || !nrow(d)) next
    # The id column is the one that actually matches this run's genes; a table for another species matches nothing.
    hits <- vapply(d, function(col) sum(col %in% ids), integer(1L))
    if (!length(hits) || max(hits) < max(20L, length(ids) %/% 20L)) next
    id_column <- names(hits)[[which.max(hits)]]
    name_columns <- setdiff(grep("_gene_name$", names(d), value = TRUE), id_column)
    if (!length(name_columns)) next
    # Prefer a name column that is not the species' own: the point is the reference species' name.
    own <- grep(sub("_[^_]*$", "", id_column), name_columns, fixed = TRUE, value = TRUE)
    column <- if (length(setdiff(name_columns, own))) setdiff(name_columns, own)[[1L]] else name_columns[[1L]]
    key <- d[[id_column]]; value <- trimws(d[[column]])
    # "NA" and "protein_coding" are how these tables spell "no name", not names.
    keep <- nzchar(key) & nzchar(value) & !toupper(value) %in% c("NA", "PROTEIN_CODING")
    if (!any(keep)) next
    map <- value[keep]; names(map) <- key[keep]
    return(map[!duplicated(names(map))])
  }
  NULL
}

# Free text (commas, semicolons, whitespace or new lines) to an upper-case vector of gene names or ids.
comet_parse_genes <- function(text) {
  if (is.null(text) || !length(text) || !nzchar(trimws(text))) return(character())
  parts <- trimws(strsplit(as.character(text), "[,;\n\r\t ]+")[[1L]])
  toupper(unique(parts[nzchar(parts)]))
}

# Each gene's set: the first of `chosen` whose members name its identifier, its own name or its ortholog's name.
# First match wins so a gene in several sets is drawn once, and the priority is the order the user chose.
comet_assign <- function(ids, names_own, ortholog, sets, chosen) {
  n <- length(ids)
  out <- rep(COMET_UNASSIGNED, n)
  if (!length(chosen)) return(out)
  keys <- toupper(as.character(ids))
  own <- toupper(ifelse(is.na(names_own), "", as.character(names_own)))
  orth <- toupper(ifelse(is.na(ortholog), "", as.character(ortholog)))
  for (set in rev(chosen)) {                     # reverse so earlier sets overwrite later ones
    members <- sets[[set]]
    if (!length(members)) next
    hit <- keys %in% members | (nzchar(own) & own %in% members) | (nzchar(orth) & orth %in% members)
    out[hit] <- set
  }
  out
}

# The comet itself, shared by every tab that has a control and a treated mean per gene. `d` needs ORF, Name,
# mean_control and mean_treatment; `overlay` is a frame of feature_id, gene_set and label_name as the tab's comet
# reactive built it. One implementation, so the plot means the same thing wherever it is drawn.
comet_ggplot <- function(d, overlay, s, base, hover = NULL, interactive = FALSE) {
  d <- d[is.finite(d$mean_control) & is.finite(d$mean_treatment), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$gene_set <- overlay$gene_set[match(as.character(d$ORF), overlay$feature_id)]
  d$gene_set[is.na(d$gene_set)] <- COMET_UNASSIGNED
  d$label_name <- overlay$label_name[match(as.character(d$ORF), overlay$feature_id)]
  # "other" is drawn first and dark so every overlay reads on top of it, whatever sets were chosen.
  named <- setdiff(sort(unique(d$gene_set)), COMET_UNASSIGNED)
  d$gene_set <- factor(d$gene_set, levels = c(COMET_UNASSIGNED, named))
  d <- d[order(d$gene_set), , drop = FALSE]
  d$x <- log10(d$mean_control + 1); d$y <- log10(d$mean_treatment + 1)
  label <- if (is.null(hover)) sprintf("%s %s", d$ORF, ifelse(is.na(d$Name), "", d$Name)) else hover(d)
  d$hover <- sprintf("%s<br>%s", label,
                     ifelse(d$gene_set == COMET_UNASSIGNED, "no set", as.character(d$gene_set)))
  # "my genes" is the user's own list and is always the accent colour, so it reads the same in every tab.
  palette <- stats::setNames(c("#1a1a1a", grDevices::hcl.colors(max(length(named), 1L), "Dark 3")),
                             c(COMET_UNASSIGNED, named))
  if ("my genes" %in% named) palette[["my genes"]] <- "#c8553d"
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, colour = gene_set, text = hover)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40") +
    ggplot2::geom_point(alpha = s$point_opacity, size = s$point_size) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::labs(x = "Control (log10 value + 1)", y = "Treated (log10 value + 1)", colour = "Gene set") + base
  rows <- which(!is.na(d$label_name) & nzchar(d$label_name))
  if (!length(rows)) return(p)
  # cmp_label_plot tags a row with its Name when that differs from its ORF, so the name chosen for the overlay
  # (the ortholog's, where there is one) is what gets drawn.
  d$Name[rows] <- d$label_name[rows]
  cmp_label_plot(p, d, rows, "x", "y", s, interactive)
}

# The overlay a tab passes to comet_ggplot: each gene's set, the user's own genes on top, and the labels.
# Shared so that Treated vs Control and Multi-condition answer "which set is this gene in" the same way.
#   genes      a frame with ORF, Name, mean_control, mean_treatment
#   sets       a GMT as comet_read_gmt returns it
#   chosen     the set names to overlay, in priority order
#   mine       free text naming the user's own genes
comet_overlay <- function(genes, sets, chosen, orthologs, mine_text, label_n = 10L,
                          direction = "both", label_min_reads = 10) {
  if (!nrow(genes)) return(data.frame(feature_id = character(), gene_set = character(),
                                      label_name = character(), stringsAsFactors = FALSE))
  ortholog <- if (is.null(orthologs)) rep("", nrow(genes)) else unname(orthologs[as.character(genes$ORF)])
  ortholog[is.na(ortholog)] <- ""
  assignment <- comet_assign(genes$ORF, genes$Name, ortholog, sets, chosen)
  mine <- comet_parse_genes(mine_text)
  if (length(mine)) {
    hit <- toupper(as.character(genes$ORF)) %in% mine |
           toupper(ifelse(is.na(genes$Name), "", as.character(genes$Name))) %in% mine |
           toupper(ortholog) %in% mine
    assignment[hit] <- "my genes"          # the user's list outranks the declared sets
  }

  # Which genes are worth naming. Distance from the diagonal is measured on the plotted scale; a gene with one
  # or two reads in both conditions has a large ratio for no reason but counting noise, so it must carry real
  # support in at least one condition to earn a label.
  off <- log2((genes$mean_treatment + 1) / (genes$mean_control + 1))
  eligible <- is.finite(off) &
    (genes$mean_control >= label_min_reads | genes$mean_treatment >= label_min_reads)
  # "Calls shown" restricts the hits being named, as it does on the volcano and the MA plot. Every member of a
  # set stays coloured whichever direction is chosen: the point of the overlay is to show where the whole
  # family sits, and hiding half of it would make a family that moved look like one that did not.
  if (identical(direction, "depleted")) eligible <- eligible & off < 0
  if (identical(direction, "enriched")) eligible <- eligible & off > 0

  n <- suppressWarnings(as.integer(label_n)); if (!length(n) || is.na(n)) n <- 0L
  label <- logical(nrow(genes))
  if (n > 0L) {
    # Top hits **per set**, not across all of them: with one pool a large family crowds out every small one, so
    # a 90-gene ribosome set would take every label and leave none for the 29-gene ergosterol set beside it.
    for (set in setdiff(unique(assignment), COMET_UNASSIGNED)) {
      pool <- which(assignment == set & eligible)
      if (!length(pool)) next
      label[pool[order(-abs(off[pool]))][seq_len(min(n, length(pool)))]] <- TRUE
    }
  }
  if (length(mine)) label[assignment == "my genes" & is.finite(off)] <- TRUE
  name <- ifelse(nzchar(ortholog), ortholog,
                 ifelse(is.na(genes$Name), as.character(genes$ORF), as.character(genes$Name)))
  data.frame(feature_id = as.character(genes$ORF), gene_set = assignment,
             label_name = ifelse(label, name, ""), stringsAsFactors = FALSE)
}
