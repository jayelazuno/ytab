ytab_glabrata_annotation_lookup_file <- function(repo_root) {
  file.path(repo_root, "resources", "comparative", "orthology",
            "20260610_Cgla_CAGL_to_Scer_annotation_lookup.csv")
}

ytab_glabrata_annotation_required_columns <- function() {
  c("cagl_id", "gwk60_id_clean", "qng_id", "cgla_common_name",
    "cg_to_sc_relationship", "pre_WGD_Ancestor", "scer_gene_id",
    "scer_gene_name", "SGD_essentiality", "SGD_description")
}

ytab_glabrata_key <- function(x) {
  tolower(trimws(as.character(x %||% "")))
}

# Gene annotation of the other species (resources/species/<species>/annotation/gene_annotation.tsv, built by
# scripts/local/ytab_build_species_resources.py) in the lookup's columns. Gene IDs are species-specific, so one lookup
# serves every project: a feature ID matches only its own species' rows. These species are named by their own common
# names (name_source "own"); C. glabrata, through the ortholog lookup, by S. cerevisiae names (name_source "scer").
ytab_species_annotation_lookup <- function(repo_root) {
  paths <- Sys.glob(file.path(repo_root, "resources", "species", "*", "annotation", "gene_annotation.tsv"))
  rows <- lapply(paths, function(path) {
    x <- tryCatch(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "", na.strings = character()), error = function(e) NULL)
    if (is.null(x) || !nrow(x) || !"gene_id" %in% names(x)) return(NULL)
    col <- function(name) if (name %in% names(x)) as.character(x[[name]]) else rep("", nrow(x))
    data.frame(cagl_id = col("systematic_id"), gwk60_id_clean = col("gene_id"), qng_id = "", cgla_common_name = col("gene_name"),
               cg_to_sc_relationship = col("scer_relationship"), pre_WGD_Ancestor = "", scer_gene_id = col("scer_gene_id"),
               scer_gene_name = col("scer_gene_name"), SGD_essentiality = col("scer_essentiality"), SGD_description = col("description"),
               name_source = "own", stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), rows)) %||% data.frame(stringsAsFactors = FALSE)
}

ytab_load_glabrata_annotation_lookup <- function(repo_root) {
  path <- ytab_glabrata_annotation_lookup_file(repo_root)
  empty <- data.frame(stringsAsFactors = FALSE)
  required <- ytab_glabrata_annotation_required_columns()
  data <- if (file.exists(path)) tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) empty) else empty
  if (nrow(data) && !all(required %in% names(data))) data <- empty
  if (nrow(data)) data$name_source <- "scer"
  other <- ytab_species_annotation_lookup(repo_root)
  if (nrow(other)) data <- if (nrow(data)) { other[setdiff(names(data), names(other))] <- ""; rbind(data, other[names(data)]) } else other
  if (!nrow(data)) return(empty)
  data$gwk60_key <- ytab_glabrata_key(data$gwk60_id_clean)
  data$cagl_key <- ytab_glabrata_key(data$cagl_id)
  data$qng_key <- ytab_glabrata_key(data$qng_id)
  data$scer_gene_id_key <- ytab_glabrata_key(data$scer_gene_id)
  data$scer_gene_name_key <- ytab_glabrata_key(data$scer_gene_name)
  data$cgla_common_name_key <- ytab_glabrata_key(data$cgla_common_name)
  data
}

ytab_first_nonempty <- function(...) {
  values <- list(...)
  if (!length(values)) return(character())
  n <- max(vapply(values, length, integer(1)))
  out <- rep("", n)
  for (value in values) {
    value <- as.character(value %||% "")
    if (length(value) == 1L && n > 1L) value <- rep(value, n)
    value[is.na(value)] <- ""
    take <- !nzchar(out) & nzchar(trimws(value))
    out[take] <- value[take]
  }
  out
}

ytab_join_glabrata_display <- function(data, repo_root, id_col,
                                       original_id_col = "original_feature_id") {
  if (!is.data.frame(data) || !nrow(data) || !nzchar(id_col %||% "") ||
      !id_col %in% names(data)) return(data)
  lookup <- ytab_load_glabrata_annotation_lookup(repo_root)
  original_id <- as.character(data[[id_col]])
  if (!nrow(lookup)) {
    data[[original_id_col]] <- if (original_id_col %in% names(data)) data[[original_id_col]] else original_id
    data$cagl_display_id <- original_id
    data$gene_display_name <- original_id
    data$cg_to_sc_relationship_display <- ""
    return(data)
  }
  lookup <- lookup[!duplicated(lookup$gwk60_key) & nzchar(lookup$gwk60_key), ,
                   drop = FALSE]
  match_ix <- match(ytab_glabrata_key(original_id), lookup$gwk60_key)
  joined <- data
  if (!original_id_col %in% names(joined)) joined[[original_id_col]] <- original_id
  raw_cols <- setdiff(names(lookup), grep("_key$", names(lookup), value = TRUE))
  for (col in raw_cols) {
    joined[[col]] <- ifelse(is.na(match_ix), "", as.character(lookup[[col]][match_ix]))
  }
  joined$cagl_display_id <- ytab_first_nonempty(joined$cagl_id, original_id)
  # Common name first, the gene ID when there is none: S. cerevisiae names for C. glabrata, the species' own otherwise.
  own <- joined$name_source == "own"
  joined$gene_display_name <- ifelse(own,
                                     ytab_first_nonempty(joined$cgla_common_name, joined$cagl_display_id, original_id),
                                     ytab_first_nonempty(joined$scer_gene_name, joined$cgla_common_name,
                                                         joined$cagl_display_id, original_id))
  joined$name_source <- NULL
  joined$cg_to_sc_relationship_display <- ytab_first_nonempty(joined$cg_to_sc_relationship, "")
  joined
}

ytab_glabrata_annotation_name_map <- function() {
  c(
    cagl_id = "Gene ID",
    cagl_display_id = "Gene ID",
    gwk60_id_clean = "Feature ID",
    original_feature_id = "Original feature ID",
    qng_id = "QNG ID",
    cgla_common_name = "Common name",
    cg_to_sc_relationship = "S. cerevisiae orthology",
    cg_to_sc_relationship_display = "S. cerevisiae orthology",
    pre_WGD_Ancestor = "Pre-WGD ancestor",
    scer_gene_id = "S. cerevisiae gene ID",
    scer_gene_name = "S. cerevisiae gene name",
    gene_display_name = "Gene name",
    SGD_essentiality = "SGD essentiality",
    SGD_description = "SGD description"
  )
}

ytab_standardize_glabrata_annotation_names <- function(data) {
  if (!is.data.frame(data) || !ncol(data)) return(data)
  data <- data[, setdiff(names(data), "cgla_gene_name_from_deseq"), drop = FALSE]
  map <- ytab_glabrata_annotation_name_map()
  hits <- intersect(names(data), names(map))
  names(data)[match(hits, names(data))] <- unname(map[hits])
  data
}

ytab_glabrata_gene_detail_columns <- function(data) {
  if (!is.data.frame(data) || !nrow(data)) {
    return(data.frame(
      .ytab_gene_detail_name = character(),
      .ytab_sgd_description = character(),
      .ytab_sgd_essentiality = character(),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  gene_name <- ytab_first_nonempty(
    if ("gene_display_name" %in% names(data)) data$gene_display_name else "",
    if ("Gene name" %in% names(data)) data[["Gene name"]] else "",
    if ("scer_gene_name" %in% names(data)) data$scer_gene_name else "",
    if ("cagl_display_id" %in% names(data)) data$cagl_display_id else "",
    if ("Gene ID" %in% names(data)) data[["Gene ID"]] else ""
  )
  description <- ytab_first_nonempty(
    if ("SGD_description" %in% names(data)) data$SGD_description else "",
    if ("SGD description" %in% names(data)) data[["SGD description"]] else ""
  )
  essentiality <- ytab_first_nonempty(
    if ("SGD_essentiality" %in% names(data)) data$SGD_essentiality else "",
    if ("SGD essentiality" %in% names(data)) data[["SGD essentiality"]] else ""
  )
  data.frame(
    .ytab_gene_detail_name = gene_name,
    .ytab_sgd_description = description,
    .ytab_sgd_essentiality = essentiality,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

ytab_glabrata_regular_display_columns <- function(data, repo_root, id_col,
                                                  include_details = FALSE) {
  data <- ytab_join_glabrata_display(data, repo_root, id_col)
  out <- data.frame(
    `Gene ID` = data$cagl_display_id,
    `Gene name` = data$gene_display_name,
    `S. cerevisiae orthology` = data$cg_to_sc_relationship_display,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (isTRUE(include_details)) out <- cbind(out, ytab_glabrata_gene_detail_columns(data))
  out
}
