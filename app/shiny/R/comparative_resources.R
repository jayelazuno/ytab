comparative_empty_gene_lookup <- function() {
  data.frame(
    gene_id = character(), species = character(), expression_id = character(),
    id_type = character(), source_info = character(), gene_name = character(),
    hog_id = character(), og_id = character(), stringsAsFactors = FALSE
  )
}

comparative_empty_orthogroups <- function() {
  data.frame(
    orthogroup_id = character(), species = character(), gene_id = character(),
    gene_name = character(), source_column = character(), stringsAsFactors = FALSE
  )
}

comparative_resource_path <- function(repo_root, ...) {
  file.path(repo_root, "resources", "comparative", ...)
}

comparative_species_defaults <- function() {
  data.frame(
    species = c("glabrata", "albicans", "cerevisiae", "lactis", "pombe"),
    label = c("Candida glabrata", "Candida albicans",
              "Saccharomyces cerevisiae", "Kluyveromyces lactis",
              "Schizosaccharomyces pombe"),
    enabled = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    ytab_species_key = c("glabrata", "albicans", "cerevisiae", "lactis", ""),
    placeholder = c(FALSE, FALSE, FALSE, FALSE, TRUE),
    note = c("", "", "", "",
             "Placeholder only; no orthology mapping is currently migrated."),
    aliases = I(list(
      c("cg", "cgla", "C_glabrata", "Candida glabrata"),
      c("ca", "calb", "C_albicans", "Candida albicans"),
      c("sc", "scer", "S_cerevisiae", "Saccharomyces cerevisiae"),
      c("kl", "klac", "K_lactis", "Kluyveromyces lactis"),
      character()
    )),
    stringsAsFactors = FALSE
  )
}

load_comparative_species_manifest <- function(repo_root) {
  path <- comparative_resource_path(
    repo_root, "manifest", "comparative_species_manifest.yaml"
  )
  if (!file.exists(path) || !requireNamespace("yaml", quietly = TRUE))
    return(comparative_species_defaults())
  data <- tryCatch(yaml::read_yaml(path), error = function(e) NULL)
  if (!is.list(data) || !length(data)) return(comparative_species_defaults())
  rows <- lapply(names(data), function(key) {
    item <- data[[key]] %||% list()
    data.frame(
      species = key,
      label = as.character(item$label %||% key),
      enabled = isTRUE(item$enabled),
      ytab_species_key = as.character(item$ytab_species_key %||% ""),
      placeholder = isTRUE(item$placeholder),
      note = as.character(item$note %||% ""),
      aliases = I(list(as.character(item$rnacross_aliases %||% character()))),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

load_comparative_gene_lookup <- function(repo_root) {
  path <- comparative_resource_path(
    repo_root, "orthology", "gene_lookup.csv"
  )
  if (!file.exists(path)) return(comparative_empty_gene_lookup())
  data <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                   error = function(e) comparative_empty_gene_lookup())
  required <- names(comparative_empty_gene_lookup())
  for (name in setdiff(required, names(data))) data[[name]] <- ""
  data[, required, drop = FALSE]
}

load_comparative_orthogroups <- function(repo_root) {
  path <- comparative_resource_path(
    repo_root, "orthology", "orthogroups_long.csv"
  )
  if (!file.exists(path)) return(comparative_empty_orthogroups())
  data <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                   error = function(e) comparative_empty_orthogroups())
  required <- names(comparative_empty_orthogroups())
  for (name in setdiff(required, names(data))) data[[name]] <- ""
  data[, required, drop = FALSE]
}

comparative_query_matches <- function(data, query, columns) {
  if (!nrow(data)) return(rep(FALSE, 0L))
  query <- trimws(tolower(as.character(query %||% "")))
  if (!nzchar(query)) return(rep(FALSE, nrow(data)))
  columns <- intersect(columns, names(data))
  if (!length(columns)) return(rep(FALSE, nrow(data)))
  Reduce(`|`, lapply(data[, columns, drop = FALSE], function(value) {
    tolower(trimws(as.character(value))) == query
  }))
}

find_orthogroup_for_gene <- function(gene_query, orthology, gene_lookup) {
  if (!is.data.frame(orthology)) orthology <- comparative_empty_orthogroups()
  if (!is.data.frame(gene_lookup)) gene_lookup <- comparative_empty_gene_lookup()
  query <- trimws(as.character(gene_query %||% ""))
  if (!nzchar(query)) return(data.frame())
  lookup_hits <- gene_lookup[comparative_query_matches(
    gene_lookup, query, c("gene_id", "gene_name", "expression_id", "og_id", "hog_id")
  ), , drop = FALSE]
  og_ids <- unique(as.character(lookup_hits$og_id %||% character()))
  og_ids <- og_ids[nzchar(og_ids)]
  orth_hits <- orthology[comparative_query_matches(
    orthology, query, c("gene_id", "gene_name", "orthogroup_id")
  ), , drop = FALSE]
  og_ids <- unique(c(og_ids, as.character(orth_hits$orthogroup_id %||% character())))
  og_ids <- og_ids[nzchar(og_ids)]
  if (!length(og_ids)) return(data.frame())
  data.frame(orthogroup_id = og_ids, stringsAsFactors = FALSE)
}

map_gene_across_species <- function(gene_query, orthology, gene_lookup) {
  groups <- find_orthogroup_for_gene(gene_query, orthology, gene_lookup)
  if (!nrow(groups) || !nrow(orthology)) return(comparative_empty_orthogroups())
  result <- orthology[orthology$orthogroup_id %in% groups$orthogroup_id, ,
                      drop = FALSE]
  result[order(result$orthogroup_id, result$species, result$gene_id), ,
         drop = FALSE]
}

comparative_resources_available <- function(repo_root) {
  manifest <- comparative_resource_path(
    repo_root, "manifest", "comparative_species_manifest.yaml"
  )
  gene_lookup <- comparative_resource_path(
    repo_root, "orthology", "gene_lookup.csv"
  )
  orthogroups <- comparative_resource_path(
    repo_root, "orthology", "orthogroups_long.csv"
  )
  threeway <- comparative_resource_path(
    repo_root, "orthology", "qng_gwk_cagl_threeway_map_complete.csv"
  )
  data.frame(
    resource = c("species_manifest", "gene_lookup", "orthogroups", "threeway_map"),
    available = file.exists(c(manifest, gene_lookup, orthogroups, threeway)),
    path = c(manifest, gene_lookup, orthogroups, threeway),
    stringsAsFactors = FALSE
  )
}
