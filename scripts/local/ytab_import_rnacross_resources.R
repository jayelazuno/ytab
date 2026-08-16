#!/usr/bin/env Rscript

`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left

parse_args <- function(args) {
  reference_name <- paste0("RNA", "cross")
  out <- list(rnacross_dir = paste0("codex/", reference_name),
              out_dir = "resources/comparative")
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--rnacross-dir", "--out-dir")) {
      if (i == length(args)) stop("Missing value for ", key, call. = FALSE)
      out[[sub("^--", "", gsub("-", "_", key))]] <- args[[i + 1L]]
      i <- i + 2L
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
  }
  out
}

script_path <- function() {
  all_args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", all_args, value = TRUE)
  if (length(hit)) normalizePath(gsub("~\\+~", " ", sub("^--file=", "", hit[[1L]])),
                                 winslash = "/", mustWork = FALSE)
  else normalizePath("scripts/local/ytab_import_rnacross_resources.R",
                     winslash = "/", mustWork = FALSE)
}

repo_root <- normalizePath(file.path(dirname(script_path()), "../.."),
                           winslash = "/", mustWork = TRUE)

resolve_path <- function(path) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) path else file.path(repo_root, path)
}

resolve_rnacross_dir <- function(path) {
  reference_name <- paste0("RNA", "cross")
  candidate <- resolve_path(path)
  if (dir.exists(candidate)) return(normalizePath(candidate, winslash = "/"))
  fallback <- file.path(repo_root, "docs", "codex", reference_name)
  if (identical(gsub("\\\\", "/", path), paste0("codex/", reference_name)) &&
      dir.exists(fallback))
    return(normalizePath(fallback, winslash = "/"))
  candidate
}

normalize_species <- function(value) {
  text <- tolower(trimws(as.character(value %||% "")))
  text <- gsub("[ .-]+", "_", text)
  dplyr_free <- c(
    cg = "glabrata", cgla = "glabrata", c_glabrata = "glabrata",
    candida_glabrata = "glabrata", glabrata = "glabrata",
    ca = "albicans", calb = "albicans", c_albicans = "albicans",
    candida_albicans = "albicans", albicans = "albicans",
    sc = "cerevisiae", scer = "cerevisiae", s_cerevisiae = "cerevisiae",
    saccharomyces_cerevisiae = "cerevisiae", cerevisiae = "cerevisiae",
    kl = "lactis", klac = "lactis", k_lactis = "lactis",
    kluyveromyces_lactis = "lactis", lactis = "lactis"
  )
  unname(ifelse(text %in% names(dplyr_free), dplyr_free[text], text))
}

safe_read_csv <- function(path, warnings) {
  if (!file.exists(path)) {
    warnings <<- c(warnings, paste("Missing source file:", path))
    return(data.frame())
  }
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) {
             warnings <<- c(warnings, paste("Could not read", path, conditionMessage(e)))
             data.frame()
           })
}

safe_read_rds <- function(path, warnings) {
  if (!file.exists(path)) {
    warnings <<- c(warnings, paste("Missing source file:", path))
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) {
    warnings <<- c(warnings, paste("Could not read", path, conditionMessage(e)))
    NULL
  })
}

write_csv_record <- function(data, source, output, manifest, warnings) {
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  if (!is.data.frame(data)) data <- data.frame()
  write.csv(data, output, row.names = FALSE, na = "")
  source_path <- normalizePath(source, winslash = "/", mustWork = FALSE)
  source_root <- paste0(normalizePath(rnacross_dir, winslash = "/",
                                      mustWork = FALSE), "/")
  source_record <- if (startsWith(source_path, source_root))
    substring(source_path, nchar(source_root) + 1L) else basename(source_path)
  manifest$resources[[length(manifest$resources) + 1L]] <- list(
    source_file = source_record,
    output_file = normalizePath(output, winslash = "/", mustWork = FALSE),
    row_count = nrow(data),
    column_names = names(data),
    warnings = unname(warnings)
  )
  manifest
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
rnacross_dir <- resolve_rnacross_dir(args$rnacross_dir)
out_dir <- resolve_path(args$out_dir)
if (!dir.exists(rnacross_dir))
  stop(paste0("RNA", "cross"), " source directory does not exist: ", rnacross_dir, call. = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE))
  stop("The jsonlite package is required.", call. = FALSE)

warnings <- character()
manifest <- list(
  import_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  rnacross_source = paste(paste0("RNA", "cross"), "reference clone"),
  out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
  resources = list(),
  warnings = list()
)

gene_csv <- file.path(rnacross_dir, "features", "gene_lookup_table", "gene_lookup.csv")
gene_rds <- file.path(rnacross_dir, "features", "gene_lookup_table", "gene_lookup.rds")
gene_lookup <- safe_read_csv(gene_csv, warnings)
gene_source <- gene_csv
if (!nrow(gene_lookup)) {
  rds <- safe_read_rds(gene_rds, warnings)
  if (is.data.frame(rds)) {
    gene_lookup <- as.data.frame(rds, stringsAsFactors = FALSE)
    gene_source <- gene_rds
  }
}
if (nrow(gene_lookup)) {
  if ("species" %in% names(gene_lookup))
    gene_lookup$species <- normalize_species(gene_lookup$species)
  required <- c("gene_id", "species", "expression_id", "id_type",
                "source_info", "gene_name", "hog_id", "og_id")
  for (name in setdiff(required, names(gene_lookup))) gene_lookup[[name]] <- ""
  gene_lookup <- gene_lookup[, required, drop = FALSE]
}
manifest <- write_csv_record(
  gene_lookup, gene_source,
  file.path(out_dir, "orthology", "gene_lookup.csv"),
  manifest, warnings
)

orthogroups_rds <- file.path(rnacross_dir, "features", "orthology", "orthogroups.rds")
orthogroups <- safe_read_rds(orthogroups_rds, warnings)
orthogroups_long <- data.frame(
  orthogroup_id = character(), species = character(), gene_id = character(),
  gene_name = character(), source_column = character(), stringsAsFactors = FALSE
)
if (is.data.frame(orthogroups)) {
  orthogroups <- as.data.frame(orthogroups, stringsAsFactors = FALSE)
  if (all(c("gene_id", "og_id") %in% names(orthogroups))) {
    lookup <- gene_lookup
    lookup <- lookup[!duplicated(lookup$gene_id), , drop = FALSE]
    joined <- merge(orthogroups, lookup[, c("gene_id", "species", "gene_name"),
                                        drop = FALSE],
                    by = "gene_id", all.x = TRUE, sort = FALSE)
    joined$species <- normalize_species(joined$species)
    orthogroups_long <- data.frame(
      orthogroup_id = as.character(joined$og_id),
      species = as.character(joined$species %||% ""),
      gene_id = as.character(joined$gene_id),
      gene_name = as.character(joined$gene_name %||% ""),
      source_column = "orthogroups.rds:og_id",
      stringsAsFactors = FALSE
    )
  } else {
    warnings <- c(warnings, "orthogroups.rds did not contain gene_id and og_id columns.")
  }
} else {
  warnings <- c(warnings, "orthogroups.rds was not a data frame and was not converted.")
}
manifest <- write_csv_record(
  orthogroups_long, orthogroups_rds,
  file.path(out_dir, "orthology", "orthogroups_long.csv"),
  manifest, warnings
)

threeway_source <- file.path(
  rnacross_dir, "features", "ID_maps", "qng_gwk_cagl_threeway_map_complete.csv"
)
threeway <- safe_read_csv(threeway_source, warnings)
manifest <- write_csv_record(
  threeway, threeway_source,
  file.path(out_dir, "orthology", "qng_gwk_cagl_threeway_map_complete.csv"),
  manifest, warnings
)

annotation_paths <- c(
  list.files(file.path(rnacross_dir, "features", "Annotations", "annotations_rds"),
             full.names = TRUE, pattern = "\\.rds$"),
  list.files(file.path(rnacross_dir, "features", "Annotations", "gff3"),
             full.names = TRUE)
)
annotation_root <- paste0(normalizePath(rnacross_dir, winslash = "/",
                                        mustWork = FALSE), "/")
annotation_records <- normalizePath(annotation_paths, winslash = "/",
                                    mustWork = FALSE)
manifest$annotations_inspected <- ifelse(
  startsWith(annotation_records, annotation_root),
  substring(annotation_records, nchar(annotation_root) + 1L),
  basename(annotation_records)
)
manifest$warnings <- as.list(unique(warnings))

manifest_path <- file.path(out_dir, "manifest", "resource_import_manifest.json")
dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)

cat("Imported ", paste0("RNA", "cross"), " comparative resources\n", sep = "")
cat("Source:", rnacross_dir, "\n")
cat("Output:", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n")
if (length(unique(warnings))) {
  cat("Warnings:\n")
  for (warning in unique(warnings)) cat("-", warning, "\n")
}
