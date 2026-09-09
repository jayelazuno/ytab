#!/usr/bin/env Rscript

`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      i <- i + 1L
      next
    }
    key <- sub("^--", "", key)
    if (key %in% c("semantic-reduction")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      if (i == length(args)) stop("Missing value for --", key, call. = FALSE)
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg, call. = FALSE)
  }
}

clean_pvalue <- function(x) {
  x <- as.character(x)
  x <- gsub("<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

split_orfs <- function(x) {
  x <- as.character(x %||% "")
  x <- unlist(strsplit(x, "[;,[:space:]]+"), use.names = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x) & grepl("^Y[A-Z]{2}[0-9]{3}[A-Z](-[A-Z])?$", x, ignore.case = TRUE)]
  unique(toupper(x))
}

expand_orf_table <- function(data) {
  rows <- lapply(seq_len(nrow(data)), function(i) {
    orfs <- split_orfs(data$scer_gene_id[[i]])
    if (!length(orfs)) return(NULL)
    data.frame(
      input_gene = as.character(data$input_gene[[i]] %||% ""),
      cagl_id = as.character(data$cagl_id[[i]] %||% ""),
      gene_name = as.character(data$gene_name[[i]] %||% ""),
      scer_gene_id = orfs,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(data.frame(input_gene = character(), cagl_id = character(),
                      gene_name = character(), scer_gene_id = character(),
                      stringsAsFactors = FALSE))
  }
  unique(out)
}

build_gene2go <- function(universe_orfs, ontology_code) {
  require_package("AnnotationDbi")
  require_package("org.Sc.sgd.db")
  universe_orfs <- unique(universe_orfs[nzchar(universe_orfs)])
  if (!length(universe_orfs)) stop("No S. cerevisiae ORFs are available for the GO universe.", call. = FALSE)
  go_map <- AnnotationDbi::select(
    org.Sc.sgd.db::org.Sc.sgd.db,
    keys = universe_orfs,
    columns = c("GOALL", "ONTOLOGYALL"),
    keytype = "ORF"
  )
  go_map <- go_map[!is.na(go_map$GOALL) & go_map$ONTOLOGYALL == ontology_code, c("ORF", "GOALL"), drop = FALSE]
  go_map <- unique(go_map)
  gene2go <- split(go_map$GOALL, go_map$ORF)
  lapply(gene2go, unique)
}

run_topgo_one <- function(target_orfs, universe_orfs, gene2go, ontology_code,
                          ontology_label, node_size, p_cutoff) {
  require_package("topGO")
  target_orfs <- unique(target_orfs[nzchar(target_orfs)])
  universe_orfs <- unique(universe_orfs[nzchar(universe_orfs)])
  all_genes <- universe_orfs[universe_orfs %in% names(gene2go)]
  target_annotated <- target_orfs[target_orfs %in% all_genes]
  if (!length(all_genes) || !length(target_annotated)) {
    return(data.frame(
      GO.ID = character(), Term = character(), ontology_code = character(),
      ontology_label = character(), classic = character(), weight01 = character(),
      classic_numeric = numeric(), weight01_numeric = numeric(),
      target_genes_total = integer(), target_genes_with_GO = integer(),
      universe_genes_with_GO = integer(), node_size = integer(),
      p_cutoff = numeric(), stringsAsFactors = FALSE
    ))
  }
  gene_list <- factor(as.integer(all_genes %in% target_annotated))
  names(gene_list) <- all_genes
  go_data <- methods::new(
    "topGOdata",
    ontology = ontology_code,
    allGenes = gene_list,
    annot = annFUN.gene2GO,
    gene2GO = gene2go,
    nodeSize = node_size
  )
  if (!length(usedGO(go_data))) {
    return(data.frame())
  }
  result_weight01 <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
  result_classic <- runTest(go_data, algorithm = "classic", statistic = "fisher")
  res <- GenTable(
    go_data,
    classic = result_classic,
    weight01 = result_weight01,
    orderBy = "weight01",
    ranksOf = "weight01",
    topNodes = length(usedGO(go_data)),
    numChar = 1000
  )
  res$classic_numeric <- clean_pvalue(res$classic)
  res$weight01_numeric <- clean_pvalue(res$weight01)
  res$ontology_code <- ontology_code
  res$ontology_label <- ontology_label
  res$target_genes_total <- length(target_orfs)
  res$target_genes_with_GO <- length(target_annotated)
  res$universe_genes_with_GO <- length(all_genes)
  res$node_size <- node_size
  res$p_cutoff <- p_cutoff
  res[order(res$weight01_numeric, res$classic_numeric), , drop = FALSE]
}

semantic_reduce_one <- function(go_df, ontology_code, threshold = 0.70,
                                method = "Rel") {
  require_package("rrvgo")
  require_package("GOSemSim")
  require_package("org.Sc.sgd.db")
  ont_df <- go_df[
    go_df$ontology_code == ontology_code &
      is.finite(go_df$weight01_numeric) &
      go_df$weight01_numeric > 0,
    ,
    drop = FALSE
  ]
  if (nrow(ont_df) < 2L) return(data.frame())
  scores_tbl <- aggregate(
    -log10(pmax(ont_df$weight01_numeric, 1e-300)),
    list(GO.ID = ont_df$GO.ID),
    max
  )
  names(scores_tbl)[2L] <- "score"
  terms <- ont_df[match(scores_tbl$GO.ID, ont_df$GO.ID), c("GO.ID", "Term", "ontology_label"), drop = FALSE]
  scores_tbl <- merge(scores_tbl, terms, by = "GO.ID", all.x = TRUE, sort = FALSE)
  scores <- stats::setNames(scores_tbl$score, scores_tbl$GO.ID)
  semdata <- GOSemSim::godata(annoDb = "org.Sc.sgd.db", ont = ontology_code, computeIC = TRUE)
  sim_matrix <- rrvgo::calculateSimMatrix(
    scores_tbl$GO.ID,
    orgdb = "org.Sc.sgd.db",
    ont = ontology_code,
    method = method,
    semdata = semdata
  )
  usable <- intersect(rownames(sim_matrix), names(scores))
  sim_matrix <- sim_matrix[usable, usable, drop = FALSE]
  scores <- scores[usable]
  if (nrow(sim_matrix) < 2L) return(data.frame())
  reduced <- rrvgo::reduceSimMatrix(
    sim_matrix,
    scores = scores,
    threshold = threshold,
    orgdb = "org.Sc.sgd.db"
  )
  reduced <- as.data.frame(reduced, stringsAsFactors = FALSE)
  reduced$ontology_code <- ontology_code
  reduced$ontology_label <- unique(ont_df$ontology_label)[[1L]]
  reduced$semantic_similarity_threshold <- threshold
  reduced$semantic_similarity_method <- method
  reduced
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
targets_path <- args[["targets"]] %||% ""
universe_path <- args[["universe"]] %||% ""
out_dir <- args[["out-dir"]] %||% ""
ontologies <- strsplit(args[["ontologies"]] %||% "BP,MF,CC", "[,;[:space:]]+")[[1]]
ontologies <- intersect(toupper(ontologies[nzchar(ontologies)]), c("BP", "MF", "CC"))
node_size <- as.integer(args[["node-size"]] %||% "10")
p_cutoff <- as.numeric(args[["p-cutoff"]] %||% "0.05")
semantic_reduction <- isTRUE(args[["semantic-reduction"]])
rrvgo_threshold <- as.numeric(args[["rrvgo-threshold"]] %||% "0.70")

for (pkg in c("topGO", "GO.db", "org.Sc.sgd.db", "AnnotationDbi")) require_package(pkg)
if (semantic_reduction) for (pkg in c("rrvgo", "GOSemSim")) require_package(pkg)
suppressPackageStartupMessages({
  library(topGO)
  library(AnnotationDbi)
  library(org.Sc.sgd.db)
})
if (!file.exists(targets_path)) stop("Target gene file is missing: ", targets_path, call. = FALSE)
if (!file.exists(universe_path)) stop("Universe gene file is missing: ", universe_path, call. = FALSE)
if (!nzchar(out_dir)) stop("--out-dir is required.", call. = FALSE)

tables_dir <- file.path(out_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

targets <- read.csv(targets_path, stringsAsFactors = FALSE, check.names = FALSE)
universe <- read.csv(universe_path, stringsAsFactors = FALSE, check.names = FALSE)
targets_expanded <- expand_orf_table(targets)
universe_expanded <- expand_orf_table(universe)
target_orfs <- unique(targets_expanded$scer_gene_id)
universe_orfs <- unique(universe_expanded$scer_gene_id)
universe_orfs <- union(universe_orfs, target_orfs)

write.csv(targets_expanded, file.path(tables_dir, "target_genes_expanded_scer_orfs.csv"), row.names = FALSE)
write.csv(universe_expanded, file.path(tables_dir, "universe_genes_expanded_scer_orfs.csv"), row.names = FALSE)

ontology_labels <- c(BP = "Biological Process", MF = "Molecular Function", CC = "Cellular Component")
topgo_results <- list()
for (ont in ontologies) {
  gene2go <- build_gene2go(universe_orfs, ont)
  topgo_results[[ont]] <- run_topgo_one(
    target_orfs = target_orfs,
    universe_orfs = universe_orfs,
    gene2go = gene2go,
    ontology_code = ont,
    ontology_label = unname(ontology_labels[[ont]]),
    node_size = node_size,
    p_cutoff = p_cutoff
  )
}
topgo_combined <- do.call(rbind, topgo_results)
if (is.null(topgo_combined) || !nrow(topgo_combined)) {
  topgo_combined <- data.frame(
    GO.ID = character(), Term = character(), ontology_code = character(),
    ontology_label = character(), weight01_numeric = numeric(),
    classic_numeric = numeric(), stringsAsFactors = FALSE
  )
}
write.csv(topgo_combined, file.path(tables_dir, "topgo_results.csv"), row.names = FALSE)
sig <- topgo_combined[
  is.finite(topgo_combined$weight01_numeric) & topgo_combined$weight01_numeric < p_cutoff,
  ,
  drop = FALSE
]
write.csv(sig, file.path(tables_dir, "topgo_significant.csv"), row.names = FALSE)

if (semantic_reduction) {
  reduced_list <- lapply(ontologies, function(ont) semantic_reduce_one(sig, ont, rrvgo_threshold, "Rel"))
  reduced <- do.call(rbind, reduced_list)
  if (is.null(reduced) || !nrow(reduced)) {
    reduced <- data.frame()
  }
  write.csv(reduced, file.path(tables_dir, "rrvgo_reduced_terms.csv"), row.names = FALSE)
  if (nrow(reduced) && all(c("parent", "parentTerm", "score", "ontology_code") %in% names(reduced))) {
    hm <- aggregate(
      reduced$score,
      list(
        representative_GO_ID = reduced$parent,
        representative_term = reduced$parentTerm,
        ontology_code = reduced$ontology_code
      ),
      max
    )
    names(hm)[4L] <- "score"
  } else {
    hm <- data.frame(representative_GO_ID = character(), representative_term = character(),
                     ontology_code = character(), score = numeric(), stringsAsFactors = FALSE)
  }
  write.csv(hm, file.path(tables_dir, "semantic_reduced_heatmap.csv"), row.names = FALSE)
} else {
  write.csv(data.frame(), file.path(tables_dir, "rrvgo_reduced_terms.csv"), row.names = FALSE)
  write.csv(data.frame(), file.path(tables_dir, "semantic_reduced_heatmap.csv"), row.names = FALSE)
}

summary <- data.frame(
  ontology_code = ontologies,
  ontology_label = unname(ontology_labels[ontologies]),
  target_orfs = length(target_orfs),
  universe_orfs = length(universe_orfs),
  significant_terms = vapply(ontologies, function(ont) {
    sum(sig$ontology_code == ont, na.rm = TRUE)
  }, integer(1)),
  stringsAsFactors = FALSE
)
write.csv(summary, file.path(tables_dir, "go_analysis_summary.csv"), row.names = FALSE)
message("GO analysis complete. Tables: ", tables_dir)
