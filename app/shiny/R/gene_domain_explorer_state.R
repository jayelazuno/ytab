gene_domain_cli_script <- function(repo_root) {
  file.path(repo_root, "scripts", "local", "ytab_draw_gene_domain_insertions.py")
}

gene_domain_json <- function(python_bin, repo_root, project_config, args = character()) {
  script <- gene_domain_cli_script(repo_root)
  full <- c(script, "--project-config", project_config, args, "--json")
  result <- run_process_sync(python_bin(), full, wd = repo_root)
  parsed <- tryCatch(jsonlite::fromJSON(result$stdout, simplifyVector = FALSE),
                     error = function(e) NULL)
  if (is.null(parsed)) {
    parsed <- list(status = "failure",
                   error = paste(c(result$stdout, result$stderr), collapse = "\n"))
  }
  parsed$exit_status <- as.integer(result$status)
  parsed$stdout <- result$stdout
  parsed$stderr <- result$stderr
  parsed$command <- format_command_for_display(python_bin(), full)
  parsed
}

gene_domain_tracks <- function(python_bin, repo_root, project_config, track_source = "raw") {
  data <- gene_domain_json(
    python_bin, repo_root, project_config,
    c("--list-tracks", "--track-source", track_source)
  )
  tracks <- data$tracks %||% list()
  if (!length(tracks)) return(data.frame())
  data.frame(
    sample = vapply(tracks, function(x) as.character(x$sample %||% ""), ""),
    track_name = vapply(tracks, function(x) as.character(x$track_name %||% x$sample %||% ""), ""),
    track_source = vapply(tracks, function(x) as.character(x$track_source %||% ""), ""),
    source_file = vapply(tracks, function(x) as.character(x$source_file %||% ""), ""),
    stringsAsFactors = FALSE
  )
}

gene_domain_candidate_table <- function(candidates) {
  if (is.null(candidates) || !length(candidates)) return(data.frame())
  rows <- lapply(candidates, function(x) {
    data.frame(
      gene_id = as.character(x$gene_id %||% ""),
      display_name = as.character(x$display_name %||% ""),
      standard_name = as.character(x$standard_name %||% ""),
      common_name = as.character(x$common_name %||% ""),
      systematic_name = as.character(x$systematic_name %||% ""),
      chromosome = as.character(x$chromosome %||% ""),
      start = as.integer(x$start %||% NA_integer_),
      end = as.integer(x$end %||% NA_integer_),
      strand = as.character(x$strand %||% ""),
      product = as.character(x$product %||% ""),
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)
  keep <- vapply(data, function(x) any(nzchar(as.character(x)) | !is.na(x)), logical(1))
  data[, keep, drop = FALSE]
}

gene_domain_manifest_relative_src <- function(manifest, project_root) {
  figure <- as.character(manifest$figure_path %||% "")
  if (!nzchar(figure) || !file.exists(figure)) return("")
  rel <- tryCatch(
    {
      normalized <- normalizePath(figure, winslash = "/", mustWork = TRUE)
      root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
      if (startsWith(normalized, paste0(root, "/"))) substring(normalized, nchar(root) + 2L) else ""
    },
    error = function(e) ""
  )
  if (!nzchar(rel)) return("")
  paste0("ytab-project-output/", rel)
}

gene_domain_read_insertions <- function(manifest) {
  table_path <- as.character(manifest$table_path %||% "")
  if (!nzchar(table_path) || !file.exists(table_path)) return(data.frame())
  tryCatch(read.csv(table_path, stringsAsFactors = FALSE), error = function(e) data.frame())
}

gene_domain_infer_track_role <- function(sample, sample_sheet = data.frame()) {
  row <- if (is.data.frame(sample_sheet) && "sample" %in% names(sample_sheet))
    sample_sheet[as.character(sample_sheet$sample) == sample, , drop = FALSE] else data.frame()
  values <- c()
  for (field in c("library_role", "condition", "treatment", "guessed_condition"))
    if (nrow(row) && field %in% names(row)) values <- c(values, as.character(row[[field]][[1]]))
  values <- values[!is.na(values)]
  text <- tolower(paste(c(values, sample), collapse = " "))
  if (grepl("parent|control", text)) return("parent")
  if (grepl("treated", text)) return("treated")
  ""
}

gene_domain_infer_track_pool <- function(sample, sample_sheet = data.frame()) {
  row <- if (is.data.frame(sample_sheet) && "sample" %in% names(sample_sheet))
    sample_sheet[as.character(sample_sheet$sample) == sample, , drop = FALSE] else data.frame()
  for (field in c("pool", "pool_id", "replicate_id", "guessed_pool")) {
    if (nrow(row) && field %in% names(row)) {
      value <- as.character(row[[field]][[1]])
      if (is.na(value) || !nzchar(value)) next
      hit <- regmatches(value, regexpr("[0-9]+", value))
      if (length(hit) && !is.na(hit) && nzchar(hit)) return(hit)
    }
  }
  hit <- regmatches(sample, regexpr("pool[_ -]?[0-9]+", sample, ignore.case = TRUE))
  if (length(hit) && !is.na(hit) && nzchar(hit)) sub(".*?([0-9]+)$", "\\1", hit) else ""
}

gene_domain_order_tracks <- function(data, sample_sheet = data.frame()) {
  if (!is.data.frame(data) || !nrow(data)) return(data.frame())
  data$role <- vapply(data$sample, gene_domain_infer_track_role, character(1), sample_sheet = sample_sheet)
  data$pool <- vapply(data$sample, gene_domain_infer_track_pool, character(1), sample_sheet = sample_sheet)
  data$order_index <- seq_len(nrow(data))
  if (any(nzchar(data$role))) {
    role_rank <- ifelse(data$role == "parent", 0L, ifelse(data$role == "treated", 1L, 2L))
    pool_rank <- suppressWarnings(as.integer(data$pool))
    pool_rank[is.na(pool_rank)] <- 9999L
    data <- data[order(role_rank, pool_rank, data$order_index), , drop = FALSE]
  }
  data
}

gene_domain_preset_choices <- function(data) {
  if (!all(c("role", "pool") %in% names(data))) data <- gene_domain_order_tracks(data)
  choices <- c("All tracks" = "all")
  if (any(data$role == "parent")) choices <- c(choices, "Parents only" = "parents")
  if (any(data$role == "treated")) choices <- c(choices, "Treated only" = "treated")
  for (pool in as.character(1:4)) {
    if (any(data$role == "parent" & data$pool == pool) &&
        any(data$role == "treated" & data$pool == pool)) {
      choices <- c(choices, setNames(paste0("pool", pool, "_pair"),
                                     paste0("Parent pool ", pool, " vs treated pool ", pool)))
    }
  }
  c(choices, "Custom" = "custom")
}

gene_domain_preset_track_rows <- function(preset, data, current_selection = character()) {
  if (!all(c("role", "pool") %in% names(data))) data <- gene_domain_order_tracks(data)
  if (!nrow(data)) return(data.frame())
  preset <- preset %||% "all"
  if (identical(preset, "custom")) return(data)
  if (identical(preset, "all")) return(data)
  if (identical(preset, "parents")) return(data[data$role == "parent", , drop = FALSE])
  if (identical(preset, "treated")) return(data[data$role == "treated", , drop = FALSE])
  pool <- sub("^pool([1-4])_pair$", "\\1", preset)
  if (identical(pool, preset)) return(data.frame())
  data[(data$role == "parent" & data$pool == pool) |
         (data$role == "treated" & data$pool == pool), , drop = FALSE]
}

gene_domain_valid_output <- function(state) {
  manifest <- state$manifest %||% NULL
  if (is.null(manifest)) return(FALSE)
  figure <- as.character(manifest$figure_path %||% "")
  table <- as.character(manifest$table_path %||% "")
  nzchar(figure) && file.exists(figure) && nzchar(table) && file.exists(table)
}
