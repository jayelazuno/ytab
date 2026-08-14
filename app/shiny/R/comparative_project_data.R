comparative_empty_projects <- function() {
  data.frame(
    species = character(), project_id = character(), display_name = character(),
    project_config = character(), project_type = character(),
    has_classifier = logical(), has_fitness = logical(),
    has_summary_stats = logical(), stringsAsFactors = FALSE
  )
}

comparative_species_key <- function(value) {
  text <- tolower(trimws(as.character(value %||% "")))
  text <- gsub("[ .-]+", "_", text)
  aliases <- c(
    cg = "glabrata", cgla = "glabrata", c_glabrata = "glabrata",
    candida_glabrata = "glabrata", glabrata = "glabrata",
    ca = "albicans", calb = "albicans", c_albicans = "albicans",
    candida_albicans = "albicans", albicans = "albicans",
    sc = "cerevisiae", scer = "cerevisiae", s_cerevisiae = "cerevisiae",
    saccharomyces_cerevisiae = "cerevisiae", cerevisiae = "cerevisiae",
    kl = "lactis", klac = "lactis", k_lactis = "lactis",
    kluyveromyces_lactis = "lactis", lactis = "lactis"
  )
  unname(ifelse(text %in% names(aliases), aliases[text], text))
}

infer_repo_root_from_project_config <- function(project_config) {
  path <- normalizePath(project_config, winslash = "/", mustWork = FALSE)
  if (basename(path) == "project.yaml" && basename(dirname(path)) == "config")
    return(normalizePath(file.path(dirname(path), "../../../.."),
                         winslash = "/", mustWork = FALSE))
  normalizePath(".", winslash = "/", mustWork = FALSE)
}

comparative_project_summary_stats_path <- function(project_root) {
  path <- file.path(project_root, "summary", "summary_stats.all_samples.csv")
  if (file.exists(path)) path else ""
}

list_ytab_projects_by_species <- function(repo_root) {
  projects <- tryCatch(discover_ytab_projects(repo_root, include_tests = TRUE),
                       error = function(e) data.frame())
  if (!nrow(projects)) return(comparative_empty_projects())
  rows <- lapply(seq_len(nrow(projects)), function(index) {
    project <- tryCatch(read_project_summary(projects$project_config[[index]],
                                             repo_root),
                        error = function(e) NULL)
    if (is.null(project)) return(NULL)
    classifiers <- tryCatch(discover_classifier_results(project$project_root),
                            error = function(e) list())
    fitness <- tryCatch(discover_fitness_results(project),
                        error = function(e) list())
    data.frame(
      species = comparative_species_key(project$species),
      project_id = project$project_id,
      display_name = project$display_name,
      project_config = project$project_config,
      project_type = project$project_type,
      has_classifier = length(classifiers) > 0L,
      has_fitness = length(fitness) > 0L,
      has_summary_stats = nzchar(comparative_project_summary_stats_path(project$project_root)),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) comparative_empty_projects() else do.call(rbind, rows)
}

load_project_classifier_results <- function(project_config) {
  repo_root <- infer_repo_root_from_project_config(project_config)
  project <- tryCatch(read_project_summary(project_config, repo_root),
                      error = function(e) NULL)
  if (is.null(project)) return(data.frame())
  results <- tryCatch(discover_classifier_results(project$project_root),
                      error = function(e) list())
  rows <- lapply(results, function(result) {
    data <- read_essentiality_predictions(result$path)
    if (!nrow(data)) return(NULL)
    data$ytab_project_id <- project$project_id
    data$ytab_species <- comparative_species_key(project$species)
    data$ytab_result_type <- "classifier_predictions"
    data$ytab_target_tag <- result$tag
    data$ytab_source_path <- result$path
    data
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

load_project_fitness_results <- function(project_config) {
  repo_root <- infer_repo_root_from_project_config(project_config)
  project <- tryCatch(read_project_summary(project_config, repo_root),
                      error = function(e) NULL)
  if (is.null(project)) return(data.frame())
  results <- tryCatch(discover_fitness_results(project), error = function(e) list())
  rows <- lapply(results, function(result) {
    data <- fitness_result_data(result)
    if (!nrow(data)) return(NULL)
    data$ytab_project_id <- project$project_id
    data$ytab_species <- comparative_species_key(project$species)
    data$ytab_result_type <- "treated_vs_parent_fitness"
    data$ytab_analysis_id <- result$analysis_id
    data$ytab_source_path <- result$table
    data
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

load_project_summary_stats <- function(project_config) {
  repo_root <- infer_repo_root_from_project_config(project_config)
  project <- tryCatch(read_project_summary(project_config, repo_root),
                      error = function(e) NULL)
  if (is.null(project)) return(data.frame())
  path <- comparative_project_summary_stats_path(project$project_root)
  if (!nzchar(path)) return(data.frame())
  data <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                   error = function(e) data.frame())
  if (!nrow(data)) return(data.frame())
  data$ytab_project_id <- project$project_id
  data$ytab_species <- comparative_species_key(project$species)
  data$ytab_result_type <- "raw_summary_stats"
  data$ytab_source_path <- path
  data
}

summarize_species_project_availability <- function(repo_root) {
  species <- load_comparative_species_manifest(repo_root)
  projects <- list_ytab_projects_by_species(repo_root)
  rows <- lapply(seq_len(nrow(species)), function(index) {
    key <- species$species[[index]]
    subset <- projects[projects$species == key, , drop = FALSE]
    data.frame(
      species = key,
      label = species$label[[index]],
      enabled = species$enabled[[index]],
      placeholder = species$placeholder[[index]],
      project_count = nrow(subset),
      classifier_projects = if (nrow(subset)) sum(subset$has_classifier) else 0L,
      fitness_projects = if (nrow(subset)) sum(subset$has_fitness) else 0L,
      summary_projects = if (nrow(subset)) sum(subset$has_summary_stats) else 0L,
      status = if (isTRUE(species$placeholder[[index]])) "Placeholder" else
        if (nrow(subset)) "YTAB outputs available" else
          "No YTAB project outputs are available for this species yet.",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
