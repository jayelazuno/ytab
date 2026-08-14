fitness_generated_plot_inventory <- function(project) {
  root <- file.path(project$project_root, "treated_vs_parent")
  files <- if (dir.exists(root)) list.files(root, pattern = "\\.png$", recursive = TRUE, full.names = TRUE) else character()
  files <- files[grepl("/plots/", files) & file.exists(files)]
  files <- files[basename(files) != "MA_treated_vs_parent_all_pools.png"]
  if (!length(files)) return(data.frame())
  project_root <- normalizePath(project$project_root, winslash = "/", mustWork = FALSE)
  data.frame(
    file = files,
    filename = basename(files),
    title = tools::file_path_sans_ext(gsub("_", " ", basename(files))),
    served_url = vapply(files, function(path) {
      path <- normalizePath(path, winslash = "/", mustWork = FALSE)
      rel <- substring(path, nchar(project_root) + 2L)
      paste0("ytab-project-output/", diagnostic_encode_relative_path(rel), fitness_plot_cache_token(path))
    }, ""),
    stringsAsFactors = FALSE
  )
}

fitness_plot_cache_token <- function(path) {
  info <- file.info(path)
  if (!nrow(info) || is.na(info$mtime[[1]])) return("")
  paste0("?v=", as.integer(as.POSIXct(info$mtime[[1]])), "-", as.integer(info$size[[1]] %||% 0L))
}

fitness_combined_ma_plot_file <- function(result) {
  if (is.null(result) || !nzchar(result$output_dir %||% "")) return("")
  candidates <- file.path(result$output_dir, "plots",
                          c("MA_treated_vs_parent_combined.png", "MA_treated_vs_parent_all_pools.png"))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) candidates[[1]] else ""
}

fitness_generated_plot_inventory_for_result <- function(project, result) {
  plots <- fitness_generated_plot_inventory(project)
  if (!nrow(plots) || is.null(result) || !nzchar(result$output_dir %||% "")) return(plots)
  ma_candidates <- file.path(result$output_dir, "plots",
                             c("MA_treated_vs_parent_combined.png", "MA_treated_vs_parent_all_pools.png"))
  dedicated_candidates <- c(ma_candidates, fitness_condition_control_plot_file(result))
  dedicated_candidates <- normalizePath(dedicated_candidates[file.exists(dedicated_candidates)], winslash = "/", mustWork = FALSE)
  if (!length(dedicated_candidates)) return(plots)
  plots[!normalizePath(plots$file, winslash = "/", mustWork = FALSE) %in% dedicated_candidates, , drop = FALSE]
}

fitness_project_served_plot <- function(path, project) {
  if (!nzchar(path %||% "") || !file.exists(path)) return("")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  project_root <- normalizePath(project$project_root, winslash = "/", mustWork = FALSE)
  if (!startsWith(path, paste0(project_root, "/"))) return("")
  rel <- substring(path, nchar(project_root) + 2L)
  paste0("ytab-project-output/", diagnostic_encode_relative_path(rel), fitness_plot_cache_token(path))
}

fitness_combined_ma_plot_card <- function(result, project) {
  path <- fitness_combined_ma_plot_file(result)
  served_url <- fitness_project_served_plot(path, project)
  if (!nzchar(served_url)) {
    return(tags$p(class = "text-muted", "The combined treated-versus-parent MA plot is not available for this selected fitness result."))
  }
  tags$article(
    class = "ytab-plot-card",
    tags$h4("Combined treated-versus-parent MA plot"),
    tags$div(
      class = "ytab-diagnostic-preview",
      tags$img(src = served_url, alt = basename(path), loading = "lazy",
               style = "width:100%;max-height:620px;object-fit:contain")
    ),
    tags$details(tags$summary("Show filename"), tags$code(relative_project_path(path, project$project_root)))
  )
}

fitness_condition_control_plot_card <- function(result, project) {
  path <- fitness_condition_control_plot_file(result)
  if (!nzchar(path) || !file.exists(path)) {
    path <- tryCatch(save_fitness_condition_control_scatter(result), error = function(e) "")
  }
  served_url <- fitness_project_served_plot(path, project)
  if (!nzchar(served_url)) {
    return(tags$p(class = "text-muted", "The condition-versus-control log-log scatter plot is not available for this selected fitness result."))
  }
  tags$article(
    class = "ytab-plot-card",
    tags$h4("Condition-versus-control log-log scatter"),
    tags$div(
      class = "ytab-diagnostic-preview",
      tags$img(src = served_url, alt = basename(path), loading = "lazy",
               style = "width:100%;max-height:620px;object-fit:contain")
    ),
    tags$details(tags$summary("Show filename"), tags$code(relative_project_path(path, project$project_root)))
  )
}

fitness_generated_plot_cards <- function(project) {
  plots <- fitness_generated_plot_inventory(project)
  if (!nrow(plots)) return(tags$p(class = "text-muted", "No treated-vs-parent plot files are available."))
  plots <- plots[order(plots$filename), , drop = FALSE]
  tags$div(class = "ytab-plot-grid", lapply(seq_len(nrow(plots)), function(i) {
    tags$article(
      class = "ytab-plot-card",
      title = plots$filename[[i]],
      tags$h4(plots$title[[i]]),
      tags$div(class = "ytab-diagnostic-preview",
        tags$img(src = plots$served_url[[i]], alt = plots$filename[[i]],
                 loading = "lazy", style = "width:100%;height:180px;object-fit:contain")
      ),
      tags$details(tags$summary("Show filename"), tags$code(plots$filename[[i]]))
    )
  }))
}
