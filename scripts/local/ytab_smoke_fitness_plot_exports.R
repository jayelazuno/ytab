#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]]))
root <- normalizePath(file.path(dirname(self), "../.."))
config <- normalizePath(args[[1]])
env <- new.env(parent=.GlobalEnv)
source(file.path(root, "app/shiny/app.R"), local=env)
shiny::testServer(env$server, {
  session$setInputs(project_selector=config, open_project=1)
  session$flushReact()
  result <- current_fitness_result()
  stopifnot(!is.null(result))
  pairs <- env$fitness_ma_pair_choices(result)
  stopifnot(length(pairs) > 0L)
  for (mode in c("combined", "individual")) {
    session$setInputs(fitness_ma_mode=mode, fitness_ma_pair=unname(pairs[[1]]),
      fitness_ma_hit_direction="both", fitness_ma_top_n=10,
      fitness_ma_annotation_mode="top")
    for (choice in c("combined_ma", "condition_control")) {
      session$setInputs(fitness_visualization_choice=choice)
      session$flushReact()
      data <- fitness_current_plot_data()
      exported <- fitness_current_plot_download_table(data, choice)
      stopifnot(nrow(data) > 0L, nrow(exported) == nrow(data), !anyDuplicated(names(exported)))
      if (choice == "condition_control") {
        stopifnot(isTRUE(all.equal(data$x_log10_control_plus_1, log10(data$control_abundance + 1))),
          isTRUE(all.equal(data$y_log10_treated_plus_1, log10(data$treated_abundance + 1))))
      }
      image <- tempfile(fileext=".png")
      fitness_render_current_plot(image)
      stopifnot(file.info(image)$size > 1000)
      cat("PASS:", mode, choice, "plot and data export\n")
    }
  }
})
