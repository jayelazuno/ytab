#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(
  sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1L]]),
  winslash = "/"
)
root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/")
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
source(file.path(root, "app/shiny/R/project_discovery.R"))
source(file.path(root, "app/shiny/R/essentiality_targets.R"))
source(file.path(root, "app/shiny/R/essentiality_results.R"))

project <- read_project_summary(args[[1L]], root)
stored <- essentiality_target_summary_data(project$project_root)
if (!nrow(stored)) stored <- data.frame(
  target = 6400, target_tag = "T6400",
  min_hit_site_retention_fraction = .95,
  min_feature_retention_fraction = .95,
  stringsAsFactors = FALSE
)
one <- stored[1L, , drop = FALSE]
two <- rbind(one, one)
three <- rbind(one, one, one)

zero_state <- essentiality_target_visual_state(data.frame())
one_state <- essentiality_target_visual_state(one)
two_state <- essentiality_target_visual_state(two)
three_state <- essentiality_target_visual_state(three)
stopifnot(
  identical(zero_state$mode, "empty"), !zero_state$show_table,
  identical(one_state$mode, "single"), !one_state$show_table,
  !one_state$show_site_plot, !one_state$show_feature_plot,
  grepl("Only one target", one_state$message, fixed = TRUE),
  identical(two_state$mode, "compact"), two_state$show_table,
  two_state$show_combined_plot, !two_state$show_site_plot,
  identical(three_state$mode, "trends"), three_state$show_table,
  three_state$show_site_plot, three_state$show_feature_plot
)

objects <- list(
  missing = NULL,
  data_frame = data.frame(
    target = 6400, target_tag = "T6400",
    recommendation_type = "Feature-evaluated", stringsAsFactors = FALSE
  ),
  named_list = list(
    target = 6400, target_tag = "T6400",
    recommendation_type = "Feature-evaluated"
  ),
  legacy_vector = c(
    target = "6400", target_tag = "T6400",
    recommendation_type = "Feature-evaluated"
  ),
  malformed = c("not", "named", "recommendation")
)
normalized <- lapply(objects, function(object)
  tryCatch(normalize_recommendation_state(object), error = identity))
stopifnot(!any(vapply(normalized, inherits, logical(1), "error")))
required <- c(
  "target", "target_tag", "recommendation_type", "site_retention",
  "feature_retention", "parents_passing", "reason"
)
for (item in normalized) stopifnot(all(required %in% names(item)))
stopifnot(!normalized$missing$available,
          normalized$data_frame$available,
          normalized$named_list$available,
          normalized$legacy_vector$available,
          !normalized$malformed$available)
evaluation_root<-file.path(tempdir(),"ytab-evaluation-state");base<-file.path(evaluation_root,"sample_normalization");dir.create(base,recursive=TRUE,showWarnings=FALSE);jsonlite::write_json(list(target=6400,target_tag="T6400"),file.path(base,"normalization_recommendation.json"),auto_unbox=TRUE);stopifnot(essentiality_feature_evaluation_state(evaluation_root)=="preliminary_recommendation_only");dir.create(file.path(base,"T6400","parent"),recursive=TRUE,showWarnings=FALSE);writeLines("x",file.path(base,"T6400","parent","parent_normalized_hits.txt"));stopifnot(essentiality_feature_evaluation_state(evaluation_root)=="feature_evaluation_ready");write.csv(data.frame(target=6400,target_tag="T6400",min_feature_retention_fraction=.95),file.path(base,"normalization_target_evaluation.csv"),row.names=FALSE);stopifnot(essentiality_feature_evaluation_state(evaluation_root)=="feature_evaluation_complete")
cat("PASS\n")
