# Reading the legacy classifier's gene predictions.
#
# The classifier itself is retired: its pipeline still runs (`ytab classifier`) and its results are still
# compared in the manuscript, but the app's classifier tab was dissolved by the module reorganisation and the
# live Essentiality tab is the window HMM and the hybrid. What survives here is the little that live code still
# needs -- the Across Species view, which pools predictions across projects, and the fitness run control, which
# offers a classifier result to annotate calls with.
#
# A run is named after whatever tag the pipeline was given, which is the normalization run that produced its
# feature table (betageom__norm_...). It used to be a numeric target (T030) from the legacy MidLC thinning;
# nothing writes those any more, and filtering for them hid every result on disk.

read_essentiality_predictions <- function(path) {
  if (!nzchar(path %||% "") || !file.exists(path) || dir.exists(path)) return(data.frame())
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) data.frame())
}

# Every classifier run under a project: its tag and the predictions file, newest first.
discover_classifier_results <- function(project_root) {
  root <- file.path(project_root, "classifier")
  if (!dir.exists(root)) return(list())
  files <- list.files(root, pattern = "^essentiality_predictions\\..+\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(list())
  results <- lapply(files, function(path) {
    parent <- basename(dirname(path))
    # The directory names the run; fall back to the tag inside the file name for a flat layout.
    tag <- if (nzchar(parent) && !identical(parent, "classifier")) parent else
      sub("^essentiality_predictions\\.(.+)\\.csv$", "\\1", basename(path))
    if (!nzchar(tag) || identical(tag, basename(path))) return(NULL)
    list(tag = tag, path = path, modified = file.info(path)$mtime)
  })
  results <- Filter(Negate(is.null), results)
  if (!length(results)) return(list())
  # One result per tag: a tag written twice keeps the newer file.
  results <- lapply(split(results, vapply(results, `[[`, "", "tag")), function(items)
    items[[which.max(vapply(items, function(x) as.numeric(x$modified), numeric(1)))]])
  unname(results[order(-vapply(results, function(x) as.numeric(x$modified), numeric(1)))])
}
