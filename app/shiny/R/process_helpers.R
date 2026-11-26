validate_executable <- function(path) {
  is.character(path) && length(path) == 1L && nzchar(path) && file.exists(path) && !dir.exists(path)
}

locate_python_executable <- function() {
  path <- Sys.which("python")
  if (!nzchar(path)) path <- Sys.which("python3")
  if (!validate_executable(path)) stop("Python was not found on PATH.", call. = FALSE)
  unname(path)
}

locate_rscript_executable <- function() {
  path <- Sys.which("Rscript")
  if (!validate_executable(path)) stop("Rscript was not found on PATH.", call. = FALSE)
  unname(path)
}

validate_process_arguments <- function(command, args) {
  if (!validate_executable(command)) stop("Executable does not exist or is not a regular file: ", command, call. = FALSE)
  if (!is.character(args) || anyNA(args)) stop("Process arguments must be a character vector without missing values.", call. = FALSE)
  forbidden <- c("2>&1", ">", "|", "&&", ";")
  if (any(args %in% forbidden)) stop("Shell operators are not valid process arguments.", call. = FALSE)
  invisible(TRUE)
}

format_command_for_display <- function(command, args = character()) {
  paste(c(shQuote(command), vapply(args, shQuote, character(1))), collapse = " ")
}

run_process_sync <- function(command, args = character(), wd = NULL, env = character()) {
  if (!requireNamespace("processx", quietly = TRUE)) stop("The R package processx is required.", call. = FALSE)
  validate_process_arguments(command, args)
  call_args <- list(command=command,args=args,wd=wd,echo=FALSE,error_on_status=FALSE)
  if(length(env))call_args$env<-env
  do.call(processx::run,call_args)
}

start_process_async <- function(command, args = character(), wd = NULL, env = character()) {
  if (!requireNamespace("processx", quietly = TRUE)) stop("The R package processx is required.", call. = FALSE)
  validate_process_arguments(command, args)
  call_args<-list(command=command,args=args,wd=wd,stdout="|",stderr="|",cleanup=TRUE,cleanup_tree=TRUE)
  if(length(env))call_args$env<-env
  do.call(processx::process$new,call_args)
}

git_commit_or_unavailable <- function(repo_root) {
  git <- Sys.which("git")
  if (!validate_executable(git)) return("unavailable")
  result <- tryCatch(run_process_sync(unname(git), c("-C", repo_root, "rev-parse", "--short", "HEAD")), error = function(e) NULL)
  if (is.null(result) || result$status != 0L || !nzchar(trimws(result$stdout))) "unavailable" else trimws(result$stdout)
}

strip_matching_outer_quotes <- function(path) {
  path <- trimws(path %||% "")
  if (nchar(path) >= 2L) {
    first <- substr(path, 1L, 1L); last <- substr(path, nchar(path), nchar(path))
    if (first == last && first %in% c("'", "\"")) path <- substr(path, 2L, nchar(path) - 1L)
  }
  path
}

normalize_local_directory <- function(path) {
  path <- strip_matching_outer_quotes(path)
  if (!nzchar(path) || !dir.exists(path)) stop("FASTQ directory does not exist: ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# A command-line flag and its value travel together.
#
# Argument vectors are assembled with c(), and c() silently drops NULL and character(0). So
# c("--hmm-run", value) with an empty value yields a bare "--hmm-run", which then swallows whatever flag
# comes next as its value -- the job starts with arguments nobody wrote, instead of failing. Optional flags
# are therefore omitted whole when the user left the field empty; required ones stop with a message naming
# the flag.
ytab_flag <- function(flag, value, required = FALSE) {
  value <- value %||% character(0)
  value <- as.character(value)
  value <- value[!is.na(value) & nzchar(trimws(value))]
  if (!length(value)) {
    if (isTRUE(required)) stop(sprintf("No value for %s.", flag), call. = FALSE)
    return(character(0))
  }
  c(flag, value)
}
