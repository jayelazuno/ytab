ytab_datatable <- function(data, page_length = 10L, filter = "top",
                           scroll_x = TRUE, compact = TRUE) {
  classes <- if (compact) "compact stripe hover" else "stripe hover"
  DT::datatable(
    data,
    rownames = FALSE,
    filter = filter,
    class = classes,
    selection = "none",
    options = list(pageLength = page_length, scrollX = scroll_x)
  )
}

ytab_file_label <- function(path) {
  path <- as.character(path %||% "")
  if (!nzchar(path)) return("")
  basename(path)
}

ytab_file_table <- function(paths) {
  paths <- as.character(paths %||% character())
  paths <- paths[nzchar(paths)]
  if (!length(paths)) return(data.frame())
  info <- file.info(paths)
  data.frame(
    File = basename(paths),
    Type = tools::file_ext(paths),
    Size = ifelse(is.na(info$size), "", human_file_size(info$size)),
    Modified = ifelse(is.na(info$mtime), "", format(info$mtime, "%Y-%m-%d %H:%M")),
    Path = paths,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

