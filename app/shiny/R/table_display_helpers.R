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

ytab_gene_details_datatable <- function(data, page_length = 10L, filter = "none",
                                        options = list(), class = "compact stripe hover",
                                        order = NULL) {
  if (!is.data.frame(data)) data <- data.frame()
  detail_cols <- c(".ytab_gene_detail_name", ".ytab_sgd_description", ".ytab_sgd_essentiality")
  for (col in detail_cols) if (!col %in% names(data)) data[[col]] <- ""
  clickable <- which(names(data) %in% c("CAGL ID", "Gene name")) - 1L
  hidden <- match(detail_cols, names(data)) - 1L
  if (!length(clickable)) {
    return(DT::datatable(data[, setdiff(names(data), detail_cols), drop = FALSE],
                         rownames = FALSE, filter = filter, selection = "none",
                         class = class, options = options))
  }
  column_defs <- options$columnDefs %||% list()
  column_defs <- c(
    column_defs,
    list(
      list(className = "ytab-gene-detail-toggle", targets = clickable),
      list(visible = FALSE, searchable = FALSE, targets = hidden)
    )
  )
  if (!is.null(order)) options$order <- order
  options$columnDefs <- column_defs
  options$pageLength <- options$pageLength %||% page_length
  options$scrollX <- options$scrollX %||% TRUE
  callback <- DT::JS(sprintf(
    "var table = this.api();
     var detailNameCol = %d;
     var detailDescriptionCol = %d;
     var detailEssentialityCol = %d;
     function ytabEscapeGeneDetail(value) {
       if (value === null || value === undefined) return '';
       value = String(value);
       if (!value.trim() || /^(NA|NaN|NULL)$/i.test(value.trim())) return '';
       return $('<div/>').text(value).html();
     }
     function ytabGeneDetailValue(rowData, index, fallback) {
       var value = ytabEscapeGeneDetail(rowData[index]);
       return value || fallback;
     }
     function ytabGeneDetailHtml(rowData) {
       var gene = ytabGeneDetailValue(rowData, detailNameCol, 'Not available');
       var description = ytabGeneDetailValue(rowData, detailDescriptionCol, 'Not available');
       var essentiality = ytabGeneDetailValue(rowData, detailEssentialityCol, 'Not available');
       return '<div class=\"ytab-gene-details-child\">' +
         '<div class=\"ytab-gene-details-title\">Gene details</div>' +
         '<dl>' +
         '<dt>Gene name</dt><dd>' + gene + '</dd>' +
         '<dt>SGD description</dt><dd class=\"ytab-gene-description\">' + description + '</dd>' +
         '<dt>SGD essentiality</dt><dd>' + essentiality + '</dd>' +
         '</dl>' +
         '</div>';
     }
     table.on('click', 'tbody td.ytab-gene-detail-toggle', function() {
       var tr = $(this).closest('tr');
       var row = table.row(tr);
       if (row.child.isShown()) {
         row.child.hide();
         tr.removeClass('ytab-gene-details-open');
       } else {
         row.child(ytabGeneDetailHtml(row.data())).show();
         tr.addClass('ytab-gene-details-open');
       }
     });",
    match(".ytab_gene_detail_name", names(data)) - 1L,
    match(".ytab_sgd_description", names(data)) - 1L,
    match(".ytab_sgd_essentiality", names(data)) - 1L
  ))
  DT::datatable(data, rownames = FALSE, filter = filter, selection = "none",
                class = class, escape = TRUE, callback = callback,
                options = options)
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
