# Turning a plotly figure into a file, outside a browser session.
#
# A plotly figure is drawn by plotly.js, so any static image of one needs a browser engine. The app's download
# buttons get this for free -- the figure is already in the user's browser and Plotly.toImage renders it there.
# A command line has no browser, so it borrows one: chromote drives the Chrome already installed on the machine
# and runs the same Plotly.toImage call against the same figure.
#
# The point of going through plotly rather than redrawing in ggplot is that there is then one definition of
# each figure. Two drawings of the same thing drift -- panels, scales, names and whole layers went missing
# between the app's plots and their downloads for exactly that reason -- and no test catches it reliably.
#
# Needs chromote and a Chrome/Chromium on PATH or in the usual place; ytab_capture_available() reports whether
# this machine can do it, so callers can say so plainly instead of failing late.

ytab_capture_available <- function() {
  if (!requireNamespace("chromote", quietly = TRUE)) return(FALSE)
  !is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL))
}

# `widget` is a plotly object; `file` ends in .png or .svg. Width and height are pixels of the rendered figure,
# and `scale` multiplies the raster resolution (3 is roughly 300 dpi at 12 inches).
ytab_capture_plotly <- function(widget, file, width = 1600, height = 900, scale = 3, timeout = 180) {
  if (!ytab_capture_available())
    stop("Writing an image of a plotly figure needs the chromote package and a Chrome or Chromium installation.", call. = FALSE)
  format <- if (grepl("\\.svg$", file, ignore.case = TRUE)) "svg" else "png"
  dir <- file.path(tempdir(), paste0("ytab_figure_", as.integer(Sys.time()), "_", sample.int(1e6, 1)))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  page <- file.path(dir, "figure.html")
  # selfcontained = TRUE routes the whole widget through pandoc, which on a figure with thousands of segments
  # takes longer than the render itself. The libraries sit beside the page instead; the directory is temporary.
  htmlwidgets::saveWidget(widget, page, selfcontained = FALSE)

  session <- chromote::ChromoteSession$new()
  on.exit(try(session$close(), silent = TRUE), add = TRUE)
  session$Page$navigate(paste0("file://", page))
  # Waiting on Page$loadEventFired blocks for ever when the event has already fired, which depends on how
  # heavy the figure is -- small ones raced through and large ones hung. Poll for the figure instead.
  deadline <- Sys.time() + timeout
  repeat {
    ready <- tryCatch(isTRUE(session$Runtime$evaluate(
      "!!(document.querySelector('.js-plotly-plot') && window.Plotly)", returnByValue = TRUE)$result$value),
      error = function(e) FALSE)
    if (ready) break
    if (Sys.time() > deadline) stop("The figure did not appear in the headless browser within ", timeout, "s.", call. = FALSE)
    Sys.sleep(0.25)
  }
  # SVG comes back percent-encoded in the data URI. Decoding that in R means URLdecode over a string of
  # millions of characters, which took minutes and looked like a hang; decodeURIComponent in the browser is
  # native and returns in milliseconds, so the decode happens there and R only writes the result.
  script <- sprintf("(() => { const g = document.querySelector('.js-plotly-plot');
                              if (!g) return '';
                              return Plotly.toImage(g, {format:'%s', width:%d, height:%d, scale:%s})
                                .then(u => %s); })()",
                    format, as.integer(width), as.integer(height), if (format == "svg") 1 else scale,
                    if (format == "svg") "decodeURIComponent(u.slice(u.indexOf(',') + 1))" else "u")
  uri <- session$Runtime$evaluate(script, awaitPromise = TRUE, returnByValue = TRUE, timeout_ = timeout * 1000)$result$value
  if (is.null(uri) || !nzchar(uri)) stop("The figure did not render in the headless browser.", call. = FALSE)

  if (identical(format, "svg")) {
    writeLines(uri, file)                      # already decoded in the browser
  } else {
    writeBin(base64enc::base64decode(sub("^data:image/png;base64,", "", uri)), file)
  }
  invisible(file)
}
