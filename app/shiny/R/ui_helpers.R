project_badge <- function() tags$div(class = "ytab-project-badge", uiOutput("active_project_badge"))
# Switching projects opens the chosen one in place: the menu lists the projects, so leaving the workspace for the
# home page is one item in it rather than the only way through. Bootstrap 3's dropdown JS ships with the
# workspace's navbarPage, so the menu needs no library of its own.
project_switcher <- function() tags$div(class = "ytab-project-switch dropdown",
  tags$button(id = "switch_project_toggle", class = "btn btn-sm dropdown-toggle",
    `data-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
    title = "Choose another dataset without leaving this page",
    onclick = "Shiny.setInputValue('switch_project_opened', Date.now(), {priority: 'event'});",
    "Switch dataset", tags$span(class = "caret")),
  tags$div(class = "dropdown-menu dropdown-menu-right ytab-project-menu",
    uiOutput("switch_project_menu")))
blocked_card <- function(title, prerequisite, action, path = "") tags$div(class = "ytab-blocked",
  tags$strong(title), tags$p(prerequisite), tags$p(tags$b("Next action: "), action), if (nzchar(path)) tags$code(path))
panel_card <- function(title, ...) tags$section(class = "ytab-card", tags$h3(title), ...)
readonly_value <- function(label, output_id) tags$div(class = "ytab-readonly", tags$span(label), textOutput(output_id, inline = TRUE))
# The laboratory's mark, linking out to its site. Defined once so the header and the landing page cannot
# drift apart, and opened in a new tab so a click never loses an analysis in progress.
YTAB_LAB_URL <- "https://www.binhe-lab.org/"
# `caption` names the laboratory under the mark. The workspace header leaves it out: there the mark sits
# beside the project switcher, where a second line of text would compete with the dataset's own name.
ytab_lab_logo <- function(class = "ytab-lab-logo", caption = NULL) {
  tags$a(href = YTAB_LAB_URL, target = "_blank", rel = "noopener noreferrer", class = class,
         title = "He laboratory \u2014 binhe-lab.org",
         tags$img(src = "ytab_lab_logo.png", alt = "He laboratory"),
         if (!is.null(caption)) tags$span(class = "ytab-lab-logo-caption", caption))
}

low_memory_note <- function() tags$div(class = "ytab-note",
  "Mapping is the most memory- and CPU-intensive step. On computers with 4–8 GB RAM, use two threads. Large datasets map faster on an HPC cluster (pipeline/hpc); import their hit files here.")

# A stylesheet URL that changes whenever the app is restarted, and whenever the file itself changes.
#
# The modification time alone is not enough: it is read with a relative path, so if the app is launched from
# anywhere but app/shiny the file cannot be stat-ed, the version is dropped, and the browser then serves its
# cached copy for ever -- an edit to the stylesheet appears to do nothing, which is indistinguishable from a
# rule that does not apply. The session start time is therefore always part of the token, so a restart is
# always enough to see a change.
YTAB_UI_BOOT <- as.integer(Sys.time())

ytab_css_href <- function(name, dir = "www") {
  stamp <- tryCatch(as.integer(file.info(file.path(dir, name))$mtime), error = function(e) NA_integer_)
  if (is.na(stamp)) stamp <- 0L
  sprintf("%s?v=%d-%d", name, YTAB_UI_BOOT, stamp)
}
