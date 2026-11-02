project_badge <- function() tags$div(class = "ytab-project-badge", uiOutput("active_project_badge"))
# Switching projects opens the chosen one in place: the menu lists the projects, so leaving the workspace for the
# home page is one item in it rather than the only way through. Bootstrap 3's dropdown JS ships with the
# workspace's navbarPage, so the menu needs no library of its own.
project_switcher <- function() tags$div(class = "ytab-project-switch dropdown",
  tags$button(id = "switch_project_toggle", class = "btn btn-sm dropdown-toggle",
    `data-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
    title = "Choose another project without leaving this page",
    onclick = "Shiny.setInputValue('switch_project_opened', Date.now(), {priority: 'event'});",
    "Switch project", tags$span(class = "caret")),
  tags$div(class = "dropdown-menu dropdown-menu-right ytab-project-menu",
    uiOutput("switch_project_menu")))
blocked_card <- function(title, prerequisite, action, path = "") tags$div(class = "ytab-blocked",
  tags$strong(title), tags$p(prerequisite), tags$p(tags$b("Next action: "), action), if (nzchar(path)) tags$code(path))
panel_card <- function(title, ...) tags$section(class = "ytab-card", tags$h3(title), ...)
readonly_value <- function(label, output_id) tags$div(class = "ytab-readonly", tags$span(label), textOutput(output_id, inline = TRUE))
low_memory_note <- function() tags$div(class = "ytab-note",
  "Mapping is the most memory- and CPU-intensive step. On computers with 4–8 GB RAM, use two threads. Large projects map faster on an HPC cluster (pipeline/hpc); import their hit files here.")
