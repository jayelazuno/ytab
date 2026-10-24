project_badge <- function() tags$div(class = "ytab-project-badge", uiOutput("active_project_badge"))
blocked_card <- function(title, prerequisite, action, path = "") tags$div(class = "ytab-blocked",
  tags$strong(title), tags$p(prerequisite), tags$p(tags$b("Next action: "), action), if (nzchar(path)) tags$code(path))
panel_card <- function(title, ...) tags$section(class = "ytab-card", tags$h3(title), ...)
readonly_value <- function(label, output_id) tags$div(class = "ytab-readonly", tags$span(label), textOutput(output_id, inline = TRUE))
low_memory_note <- function() tags$div(class = "ytab-note",
  "Mapping is the most memory- and CPU-intensive step. On computers with 4–8 GB RAM, use two threads. Large projects map faster on an HPC cluster (pipeline/hpc); import their hit files here.")
