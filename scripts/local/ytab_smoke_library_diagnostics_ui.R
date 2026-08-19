#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: ytab_smoke_library_diagnostics_ui.R <project.yaml>", call. = FALSE)
all_args <- commandArgs(trailingOnly = FALSE)
self <- normalizePath(sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[[1L]]),
                      winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(self), "../.."), winslash = "/", mustWork = TRUE)

library(shiny)
library(bslib)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left

source(file.path(repo_root, "app/shiny/R/project_discovery.R"))
source(file.path(repo_root, "app/shiny/R/ui_helpers.R"))
source(file.path(repo_root, "app/shiny/R/ui_components.R"))
source(file.path(repo_root, "app/shiny/R/plot_customization_helpers.R"))
source(file.path(repo_root, "app/shiny/R/plot_display_helpers.R"))
source(file.path(repo_root, "app/shiny/R/table_display_helpers.R"))
source(file.path(repo_root, "app/shiny/R/job_progress.R"))
source(file.path(repo_root, "app/shiny/R/qc_result_state.R"))
source(file.path(repo_root, "app/shiny/R/qc_plot_utils.R"))
source(file.path(repo_root, "app/shiny/R/qc_library_diagnostics_plots.R"))
source(file.path(repo_root, "app/shiny/R/ui_qc.R"))

project <- read_project_summary(args[[1L]], repo_root)
project_id <- as.character(project$project_id %||% "")

render_plot <- function(expr, width = 1600, height = 1000) {
  input <- list(
    library_diagnostics_text_size = "large",
    library_diagnostics_label_mode = "full",
    library_diagnostics_label_angle = "0",
    library_diagnostics_grid = "show",
    library_diagnostics_bar_orientation = "horizontal",
    library_diagnostics_show_value_labels = FALSE
  )
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = width, height = height, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  ytab_with_plot_display_options(input, "library_diagnostics", expr)
  grDevices::dev.off()
  on.exit(NULL)
  stopifnot(file.exists(file), file.info(file)$size > 1000)
  invisible(file)
}

midlc <- qc_library_midlc_plot_data(project)
summary_tbl <- qc_library_summary_plot_data(project)
centromere <- qc_library_centromere_plot_data(project)
metaplot <- qc_library_metaplot_plot_data(project)
seqbias <- qc_library_sequence_bias_plot_data(project)

stopifnot(nrow(midlc) > 0L, nrow(summary_tbl) > 0L)
choices <- qc_library_group_choices(project)
choice_values <- unname(unlist(choices, use.names = FALSE))
stopifnot("all" %in% choice_values)
if (project_id %in% c("Zn_toxicity_screen", "H2O2_screen_v1")) {
  stopifnot(all(c("control", "treated", "both") %in% choice_values))
}
if (identical(project_id, "Zn_toxicity_screen")) {
  stopifnot(any(grepl("Mock control", names(choices))),
            any(grepl("1.5 mM Zn-treated", names(choices), fixed = TRUE)))
}
if (identical(project_id, "H2O2_screen_v1")) {
  stopifnot(any(grepl("Parent/control|Parents", names(choices))),
            any(grepl("H2O2-treated", names(choices), fixed = TRUE)))
}

color_choices <- qc_library_color_choices()
stopifnot(all(c("group", "pool", "sample", "none") %in% unname(color_choices)))
aesthetic_cent <- qc_library_plot_aesthetics("centromere")
aesthetic_midlc <- qc_library_plot_aesthetics("midlc")
aesthetic_seq <- qc_library_plot_aesthetics("sequence_bias")
stopifnot(
  identical(aesthetic_cent$use_points, FALSE),
  identical(aesthetic_cent$sample_shapes, FALSE),
  identical(aesthetic_cent$sample_linetype, FALSE),
  identical(aesthetic_cent$pair_connectors, FALSE),
  identical(aesthetic_cent$summary_line, TRUE),
  identical(aesthetic_midlc$line_by_sample, TRUE),
  identical(aesthetic_midlc$use_points, TRUE),
  identical(aesthetic_midlc$summary_line, FALSE),
  identical(aesthetic_seq$use_bars, TRUE),
  identical(aesthetic_seq$use_lines, FALSE),
  identical(aesthetic_seq$use_points, FALSE)
)
roles <- qc_library_sample_roles(project, unique(as.character(midlc$sample)))
stopifnot(all(c("sample", "role", "display_group", "pool") %in% names(roles)))
control_samples <- roles$sample[roles$role == "control"]
treated_samples <- roles$sample[roles$role == "treated"]
if (length(control_samples)) {
  control_midlc <- qc_library_filter_data(project, midlc, "control")
  stopifnot(nrow(control_midlc) > 0L,
            all(as.character(control_midlc$sample) %in% control_samples),
            !any(as.character(control_midlc$sample) %in% treated_samples))
}
if (length(treated_samples)) {
  treated_midlc <- qc_library_filter_data(project, midlc, "treated")
  stopifnot(nrow(treated_midlc) > 0L,
            all(as.character(treated_midlc$sample) %in% treated_samples),
            !any(as.character(treated_midlc$sample) %in% control_samples))
}
both_midlc <- qc_library_filter_data(project, midlc, "both")
stopifnot(nrow(both_midlc) <= nrow(midlc),
          all(as.character(both_midlc$role) %in% c("control", "treated")))
if (identical(project_id, "Zn_toxicity_screen")) {
  stopifnot(
    length(unique(as.character(qc_library_filter_data(project, midlc, "all")$sample))) == 8L,
    length(unique(as.character(qc_library_filter_data(project, midlc, "control")$sample))) == 4L,
    length(unique(as.character(qc_library_filter_data(project, midlc, "treated")$sample))) == 4L
  )
  if (nrow(centromere)) {
    stopifnot(
      length(unique(as.character(qc_library_filter_data(project, centromere, "all")$sample))) == 8L,
      length(unique(as.character(qc_library_filter_data(project, centromere, "control")$sample))) == 4L,
      length(unique(as.character(qc_library_filter_data(project, centromere, "treated")$sample))) == 4L
    )
  }
  if (nrow(metaplot)) {
    stopifnot(
      length(unique(as.character(qc_library_filter_data(project, metaplot, "all")$sample))) == 8L,
      length(unique(as.character(qc_library_filter_data(project, metaplot, "control")$sample))) == 4L,
      length(unique(as.character(qc_library_filter_data(project, metaplot, "treated")$sample))) == 4L
    )
  }
}

colors <- qc_library_group_color_map()
stopifnot(identical(unname(colors[["control"]]), "#2c7fb8"),
          identical(unname(colors[["treated"]]), "#d95f02"),
          identical(qc_library_group_color("treated"), "#d95f02"),
          identical(qc_library_group_color("control"), "#2c7fb8"))
if (nrow(seqbias) && length(treated_samples)) {
  treated_seq <- qc_library_filter_data(project, seqbias, "treated")
  stopifnot(nrow(treated_seq) > 0L,
            all(as.character(treated_seq$role) == "treated"),
            identical(unique(qc_library_group_color(treated_seq$role)), "#d95f02"))
}
if (nrow(seqbias) && length(control_samples)) {
  control_seq <- qc_library_filter_data(project, seqbias, "control")
  stopifnot(nrow(control_seq) > 0L,
            all(as.character(control_seq$role) == "control"),
            identical(unique(qc_library_group_color(control_seq$role)), "#2c7fb8"))
}
sample_style <- qc_library_sample_style(qc_library_filter_data(project, midlc, "treated"), "group")
if (length(treated_samples)) stopifnot(all(sample_style$key$color == "#d95f02"))
sample_style <- qc_library_sample_style(qc_library_filter_data(project, midlc, "control"), "group")
if (length(control_samples)) stopifnot(all(sample_style$key$color == "#2c7fb8"))

render_plot(plot_qc_library_midlc(project, "both", "group"))
render_plot(plot_qc_library_jackpot_depth(project, "both", "group"))
if (nrow(seqbias)) render_plot(plot_qc_library_sequence_bias(project, "both"))
if (nrow(seqbias) && length(treated_samples)) render_plot(plot_qc_library_sequence_bias(project, "treated"))
if (nrow(centromere)) render_plot(plot_qc_library_centromere_bias(project, "both", "group"))
if (nrow(metaplot)) {
  panels <- unique(as.character(metaplot$feature))
  stopifnot(all(c("tss", "tts", "trna") %in% panels))
  render_plot(plot_qc_library_metaplot(project, "tss", "both", "group"))
  render_plot(plot_qc_library_metaplot(project, "tts", "both", "group"))
  render_plot(plot_qc_library_metaplot(project, "trna", "both", "group"))
}

ui_text <- paste(readLines(file.path(repo_root, "app/shiny/R/ui_qc.R"), warn = FALSE), collapse = "\n")
app_text <- paste(readLines(file.path(repo_root, "app/shiny/app.R"), warn = FALSE), collapse = "\n")
plot_text <- paste(readLines(file.path(repo_root, "app/shiny/R/qc_library_diagnostics_plots.R"), warn = FALSE), collapse = "\n")

stopifnot(
  grepl('"MidLC saturation" = "midlc"', ui_text, fixed = TRUE),
  grepl('"Centromere bias" = "centromere"', ui_text, fixed = TRUE),
  grepl('"Feature metaplots" = "metaplots"', ui_text, fixed = TRUE),
  !grepl("App-rendered:", ui_text, fixed = TRUE),
  !grepl("Generated:", ui_text, fixed = TRUE),
  !grepl("Live plot", ui_text, fixed = TRUE),
  !grepl("Static generated image", plot_text, fixed = TRUE),
  grepl("download_library_diagnostics_plot", ui_text, fixed = TRUE),
  grepl("download_library_diagnostics_plotted_data", ui_text, fixed = TRUE),
  grepl("download_library_diagnostics_plot", app_text, fixed = TRUE),
  grepl("download_library_diagnostics_plotted_data", app_text, fixed = TRUE),
  grepl('selectInput("library_diagnostics_plot_choice", "Visualization"', ui_text, fixed = TRUE),
  grepl('library_diagnostics_group_selector', ui_text, fixed = TRUE),
  grepl('ytab_plot_customization_controls("library_diagnostics"', ui_text, fixed = TRUE),
  !grepl('job_progress_ui("library_diagnostics_job")', ui_text, fixed = TRUE),
  !grepl('uiOutput("library_diagnostics_result")', ui_text, fixed = TRUE),
  !grepl('No matching diagnostics result', ui_text, fixed = TRUE),
  !grepl('Latest historical result', ui_text, fixed = TRUE),
  !grepl('No job progress is available', ui_text, fixed = TRUE),
  grepl("library_diagnostics_group_selector", app_text, fixed = TRUE),
  grepl("library_diagnostics_color_by", ui_text, fixed = TRUE),
  grepl("library_diagnostics_metaplot_panel", app_text, fixed = TRUE),
  grepl("qc_library_plot_cache_key", app_text, fixed = TRUE),
  grepl("sample_group", plot_text, fixed = TRUE),
  grepl("qc_library_group_color_map", plot_text, fixed = TRUE),
  grepl("qc_library_group_color", plot_text, fixed = TRUE),
  grepl("qc_library_plot_aesthetics", plot_text, fixed = TRUE),
  grepl("midlc = modifyList\\(defaults, list\\(use_points = TRUE", plot_text),
  grepl("type = if \\(isTRUE\\(aesthetic\\$use_points\\)\\) \"b\" else \"l\"", plot_text),
  grepl("pair_connectors = FALSE", plot_text, fixed = TRUE),
  grepl("aggregate\\(mean_reads ~ display_group \\+ role \\+ distance_kb", plot_text),
  grepl("aggregate\\(mean_reads_per_site ~ display_group \\+ role \\+ rel_bp", plot_text),
  !grepl("adjustcolor", plot_text, fixed = TRUE),
  grepl("aggregate\\(enrichment ~ role \\+ display_group \\+ pair", plot_text),
  grepl("Source: existing centromere_bins.tsv tables", app_text, fixed = TRUE),
  grepl("Source: existing TSS/TTS/tRNA metaplot TSV tables", app_text, fixed = TRUE),
  grepl("Archived static diagnostic PNGs are not shown as current outputs", app_text, fixed = TRUE),
  grepl("src/ytab/qc/LibraryDiagnostics.py", paste(capture.output(tools::md5sum(file.path(repo_root, "src/ytab/qc/LibraryDiagnostics.py"))), collapse = "\n")) ||
    file.exists(file.path(repo_root, "src/ytab/qc/LibraryDiagnostics.py")),
  file.exists(file.path(repo_root, "scripts/plot_gene_feature_metaplots.R")),
  file.exists(file.path(repo_root, "scripts/plot_centromere_bias.R")),
  file.exists(file.path(repo_root, "docs/codex/LibraryQC.R"))
)

cat("PASS\n")
