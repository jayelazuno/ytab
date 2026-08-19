#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
repo_root <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

must_exist <- function(path) {
  full <- file.path(repo_root, path)
  if (!file.exists(full)) fail(paste("Missing expected file:", path))
  full
}

must_contain <- function(path, pattern, label = pattern) {
  text <- paste(readLines(file.path(repo_root, path), warn = FALSE), collapse = "\n")
  if (!grepl(pattern, text, perl = TRUE)) fail(paste("Missing", label, "in", path))
  invisible(TRUE)
}

must_not_contain <- function(paths, pattern, label = pattern) {
  for (path in paths) {
    if (!file.exists(path) || dir.exists(path)) next
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    if (grepl(pattern, text, perl = TRUE)) fail(paste("Unexpected", label, "in", path))
  }
}

library(shiny)
library(bslib)
`%||%` <- function(left, right) if (is.null(left) || !length(left)) right else left
human_file_size <- function(bytes) {
  bytes <- suppressWarnings(as.numeric(bytes))
  ifelse(is.na(bytes), "", paste0(round(bytes), " B"))
}

invisible(parse(file = must_exist("app/shiny/app.R")))

helper_files <- c(
  "app/shiny/R/ui_components.R",
  "app/shiny/R/plot_display_helpers.R",
  "app/shiny/R/table_display_helpers.R",
  "app/shiny/R/plot_customization_helpers.R"
)
for (helper in helper_files) {
  full <- must_exist(helper)
  invisible(parse(file = full))
  source(full, local = TRUE)
}

invisible(must_exist("app/shiny/www/ytab_release_ui.css"))
must_contain("app/shiny/app.R", "ytab_release_ui\\.css", "release CSS reference")
must_contain("app/shiny/app.R", "ui_components\\.R", "shared UI helper source")
must_contain("app/shiny/app.R", "plot_customization_helpers\\.R", "plot customization helper source")

for (fn in c(
  "ytab_page_header", "ytab_result_card", "ytab_metric_card",
  "ytab_plot_card", "ytab_static_image_card", "ytab_table_card",
  "ytab_technical_details", "ytab_empty_state",
  "ytab_plot_output", "ytab_static_image_ui", "ytab_generated_file_gallery",
  "ytab_static_image_output_card", "ytab_plot_frame", "ytab_two_column_layout",
  "ytab_datatable", "ytab_file_table", "ytab_plot_customization_controls",
  "ytab_plot_customization_values", "ytab_plot_height_px",
  "ytab_with_plot_display_options"
)) {
  if (!exists(fn, mode = "function")) fail(paste("Missing helper function:", fn))
}

must_contain("app/shiny/R/plot_customization_helpers.R", "Display-only controls", "display-only control note")
must_contain("app/shiny/R/plot_customization_helpers.R", "plot_width", "plot width control")
must_contain("app/shiny/R/plot_customization_helpers.R", "label_angle", "label angle control")
must_contain("app/shiny/R/plot_customization_helpers.R", "small = 12L, medium = 15L, large = 19L", "distinct Shiny text-size pixel scale")
must_contain("app/shiny/R/qc_plot_utils.R", "qc_plot_text_sizes", "global QC plot text-size helper")
must_contain("app/shiny/R/qc_plot_utils.R", "small = 0.95, large = 1.55, 1.2", "distinct QC text-size scale")
must_contain("app/shiny/R/qc_plot_utils.R", "axis = 1.05 \\* scale", "QC axis text-size scaling")
must_contain("app/shiny/R/qc_plot_utils.R", "sample = 1.05 \\* scale", "QC sample-label text-size scaling")
must_contain("app/shiny/R/qc_plot_utils.R", "qc_plot_value_cex", "QC value-label text-size helper")
must_contain("app/shiny/R/qc_plot_utils.R", "qc_plot_key_cex", "QC metric-key text-size helper")
must_contain("app/shiny/www/ytab_release_ui.css", "height: auto !important", "aspect-ratio preserving image CSS")
must_contain("app/shiny/www/ytab_release_ui.css", "ytab-two-column-layout", "two-column layout CSS")
must_contain("app/shiny/www/ytab_release_ui.css", "ytab-workspace-header\\.ytab-release-header", "readable header contrast CSS")
must_contain("app/shiny/www/ytab_release_ui.css", "ytab-diagnostic-gallery-grid", "diagnostic gallery density CSS")

workspace_text <- paste(readLines(must_exist("app/shiny/R/ui_workspace.R"), warn = FALSE), collapse = "\n")
for (tab in c("Quality Control", "Essentiality", "Fitness Screen",
              "Gene & Domain Insertion Explorer", "Comparative Species",
              "Results & Exports")) {
  if (!grepl(tab, workspace_text, fixed = TRUE)) fail(paste("Missing major tab:", tab))
}

must_contain("app/shiny/R/ui_qc.R", "summary_qc_plot_choice", "Summary QC plot selector")
must_contain("app/shiny/R/ui_qc.R", "mapping_qc_plot_choice", "Mapping QC plot selector")
must_contain("app/shiny/R/ui_qc.R", "ytab_plot_customization_controls\\(\"mapping_qc\"", "Mapping QC display controls")
must_contain("app/shiny/R/ui_qc.R", "default_bar_orientation = \"horizontal\"", "Mapping QC horizontal default")
must_contain("app/shiny/R/qc_mapping_stats_plot.R", "qc_plot_display_labels", "Mapping QC label mode wiring")
must_contain("app/shiny/R/qc_mapping_stats_plot.R", "qc_plot_grid_enabled", "Mapping QC grid wiring")
must_contain("app/shiny/R/qc_mapping_stats_plot.R", "qc_plot_show_value_labels", "Mapping QC value-label wiring")
must_contain("app/shiny/R/qc_mapping_stats_plot.R", "qc_plot_label_angle", "Mapping QC label-angle wiring")
must_contain("app/shiny/R/qc_mapping_stats_plot.R", "qc_plot_bar_horizontal", "Mapping QC orientation wiring")
must_contain("app/shiny/R/qc_mapping_stats_plot.R", "legend\\(\"center\"", "Mapping QC metric key row")
must_contain("app/shiny/R/qc_mapping_stats_plot.R", "layout\\(matrix", "Mapping QC separate metric key row")
must_contain("app/shiny/R/qc_mapping_stats_plot.R", "% reads mapped", "Mapping QC readable mapped legend label")
must_contain("app/shiny/R/qc_mapping_stats_plot.R", "% HQ aligned", "Mapping QC readable HQ legend label")
must_contain("app/shiny/R/qc_mapping_stats_plot.R", "qc_plot_label_margin_lines", "Mapping QC dynamic label margins")
must_not_contain(file.path(repo_root, "app/shiny/R/qc_mapping_stats_plot.R"),
                 "% l",
                 "truncated Mapping QC legend label")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "App-rendered: Mapping summary",
                 "dominant app-rendered Mapping QC title")
must_not_contain(file.path(repo_root, "app/shiny/R/qc_mapping_stats_plot.R"),
                 "abbreviate\\(|strtrim\\(",
                 "Mapping QC silent label truncation")
must_contain("app/shiny/R/ui_qc.R", "ytab_plot_customization_controls\\(\"summary_qc\"", "Summary QC display controls")
must_contain("app/shiny/R/ui_qc.R", "default_bar_orientation = \"horizontal\"", "QC horizontal long-label default")
must_contain("app/shiny/R/ui_qc.R", "ytab_two_column_layout", "Summary QC two-column layout")
must_contain("app/shiny/R/qc_summary_library_plots.R", "qc_plot_label_margin_lines", "Summary QC dynamic label margins")
must_contain("app/shiny/R/qc_summary_library_plots.R", "qc_plot_draw_vertical_labels", "Summary QC vertical label drawing")
must_contain("app/shiny/R/qc_summary_library_plots.R", "qc_plot_draw_horizontal_labels", "Summary QC horizontal label drawing")
must_contain("app/shiny/R/qc_summary_library_plots.R", "qc_plot_metric_key_row", "Summary QC metric key row")
must_contain("app/shiny/R/qc_summary_library_plots.R", "qc_plot_show_value_labels", "Summary QC value-label control wiring")
must_contain("app/shiny/R/ui_qc.R", "summary_qc_combined_group_selector", "Summary QC combined group selector UI")
must_contain("app/shiny/app.R", "summary_qc_combined_group_selector", "Summary QC combined group selector server")
must_contain("app/shiny/app.R", "input\\$summary_qc_combined_group", "Summary QC combined group plot input")
must_contain("app/shiny/R/qc_summary_library_plots.R", "qc_summary_combined_feature_group_choices", "Summary QC combined group choices")
must_contain("app/shiny/R/qc_summary_library_plots.R", "qc_summary_combined_feature_sample_roles", "Summary QC metadata-aware combined group detection")
must_contain("app/shiny/R/qc_summary_library_plots.R", "control_or_treated", "Summary QC control/treated metadata priority")
must_contain("app/shiny/R/qc_summary_library_plots.R", "bar_half <- if \\(n == 1L\\) 0\\.22", "Summary QC single-bar width cap")
must_contain("app/shiny/R/qc_summary_library_plots.R", "lwd = max\\(1\\.6, qc_plot_lwd\\)", "Summary QC visible combined bar border")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "qc_plot_metric_key_row", "Library Diagnostics metric/sample key rows")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "qc_plot_begin_key_layout", "Library Diagnostics separate key layout")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "qc_plot_bar_horizontal", "Library Diagnostics bar orientation wiring")
must_contain("app/shiny/R/ui_qc.R", "\"MidLC saturation\" = \"midlc\"", "Library Diagnostics clean MidLC label")
must_contain("app/shiny/R/ui_qc.R", "\"Centromere bias\" = \"centromere\"", "Library Diagnostics clean centromere label")
must_contain("app/shiny/R/ui_qc.R", "\"Feature metaplots\" = \"metaplots\"", "Library Diagnostics clean metaplot label")
must_contain("app/shiny/R/ui_qc.R", "library_diagnostics_color_by", "Library Diagnostics color-by control")
must_contain("app/shiny/app.R", "library_diagnostics_group_selector", "Library Diagnostics sample-group selector")
must_contain("app/shiny/app.R", "library_diagnostics_metaplot_panel", "Library Diagnostics metaplot panel selector")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "plot_qc_library_centromere_bias", "Library Diagnostics live centromere plot")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "plot_qc_library_metaplot", "Library Diagnostics live feature metaplots")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "qc_library_sample_roles", "Library Diagnostics metadata-aware sample grouping")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "qc_library_group_color_map", "Library Diagnostics stable group color map")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "qc_library_group_color", "Library Diagnostics stable group color helper")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "qc_library_plot_aesthetics", "Library Diagnostics plot-specific aesthetic helper")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "midlc = modifyList\\(defaults, list\\(use_points = TRUE", "MidLC line-plus-point geometry")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "type = if \\(isTRUE\\(aesthetic\\$use_points\\)\\) \"b\" else \"l\"", "MidLC point/line rendering")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "pair_connectors = FALSE", "Library Diagnostics no global pair connectors")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "sample_shapes = FALSE", "Library Diagnostics no global sample-shape mapping")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "aggregate\\(mean_reads ~ display_group \\+ role \\+ distance_kb", "Centromere bias group summary curves")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "aggregate\\(mean_reads_per_site ~ display_group \\+ role \\+ rel_bp", "Feature metaplot group summary curves")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "aggregate\\(enrichment ~ role \\+ display_group \\+ pair", "Sequence bias grouped by biological group")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "out <- unique\\(out\\)", "Library Diagnostics duplicate filtered-row guard")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "droplevels\\(out\\)", "Library Diagnostics filtered factor cleanup")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "qc_library_plot_cache_key", "Library Diagnostics plot cache/reactive key helper")
must_contain("app/shiny/app.R", "qc_library_plot_cache_key\\(active\\(\\),\"sequence_bias\"", "Library Diagnostics sequence-bias reactive key")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "App-rendered:",
                 "Library Diagnostics prominent app-rendered label")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "Generated:",
                 "Library Diagnostics prominent generated label")
must_not_contain(file.path(repo_root, "app/shiny/R/qc_library_diagnostics_plots.R"),
                 "Static generated image",
                 "Library Diagnostics prominent static-image badge")
must_contain("app/shiny/R/ui_qc.R", "download_mapping_qc_plot", "Mapping QC download plot control")
must_contain("app/shiny/R/ui_qc.R", "download_mapping_qc_plotted_data", "Mapping QC download plotted-data control")
must_contain("app/shiny/R/ui_qc.R", "download_summary_qc_plot", "Summary QC download plot control")
must_contain("app/shiny/R/ui_qc.R", "download_summary_qc_plotted_data", "Summary QC download plotted-data control")
must_contain("app/shiny/R/ui_qc.R", "download_library_diagnostics_plot", "Library Diagnostics download plot control")
must_contain("app/shiny/R/ui_qc.R", "download_library_diagnostics_plotted_data", "Library Diagnostics download plotted-data control")
must_contain("app/shiny/R/ui_qc.R", "library_diagnostics_plot_choice", "Library Diagnostics visualization selector")
must_contain("app/shiny/R/ui_qc.R", "ytab_plot_customization_controls\\(\\\"library_diagnostics\\\"", "Library Diagnostics plot display options")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "job_progress_ui\\(\\\"library_diagnostics_job\\\"\\)",
                 "Library Diagnostics job-progress card in normal plot controls")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "uiOutput\\(\\\"library_diagnostics_result\\\"\\)",
                 "Library Diagnostics result/history cards in normal plot controls")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "No matching diagnostics result",
                 "Library Diagnostics no-match card in normal UI")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "Latest historical result",
                 "Library Diagnostics historical-result card in normal UI")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "No job progress is available",
                 "Library Diagnostics no-progress card in normal UI")
must_contain("app/shiny/R/ui_qc.R", "diagnostic_files_table", "Diagnostic Files table remains available")
must_contain("app/shiny/R/ui_qc.R", "diagnostic_show_archived_static", "Diagnostic Files archived/static image opt-in")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "Plot gallery",
                 "Diagnostic Files normal plot gallery option")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "diagnostic_plot_cards",
                 "Diagnostic Files normal static image cards")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "Live plot",
                 "prominent live-plot label in QC UI")
must_contain("app/shiny/app.R", "download_mapping_qc_plot", "Mapping QC plot download handler")
must_contain("app/shiny/app.R", "download_summary_qc_plot", "Summary QC plot download handler")
must_contain("app/shiny/app.R", "download_library_diagnostics_plot", "Library Diagnostics plot download handler")
must_contain("app/shiny/app.R", "!isTRUE\\(input\\$diagnostic_show_archived_static\\)", "Diagnostic Files hides archived static images by default")
must_contain("app/shiny/app.R", "Archived static diagnostic PNGs are not shown as current outputs", "stale static image fallback message")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_qc.R"),
                 "App-rendered: Summary QC",
                 "dominant app-rendered Summary QC title")
must_contain("app/shiny/R/ui_essentiality.R", "ytab_two_column_layout", "Essentiality two-column layout")
must_contain("app/shiny/R/ui_essentiality.R", "ytab_plot_customization_controls\\(\"essentiality\"", "Essentiality display controls")
must_contain("app/shiny/R/ui_fitness.R", "fitness_visualization_selector", "Fitness plot selector")
must_contain("app/shiny/R/ui_fitness.R", "ytab_plot_customization_controls\\(\"fitness\"", "Fitness display controls")
must_contain("app/shiny/app.R", "\"MA plot\"=\"combined_ma\"", "dynamic Fitness MA plot option")
must_contain("app/shiny/app.R", "Condition versus control log-log scatter", "Fitness condition-control plot option")
must_contain("app/shiny/app.R", "Top selected hits log2FC heatmap", "Fitness selected-hit heatmap option")
must_contain("app/shiny/app.R", "fitness_condition_control_plot", "Fitness condition-control live plot output")
must_contain("app/shiny/app.R", "fitness_selected_hit_heatmap", "Fitness selected-hit heatmap live plot output")
must_contain("app/shiny/app.R", "plot_fitness_condition_control_scatter", "Fitness condition-control live plot renderer")
must_contain("app/shiny/app.R", "plot_fitness_selected_hit_heatmap", "Fitness selected-hit heatmap renderer")
must_contain("app/shiny/app.R", "download_fitness_heatmap_data", "Fitness selected-hit heatmap data download")
must_contain("app/shiny/app.R", "fitness_ma_mode", "Fitness MA comparison mode")
must_contain("app/shiny/app.R", "fitness_ma_hit_direction", "Fitness MA hit direction")
must_contain("app/shiny/app.R", "fitness_ma_top_n", "Fitness MA top-hit count")
must_contain("app/shiny/app.R", "fitness_ma_annotation_mode", "Fitness MA annotation mode")
must_contain("app/shiny/app.R", "fitness_ma_min_support", "Fitness MA minimum supporting pools control")
must_contain("app/shiny/app.R", "Concordance filter", "Fitness MA support-control label")
must_contain("app/shiny/app.R", "All support classes", "Fitness MA no-support-filter default")
must_contain("app/shiny/app.R", "fitness_ma_selection_status", "Fitness MA candidate/annotation status")
must_contain("app/shiny/app.R", "Ranked candidates passing log2FC/CPM filters", "directional evidence status wording")
must_contain("app/shiny/app.R", "download_fitness_ma_plot", "Fitness MA plot download")
must_contain("app/shiny/app.R", "download_fitness_ma_data", "Fitness MA plotted-data download")
must_contain("app/shiny/app.R", "download_fitness_ma_top_hits", "Fitness MA selected top-hits download")
must_contain("app/shiny/app.R", "fitness_ma_selected_hits", "shared selected top-hits reactive")
must_contain("app/shiny/app.R", "fitness_ma_selected_hits_visible", "selected top-hits search-filtered table view")
must_contain("app/shiny/R/ui_fitness.R", "fitness_selected_top_hits_table", "selected top-hits table UI")
must_contain("app/shiny/R/ui_fitness.R", "Search gene or feature ID", "selected top-hits search control")
must_contain("app/shiny/R/ui_fitness.R", "Full result fitness-call filter", "full-result-only fitness call filter")
must_contain("app/shiny/R/ui_fitness.R", "abs\\(log2FC\\) >= 1 and CPM/read support >= 1", "selected top-hits threshold note")
must_contain("app/shiny/R/fitness_plot_utils.R", "passes_log2FC_threshold", "selected-hit log2FC eligibility")
must_contain("app/shiny/R/fitness_plot_utils.R", "passes_CPM_threshold", "selected-hit CPM eligibility")
must_contain("app/shiny/app.R", "Supporting pool IDs", "selected top-hits supporting pool column")
must_contain("app/shiny/app.R", "order = list\\(list\\(0, \"asc\"\\)\\)", "selected top-hits Rank ascending order")
must_contain("app/shiny/app.R", "scrollX = TRUE", "selected top-hits horizontal scrolling")
must_contain("app/shiny/app.R", "Showing %d selected top hit", "selected top-hits count wording")
must_contain("app/shiny/app.R", "Showing %d of %d selected top hits matching search", "selected top-hits search count wording")
must_contain("app/shiny/app.R", "names\\(x\\)\\[names\\(x\\) == \"selected_rank\"\\] <- \"Rank\"", "selected top-hits visible Rank column")
must_contain("app/shiny/R/fitness_plot_utils.R", "overall_candidate_rank", "selected top-hits CSV audit rank")
must_contain("app/shiny/R/fitness_plot_utils.R", "Depleted/enriched colors reflect stored calls from local-abundance z-scores.", "Fitness MA call explanation")
must_contain("app/shiny/R/fitness_plot_utils.R", "present <- c\\(depleted", "Fitness MA dynamic call legend")
must_contain("app/shiny/R/fitness_plot_utils.R", "selected_hit_colors", "Fitness MA selected-hit direction colors")
must_contain("app/shiny/R/fitness_plot_utils.R", "Depleted selected hit", "Fitness MA selected depleted legend")
must_contain("app/shiny/R/fitness_plot_utils.R", "Enriched selected hit", "Fitness MA selected enriched legend")
must_contain("app/shiny/R/fitness_plot_utils.R", "M: mean abundance, log10\\(CPM \\+ 1\\)", "Fitness MA x-axis abundance label")
must_contain("app/shiny/R/fitness_plot_utils.R", "A: fitness effect, log2FC treated/control", "Fitness MA y-axis effect label")
must_contain("app/shiny/R/fitness_plot_utils.R", "More than 10 labels may overlap; reduce Number of hits for cleaner labels.", "Fitness MA label crowding note")
must_contain("app/shiny/R/fitness_plot_utils.R", "fitness_ma_rank_data", "Fitness MA evidence ranking helper")
must_contain("app/shiny/R/fitness_plot_utils.R", "fitness_ma_top_hits_table", "Fitness MA selected top-hits table helper")
must_contain("app/shiny/R/fitness_plot_utils.R", "fitness_selected_hit_heatmap_data", "Fitness selected-hit heatmap data helper")
must_contain("app/shiny/R/fitness_plot_utils.R", "Rows = currently selected hits", "Fitness selected-hit heatmap note")
must_contain("app/shiny/R/fitness_plot_utils.R", "log2FC treated/control", "Fitness selected-hit heatmap value label")
must_contain("app/shiny/R/fitness_plot_utils.R", "ranked by directional log2FC magnitude, CPM/read support, local z-score support, and stable feature ID", "correct Fitness MA ranking priority")
must_contain("app/shiny/R/fitness_condition_control_plot.R", "parent_cpm", "condition-control control abundance column")
must_contain("app/shiny/R/fitness_condition_control_plot.R", "treated_cpm", "condition-control treated abundance column")
must_contain("app/shiny/R/fitness_condition_control_plot.R", "Control abundance, log10\\(CPM \\+ 1\\)", "condition-control x-axis label")
must_contain("app/shiny/R/fitness_condition_control_plot.R", "Treated abundance, log10\\(CPM \\+ 1\\)", "condition-control y-axis label")
must_contain("app/shiny/R/fitness_condition_control_plot.R", "abline\\(0, 1", "condition-control 1:1 reference line")
must_contain("app/shiny/R/fitness_condition_control_plot.R", "fitness_ma_highlight_rows", "condition-control shared selected-hit rows")
must_contain("app/shiny/R/fitness_generated_plots.R", "top_features_log2FC_heatmap.png", "legacy static Fitness heatmap excluded from generated inventory")
for (column in c("rank_z_strength", "rank_lfc_strength", "rank_cpm_support", "depleted_support_n", "enriched_support_n", "candidate_rank_order"))
  must_contain("app/shiny/R/fitness_plot_utils.R", column, paste("Fitness MA ranking column", column))
must_contain("app/shiny/R/fitness_plot_utils.R", "fitness_ma_available_pools", "dynamic pool-support choices helper")
must_not_contain(file.path(repo_root, "app/shiny/R/fitness_plot_utils.R"), "top100_depleted_in_treated", "top100 files as dynamic MA ranking source")
must_contain("app/shiny/R/fitness_plot_utils.R", "label_table\\$data", "Fitness MA current annotation layer")
must_not_contain(file.path(repo_root, "app/shiny/R/ui_fitness.R"), "Generated: Combined treated", "static MA label in Fitness UI")
must_contain("app/shiny/R/ui_fitness.R", "ytab_two_column_layout", "Fitness two-column layout")
must_contain("app/shiny/app.R", "ytab-static-image-card", "Fitness generated plot static card")
must_contain("app/shiny/app.R", "ytab_plot_frame", "App-rendered plot frame")
must_contain("app/shiny/R/ui_gene_domain_explorer.R", "gene_domain_track_preset", "Gene Explorer track preset")
must_contain("app/shiny/R/ui_gene_domain_explorer.R", "gene_domain_display_height", "Gene Explorer display height")
must_contain("app/shiny/R/gene_domain_explorer_server.R", "ytab_static_image_output_card", "Gene Explorer static image card")
must_contain("app/shiny/R/gene_domain_explorer_state.R", "gene_domain_preset_track_rows", "Gene Explorer preset logic")
must_contain("app/shiny/R/gene_domain_explorer_server.R", "Gene Explorer could not generate a figure", "Gene Explorer friendly failure state")
must_contain("app/shiny/R/gene_domain_explorer_server.R", "Search for a gene, choose tracks, then click Generate figure", "Gene Explorer friendly empty state")
must_contain("app/shiny/R/ui_comparative.R", "comparative_single_species", "Comparative species controls")
must_contain("app/shiny/R/ui_comparative.R", "comparative_ortholog_table", "Comparative orthology table")
must_contain("app/shiny/R/plot_display_helpers.R", "ytab_generated_file_gallery", "generated file gallery helper")
must_contain("app/shiny/R/fitness_generated_plots.R", "ytab-static-image-card", "Fitness generated plot static card")
must_contain("app/shiny/R/qc_library_diagnostics_plots.R", "ytab-static-image-card", "Diagnostics generated plot static card helper retained")
must_contain("app/shiny/app.R", "View · Show path", "Diagnostic Files file table action text")
must_contain("app/shiny/R/qc_summary_library_plots.R", "qc_plot_bar_horizontal", "Summary QC long-label handling")
must_contain("docs/ui_reference_research_step14F1.md", "Runtime dependency confirmation", "research manifest")

runtime_files <- file.path(repo_root, c(
  "app/shiny/app.R",
  "app/shiny/R/ui_workspace.R",
  "app/shiny/R/ui_qc.R",
  "app/shiny/R/ui_fitness.R",
  "app/shiny/R/ui_gene_domain_explorer.R",
  "app/shiny/R/gene_domain_explorer_server.R",
  "app/shiny/R/ui_comparative.R",
  helper_files
))
runtime_files <- runtime_files[file.exists(runtime_files)]
reference_name <- paste0("RNA", "cross")
external_reference_dir <- paste0("codex/", "ui_", "references")
blocked <- paste0("codex/", reference_name, "|source\\(.*", reference_name,
                  "|readRDS\\(.*codex/", reference_name, "|", external_reference_dir)
must_not_contain(runtime_files, blocked, "runtime dependency on external UI reference")

cat("PASS\n")
