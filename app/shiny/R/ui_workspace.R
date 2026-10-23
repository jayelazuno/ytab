project_module_ui <- function() navset_tab(id="project_tabs",
  nav_panel("Samples",value="samples",panel_card("Project and samples",uiOutput("project_metadata_ui"),actionButton("edit_sample_names","Edit sample names"),DT::DTOutput("sample_sheet_table"))),
  nav_panel("Logs & Provenance",value="provenance",panel_card("Logs and provenance",uiOutput("provenance_ui"),uiOutput("last_job_summary"),verbatimTextOutput("pipeline_job_text_logs"))),
  nav_panel("Help",value="help",panel_card("Quick start",uiOutput("help_ui"),low_memory_note())))

reports_exports_modal <- function() modalDialog(title="Reports and exports",size="m",easyClose=TRUE,
  uiOutput("results_readiness"),
  tags$p("Build a summary report of the project, or a compact bundle of its result tables and plots to share."),
  tags$div(class="ytab-actions",actionButton("build_project_report","Build project summary report"),actionButton("build_export_bundle","Create compact export bundle")),
  tags$details(class="ytab-technical-details",tags$summary("Technical details"),verbatimTextOutput("exports_text")),
  footer=modalButton("Close"))

gene_explorer_module_ui <- function() navset_tab(id="gene_explorer_tabs",
  nav_panel("Genome Browser",value="browser",genome_browser_ui()),
  nav_panel("Gene & Domain Insertions",value="gene_domain",gene_domain_explorer_ui()),
  nav_panel("GO Enrichment",value="go",go_analysis_ui()),
  nav_panel("Across Species",value="comparative",comparative_ui()))

workspace_ui <- function(detected_cpu) {
  tags$div(class="ytab-workspace",
    tags$header(class="ytab-workspace-header ytab-release-header",
      tags$div(tags$h2("YTAB"), tags$p(class="ytab-release-subtitle","Yeast Transposon Analysis Browser")),
      project_badge(),
      actionButton("open_reports_exports", "Reports & exports"),
      actionButton("switch_project", "Switch project")),
    job_progress_ui("global_job",compact=TRUE),
    navbarPage(id="workspace_tabs", title=NULL, collapsible=TRUE, inverse=TRUE,
      tabPanel("Project",value="project",project_module_ui()),
      tabPanel("Run",value="run",run_module_ui(detected_cpu)),
      tabPanel("Quality Control",value="quality_control",quality_control_ui()),
      tabPanel("Normalization",value="normalization",normalization_module_ui()),
      tabPanel("Condition Comparisons",value="fitness",comparisons_module_ui()),
      tabPanel("Essentiality",value="essentiality",essentiality_module_ui()),
      tabPanel("Gene Explorer",value="gene_explorer",gene_explorer_module_ui()),
      tabPanel("Screen Planning",value="screen_coverage",screen_coverage_ui()),
      bslib::nav_item(tags$div(class="ytab-nav-tools",
        actionLink("open_search_nav", tagList(icon("magnifying-glass"), "Search"), title = "Search the methods, tutorial and tabs"),
        actionLink("open_help_methods_nav", tagList(icon("book"), "Methods"), title = "Methods and help")))))
}
