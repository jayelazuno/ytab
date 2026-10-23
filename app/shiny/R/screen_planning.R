# Screen Planning module: the analytical Poisson / equal-gene coverage planner (from the project's parent libraries)
# and mutant identification. Replaces screen_coverage.R + screen_coverage_active.R; the retired Monte Carlo pool
# planner UI and server that screen_coverage_active.R overrode are kept in archive/ui_reorganization/backup/.

sc_theme <- function(values) {
  size <- values$text_size + switch(values$style,compact=-2,presentation=3,0)
  ggplot2::theme_bw(base_size=size) + ggplot2::theme(
    text=ggplot2::element_text(colour="black"),axis.text=ggplot2::element_text(face="bold",colour="black"),
    axis.title=ggplot2::element_text(face="bold"),plot.title=ggplot2::element_text(face="bold"),
    panel.grid.minor=ggplot2::element_blank(),
    panel.grid.major=if(identical(values$grid,"show"))ggplot2::element_line(colour="grey90")else ggplot2::element_blank(),
    legend.position="bottom",legend.title=ggplot2::element_blank(),plot.margin=ggplot2::margin(12,20,12,12))
}

# Screen Coverage needs scipy and PyYAML in addition to the lighter-weight
# dependencies used by some other YTAB features. Prefer the configured app
# interpreter, but fall back to the repository's declared ytab-local conda
# environment when the app was launched from a system Python/R session.
sc_python_command <- function(python_bin, repo_root) {
  probe <- function(command, prefix = character()) {
    if (!nzchar(command) || !file.exists(command)) return(FALSE)
    check <- tryCatch(processx::run(command,
      c(prefix, "-c", "import numpy, scipy, yaml"), error_on_status = FALSE),
      error = function(e) NULL)
    !is.null(check) && identical(check$status, 0L)
  }
  configured_python <- tryCatch(python_bin(), error = function(e) "")
  candidates <- unique(c(Sys.getenv("YTAB_PYTHON", unset = ""), configured_python,
    file.path(repo_root, ".venv", "bin", "python")))
  for (candidate in candidates) {
    if (probe(candidate)) return(list(command = candidate, prefix = character()))
  }
  conda_candidates <- unique(c(Sys.which("conda"),
    file.path(dirname(configured_python), "conda"),
    file.path(dirname(dirname(configured_python)), "bin", "conda")))
  env_name <- Sys.getenv("YTAB_CONDA_ENV", unset = "ytab-local")
  conda_prefix <- c("run", "--no-capture-output", "-n", env_name, "python")
  for (conda in conda_candidates) {
    if (nzchar(conda) && file.exists(conda) && probe(conda, conda_prefix)) {
      return(list(command = conda, prefix = conda_prefix))
    }
  }
  stop(paste0("Screen Coverage requires Python packages numpy, scipy, and PyYAML. ",
    "The configured Python and the '", env_name, "' conda environment were not usable. ",
    "Set YTAB_PYTHON to a compatible interpreter or YTAB_CONDA_ENV to the environment name."),
    call. = FALSE)
}

sc_poisson_ggplot <- function(n, q, k, theme) {
  xmax <- max(25000,k)*1.04
  points <- data.frame(k=seq(0,xmax,length.out=201));points$coverage<-100*(1-exp(-points$k/n))
  ggplot2::ggplot(points,ggplot2::aes(x=k,y=coverage,group=1))+
    ggplot2::geom_line(colour="#E76F64",linewidth=1.15)+
    ggplot2::geom_hline(yintercept=100*q,colour="#B34A36",linetype="dashed",linewidth=.5)+
    ggplot2::geom_vline(xintercept=k,colour="#2C7FB8",linetype="dotted",linewidth=.5)+
    ggplot2::annotate("point",x=k,y=100*(1-exp(-k/n)),colour="#2C7FB8",size=2.5)+
    ggplot2::scale_y_continuous(limits=c(0,100),labels=function(x)paste0(x,"%"))+
    ggplot2::scale_x_continuous(labels=function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE))+
    ggplot2::labs(title="Theoretical Poisson / Equal-Gene Benchmark",x="Total colonies screened",y="Represented genes recovered (%)")+
    theme+ggplot2::theme(legend.position="none")
}

# The active analytical planner must never inherit the retired, global audit.
screen_coverage_server <- function(input, output, session, repo_root, active, active_project_path, python_bin) {
  if (!requireNamespace("plotly",quietly=TRUE)) return(invisible(NULL))
  inventory <- reactiveVal(NULL)
  error <- reactiveVal("Loading this project's parent-library summary…")
  worker <- new.env(parent=emptyenv())
  worker$process <- NULL
  stop_worker <- function() {
    if (!is.null(worker$process) && worker$process$is_alive()) worker$process$kill()
    worker$process <- NULL
  }
  session$onSessionEnded(stop_worker)
  observeEvent(active_project_path(), {
    stop_worker()
    inventory(NULL)
    error("Loading this project's parent-library summary…")
    worker$project <- active_project_path()
    if (is.null(worker$project) || !length(worker$project) || !nzchar(worker$project)) return()
    work <- tempfile("ytab-project-coverage-")
    dir.create(work)
    worker$output <- file.path(work, "summary.json")
    tryCatch({
      runner <- sc_python_command(python_bin, repo_root)
      worker$process <- processx::process$new(runner$command,
        c(runner$prefix, file.path(repo_root,"scripts/local/ytab_screen_coverage_project_summary.py"),
          "--project-config", worker$project, "--output", worker$output), wd=repo_root,
        stdout=file.path(work,"stdout.log"), stderr=file.path(work,"stderr.log"), cleanup=TRUE)
    }, error=function(e) error(conditionMessage(e)))
  }, ignoreNULL=FALSE)
  observe({
    invalidateLater(300,session)
    if (is.null(worker$process) || worker$process$is_alive()) return()
    value <- if (file.exists(worker$output)) jsonlite::fromJSON(worker$output) else list(error="Could not read this project's library summary.")
    worker$process <- NULL
    if (!identical(worker$project, isolate(active_project_path()))) return()
    if (!is.null(value$error)) error(value$error) else {
      inventory(value)
      error("")
    }
  })
  current <- reactive({
    value <- inventory()
    validate(need(!is.null(value),error()))
    validate(need(identical(value$project_config,active_project_path()),"Loading the active project…"))
    validate(need(value$represented > 0,"No protein-coding genes have observed coding insertions in the selected parent libraries."))
    value
  })
  target <- reactive({
    preset <- input$sc_target_preset %||% "0.95"
    q <- if (identical(preset,"custom")) (input$sc_target %||% 95)/100 else as.numeric(preset)
    max(.000001,min(.999999,q))
  })
  fmt <- function(x) format(round(x),big.mark=",",scientific=FALSE,trim=TRUE)
  output$sc_primary_summary <- renderUI({
    d <- current()
    vals <- list("Reference protein-coding genes"=fmt(d$annotated),
      "Parent-union represented genes"=fmt(d$represented),
      "Observed representation ceiling"=sprintf("%.1f%%",100*d$represented/d$annotated),
      "Observed genomic sites"=fmt(d$unique_sites),"Coding insertion sites"=fmt(d$coding_sites))
    panel_card("Library representation — before colony sampling",
      tags$p(tags$strong(paste("Project:",d$project_id))),
      tags$div(class="ytab-stat-grid",lapply(names(vals),function(name) tags$div(tags$b(vals[[name]]),name))),
      tags$details(tags$summary("Parent libraries and gene-hit summary"),
        tags$p(paste(ytab_sample_labels(d$samples,active()),collapse="; ")),
        tags$p(sprintf("Genes with ≥2 distinct coding sites: %s; ≥3: %s; no observed coding insertion: %s.",fmt(d$ge2),fmt(d$ge3),fmt(d$annotated-d$represented))),
        tags$p(d$definition)))
  })
  output$sc_poisson_target <- renderUI({
    d <- current(); n <- d$represented; q <- target(); k <- ceiling(-n*log1p(-q)); recovered <- n*(1-exp(-k/n))
    panel_card("Poisson target requirement",
      tags$strong(paste(fmt(k),"colonies required")),
      tags$p({digits<-if(q>=.999)2L else 1L;sprintf("Expected recovery: %s of %s represented genes (%.*f%%); %.*f%% of %s reference genes.",fmt(round(recovered)),fmt(n),digits,100*recovered/n,digits,100*recovered/d$annotated,fmt(d$annotated))}),
      withMathJax(tags$p(sprintf("\\[K_q=\\left\\lceil-N\\ln(1-q)\\right\\rceil,\\quad N=%s,\\quad q=%.4f\\]",fmt(n),q))),
      tags$p("The blue marker shows required K; the red dashed line is the target. This is an equal-gene analytical benchmark, not a confidence guarantee."))
  })
  output$sc_frame_primary <- renderUI({v<-ytab_plot_customization_values(input,"sc");ytab_plot_frame(plotly::plotlyOutput("sc_plot_primary",height=v$height),width=v$width)})
  poisson_plot <- function() {
    d <- current(); n <- d$represented; q <- target(); k <- ceiling(-n*log1p(-q))
    sc_poisson_ggplot(n, q, k, sc_theme(ytab_plot_customization_values(input,"sc")))
  }
  output$sc_plot_primary <- plotly::renderPlotly(plotly::config(ytab_ggplotly(poisson_plot(), tooltip = "all"),displaylogo=FALSE))
  download_name <- function(ext) sprintf("%s_screen_planning_poisson_%s%s", active()$project_id, sub("\\.", "p", format(100*target())), ext)
  for (format in c("png", "pdf")) local({
    ext <- format
    output[[paste0("sc_download_", ext)]] <- downloadHandler(filename = function() download_name(paste0(".", ext)),
      content = function(file) ytab_save_ggplot(file, poisson_plot(), width = 9, height = 6))
  })
  output$sc_download_curve <- downloadHandler(filename = function() download_name("_curve.tsv"), content = function(file) {
    d <- current(); n <- d$represented; k <- seq(0, max(25000, ceiling(-n*log1p(-target())))*1.04, length.out = 201)
    utils::write.table(data.frame(colonies = round(k), expected_recovered_genes = n*(1-exp(-k/n)), represented_genes = n,
                                  percent_represented_recovered = 100*(1-exp(-k/n)), percent_reference_recovered = 100*n*(1-exp(-k/n))/d$annotated),
                       file, sep = "\t", row.names = FALSE, quote = FALSE)
  })
}

screen_coverage_ui <- function() {
  if(!requireNamespace("plotly",quietly=TRUE))return(panel_card("Screen Planning",p("Install plotly to display Screen Planning.")))
  panel_card("Screen Planning",ytab_two_column_layout(
    controls=ytab_control_panel("Screen planning",
      conditionalPanel("input.screen_coverage_tabs !== 'mutants'",
      tags$strong("Theoretical Poisson / Equal-Gene Benchmark"),
      selectInput("sc_target_preset","Target coverage",c("90%"=.9,"95%"=.95,"99%"=.99,"99.9%"=.999,"Custom"="custom"),selected="0.95"),
      conditionalPanel("input.sc_target_preset === 'custom'",
        numericInput("sc_target","Custom target gene coverage (%)",95,min=.001,max=99.9999,step=.1)),
      p(class="help-block","The target applies to the represented genes in the parent union. Changing it updates the required colony count immediately."),
      ytab_plot_customization_controls("sc",include_labels=FALSE,default_height="medium",default_width="standard",open=FALSE)),
      conditionalPanel("input.screen_coverage_tabs === 'mutants'",mutant_screen_controls())),
    main=navset_tab(id="screen_coverage_tabs",selected="planner",
      nav_panel("Coverage Planner",value="planner",
        uiOutput("sc_primary_summary"),
        panel_card("The theoretical benchmark: Poisson / 3N",withMathJax(tagList(
          p("\\[\\begin{aligned} C(K) & = 1-e^{-K/N} \\\\ C & = \\text{expected fraction of represented genes recovered} \\\\ K & = \\text{idealized gene-hit opportunities / picks} \\\\ N & = \\text{represented genes in the theoretical benchmark} \\\\ K_{95} & = -N\\ln(0.05) \\approx 2.9957N \\approx 3N \\end{aligned}\\]"))),
          p("Equal gene probabilities are assumed. The model does not account for insertion bias, gene length, abundance skew, growth or bottlenecks. It gives a represented-union benchmark and does not determine pool allocation.")),
        tags$h3("Expected recovery versus total colonies"),uiOutput("sc_frame_primary"),uiOutput("sc_poisson_target"),
        ytab_plot_downloads_ui("sc", c(curve = "Download curve (TSV)")),
        ytab_methods_drawer(
        panel_card("Model, equations and statistical interpretation",withMathJax(tagList(
          p("Each pick is an idealized opportunity to hit one of N equally recoverable represented genes:"),
          p("\\[\\begin{aligned} \\lambda & = K/N \\\\ X & \\sim \\operatorname{Poisson}(\\lambda) \\\\ P(X=x) & = e^{-\\lambda}\\lambda^x/x! \\\\ P(X\\ge1) & = 1-e^{-K/N} \\\\ E[G] & = N(1-e^{-K/N}) \\\\ K_q & = \\left\\lceil-N\\ln(1-q)\\right\\rceil \\\\ K_{95} & \\approx 3N \\end{aligned}\\]"),
          p("A target of 95% means 95% expected recovery of the represented-gene union. It is not a 95% probability of recovering every gene or a confidence guarantee."),
          p("Library representation is the fraction of reference protein-coding genes with observed coding insertions. Genes absent from the observed union cannot be recovered by increasing K; reference recovery is limited by this observed representation ceiling."),
          p("The model assumes independent sampling and equal gene probabilities. It does not estimate pool allocation or fit gene length, insertion architecture, read abundance, growth, bottlenecks or loss-of-function probability.")))),
        panel_card("References",
          p(tags$a(href="https://doi.org/10.1038/ncomms13270","Baym et al. (2016). Knockout Sudoku.")),
          p(tags$a(href="https://doi.org/10.3389/fmicb.2019.00280","Conte et al. (2019). Poisson coverage and the approximately 3N rule."))))),
      nav_panel("Mutant Identification",value="mutants",mutant_screen_ui()))))
}
