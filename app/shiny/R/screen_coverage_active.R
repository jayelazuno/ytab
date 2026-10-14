# Preserve the original planner/server and retire only the advanced UI.
screen_coverage_advanced_ui <- screen_coverage_ui
screen_coverage_advanced_server <- screen_coverage_server

screen_coverage_ui <- function() {
  if(!requireNamespace("plotly",quietly=TRUE))return(panel_card("Screen Coverage",p("Install plotly to display Screen Coverage.")))
  panel_card("Screen Coverage",ytab_two_column_layout(
    controls=ytab_control_panel("Screen planning",
      conditionalPanel("input.screen_coverage_tabs !== 'mutants'",
      tags$strong("Theoretical Poisson / Equal-Gene Benchmark"),
      selectInput("sc_target_preset","Target coverage",c("90%"=.9,"95%"=.95,"99%"=.99,"99.9%"=.999,"Custom"="custom"),selected="0.95"),
      conditionalPanel("input.sc_target_preset === 'custom'",
        numericInput("sc_target","Custom target gene coverage (%)",95,min=.001,max=99.9999,step=.1)),
      p(class="help-block","The target applies to the represented genes in the parent union. Changing it updates the required colony count immediately."),
      ytab_plot_customization_controls("sc",include_labels=FALSE,default_height="medium",default_width="standard")),
      conditionalPanel("input.screen_coverage_tabs === 'mutants'",mutant_screen_controls())),
    main=navset_tab(id="screen_coverage_tabs",selected="planner",
      nav_panel("Screen Planner",value="planner",
        uiOutput("sc_primary_summary"),
        panel_card("The theoretical benchmark: Poisson / 3N",withMathJax(tagList(
          p("\\[\\begin{aligned} C(K) & = 1-e^{-K/N} \\\\ C & = \\text{expected fraction of represented genes recovered} \\\\ K & = \\text{idealized gene-hit opportunities / picks} \\\\ N & = \\text{represented genes in the theoretical benchmark} \\\\ K_{95} & = -N\\ln(0.05) \\approx 2.9957N \\approx 3N \\end{aligned}\\]"))),
          p("Equal gene probabilities are assumed. The model does not account for insertion bias, gene length, abundance skew, growth or bottlenecks. It gives a represented-union benchmark and does not determine pool allocation.")),
        tags$h3("Expected recovery versus total colonies"),uiOutput("sc_frame_primary"),uiOutput("sc_poisson_target")),
      nav_panel("Mutant Identification",value="mutants",mutant_screen_ui()),
      nav_panel("Methods & References",value="methods",
        panel_card("Model, equations and statistical interpretation",withMathJax(tagList(
          p("Each pick is an idealized opportunity to hit one of N equally recoverable represented genes:"),
          p("\\[\\begin{aligned} \\lambda & = K/N \\\\ X & \\sim \\operatorname{Poisson}(\\lambda) \\\\ P(X=x) & = e^{-\\lambda}\\lambda^x/x! \\\\ P(X\\ge1) & = 1-e^{-K/N} \\\\ E[G] & = N(1-e^{-K/N}) \\\\ K_q & = \\left\\lceil-N\\ln(1-q)\\right\\rceil \\\\ K_{95} & \\approx 3N \\end{aligned}\\]"),
          p("A target of 95% means 95% expected recovery of the represented-gene union. It is not a 95% probability of recovering every gene or a confidence guarantee."),
          p("Library representation is the fraction of reference protein-coding genes with observed coding insertions. Genes absent from the observed union cannot be recovered by increasing K; reference recovery is limited by this observed representation ceiling."),
          p("The model assumes independent sampling and equal gene probabilities. It does not estimate pool allocation or fit gene length, insertion architecture, read abundance, growth, bottlenecks or loss-of-function probability.")))),
        panel_card("References",
          p(tags$a(href="https://doi.org/10.1038/ncomms13270","Baym et al. (2016). Knockout Sudoku.")),
          p(tags$a(href="https://doi.org/10.3389/fmicb.2019.00280","Conte et al. (2019). Poisson coverage and the approximately 3N rule.")))))))
}
