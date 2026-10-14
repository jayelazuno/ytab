# Screen Coverage is isolated from analysis jobs and their polling reactives.
sc_method <- function(kind) {
  content <- switch(kind,
    poisson = tagList(
      p("Expected gene representation under independent, equally likely sampling. N is the selected gene denominator; k is colonies picked."),
      p("\\(\\lambda=k/N,\\quad X\\sim\\mathrm{Poisson}(\\lambda),\\quad P(X=x)=e^{-\\lambda}\\lambda^x/x!\\)"),
      p("\\(C(k)=1-e^{-k/N},\\quad E[G]=NC(k),\\quad E[M]=Ne^{-k/N},\\quad k=\\lceil-N\\ln(1-C)\\rceil\\)"),
      p("This is an idealized baseline, not a guarantee for a real library. Baym et al. 2016; Conte et al. 2019.")),
    redundancy = tagList(
      p("\\(P(X\\ge r)=1-\\sum_{x=0}^{r-1}e^{-\\lambda}\\lambda^x/x!\\)"),
      p("\\(P(X\\ge1)=1-e^{-\\lambda};\\quad P(X\\ge2)=1-e^{-\\lambda}(1+\\lambda);\\quad P(X\\ge3)=1-e^{-\\lambda}(1+\\lambda+\\lambda^2/2)\\)"),
      p("The inverse uses a bracketed numerical root, rounded upward. Redundancy counts sampled colonies; repeated recovery of one allele is possible and does not establish independent knockout alleles.")),
    empirical = tagList(
      p("This gene-hit occupancy model is extended to independent pools. In the two Hermes empirical models, K means total insertion-positive colony/mutant picks. Each pool is sampled from its full observed insertion-coordinate catalogue P_j. Coding coordinates can recover a protein-coding gene; noncoding/intergenic/unassigned coordinates still consume picks."),
      p("For pool j, M_j is the number of full-catalogue coordinates, C_j is the coding subset, and m_gj counts coordinates assigned to gene g. Equal-Insertion sampling assigns probability 1/M_j to every coordinate. Observed-Abundance Sensitivity normalizes each read observation by its own usable depth, averages normalized coordinate support within pool, then normalizes over all of P_j."),
      p("\\[\\begin{aligned}P(g\\text{ hit in }j) & = 1-(1-m_{gj}/M_j)^{K_j} \\\\ p_{ij}^{abundance} & = a_{ij}/\\sum_{u\\in P_j}a_{uj} \\\\ P(g\\text{ hit in }j;\\,abundance) & = 1-\\left(1-\\sum_{i\\in g}p_{ij}^{abundance}\\right)^{K_j}\\end{aligned}\\]"),
      p("A gene with more observed coding insertion coordinates has greater recovery probability, while the represented-union denominator remains the observed 4,379-gene target and reference recovery uses the 5,220-gene annotation. Gene length and Hermes target-site preference are not fitted again; their realized effects are reflected in the observed catalogue and m_gj."),
      p("Read abundance is a sensitivity proxy, not measured viable-mutant frequency. Noncoding insertion coordinates remain in both empirical sampling universes and consume K without recovering protein-coding genes. Exact binomial marginals are not a Poisson model.")),
    simulation = tagList(
      p("\\((X_{1j},\\ldots,X_{M_jj})\\sim\\mathrm{Multinomial}(K_j;p_{1j},\\ldots,p_{M_jj})\\), with each independent pool j sampled separately; coding-coordinate probabilities can be summed by gene."),
      p("Hermes simulations sample insertion-positive picks with replacement from each full P_j. Noncoding/unassigned outcomes consume picks. Pool-specific outcomes are combined only after sampling; shared genomic coordinates are deduplicated for distinct-coordinate alleles."),
      p("Seeded runs report mean, SD and 2.5th–97.5th percentiles. Uniform, equal-insertion, H2O2, Zn and rarefaction use seed offsets 0–4; allele models use offsets 5–7, recorded in the output. Intervals describe sampling variability conditional on assumed probabilities, not biological uncertainty. Simulation targets are first tested grid points, not precise inverses. Baym et al. 2016.")),
    alleles = tagList(
      p("Distinct-allele recovery samples insertion coordinates with replacement, records which coordinates have ever been seen, then counts distinct observed coordinates per gene. Repeated recovery of a coordinate never increases the allele count."),
      p("\\[\\begin{aligned}D_i(k) & = \\sum_{a\\in i}\\mathbf{1}(a\\text{ observed in }k\\text{ picks}) \\\\ C_{allele,r}(k) & = N^{-1}\\sum_i\\mathbf{1}(D_i(k)\\ge r)\\end{aligned}\\]"),
      p("This differs from colony-representation redundancy. Catalogue counts bound recoverable allele redundancy. Coordinates are putative distinct alleles, not proof of independent clones or functional knockouts. The theoretical equal-gene model has no allele identities and therefore no invented allele-level curve.")),
    rarefaction = tagList(
      p("\\(X_i\\sim\\mathrm{Hypergeometric}(M,m_i,k),\\quad E[C_r]=N^{-1}\\sum_iP(X_i\\ge r)\\)"),
      p("M is the distinct genomic-site catalogue size and m_i is the number assigned to gene i. Random permutations subsample this catalogue without replacement; the hypergeometric mean provides an independent check. Coordinates shared across pools occur once in the union catalogue."),
      p("This diagnoses insertion-site representation, not viable-colony recovery. Its x-axis and target units are sites, so it is deliberately plotted separately from colony curves. No experimental catalogue of individually isolated colonies is available for colony rarefaction. Concept follows Baym et al. 2016.")),
    ceiling = tagList(
      p("\\(C_{max}=\\#\\{i:p_i>0\\}/N\\)"),
      p("Library representation describes which genes are present; screen sampling coverage describes their expected recovery after k picks. A represented-only denominator conditions on genes already present. Absent genes cannot be recovered by increasing k. The ceiling is based on observed data, not proof of biological absence.")),
    plate = tagList(p("Required plates = ceiling(desired colonies / colonies per plate). Capacity = colonies per plate × plates. Coverage uses the same selected universe and probability model; plating efficiency and colony viability are not measured.")))
  tags$details(class="ytab-technical-details",tags$summary("Method and equations"),withMathJax(content))
}

sc_plot_ui <- function(id, method) tagList(uiOutput(paste0("sc_frame_",id)),sc_method(method))

sc_pool_method<-function()tags$details(class="ytab-technical-details",tags$summary("Primary model assumptions"),withMathJax(tagList(
  p("Pools are independently constructed insertion libraries on the same genetic background. For empirical models K=Σ_jK_j total insertion-positive colony/mutant picks; each pool is sampled from its full observed catalogue P_j before pool outcomes are combined."),
  p("\\[\\begin{aligned}K & = \\sum_{j=1}^{J}K_j \\\\ R_{global} & = \\bigcup_{j=1}^{J}R_j \\\\ G_{global} & = \\bigcup_{j=1}^{J}G_j\\end{aligned}\\]"),
  p("Coding coordinates are the subset that can recover protein-coding genes; noncoding/intergenic/unassigned insertion-positive picks consume budget but do not recover those genes. A repeated coding coordinate across pools is one unique coordinate allele globally, though provenance retains pool-derived occurrences."),
  p("Complexity-normalized allocation uses full genomic coordinate counts M_j: w_j=M_j/sum M_j and K_j approximately K w_j. Balanced allocation uses each pool's whole-library gene-saturation curve. Equal K/J is retained as a comparison. Manual counts describe a chosen plan."),
  p("Observed-abundance sensitivity normalizes each sequencing observation by its depth and combines observations only within their pool. Read abundance is an imperfect proxy for viable-colony abundance; PCR, extraction, mapping, growth and sequencing affect it."),
  p("A coding insertion does not guarantee complete loss of function. Individual colonies may contain multiple linked insertions not resolved by this model. Simulation intervals quantify sampling variability under the stated pool probabilities; they do not include all biological uncertainty."))))

sc_math_equations<-function(){
  math_section<-function(title,body)tags$details(class="ytab-technical-details",tags$summary(title),body)
  panel_card("Math & equations",withMathJax(tagList(
    math_section("A. Theoretical equal-gene Poisson benchmark",tagList(
      p("For N represented genes and K colony picks:"),
      p("\\[\\begin{aligned}\\lambda & = K/N \\\\ P(X\\ge1) & = 1-e^{-K/N} \\\\ E[G] & = N[1-e^{-K/N}]\\end{aligned}\\]"),
      p("For target fraction q, the colony requirement is:"),
      p("\\[K_q=-N\\ln(1-q)\\]"),
      p("At 95%, \\(K_{0.95}=-N\\ln(0.05)=2.9957N\\approx3N\\), the origin of the conventional 3N rule."))),
    math_section("B. Independent pools and observed complexity",tagList(
      p("For J independently constructed pools, each pool retains its full observed insertion-coordinate catalogue P_j; C_j is the coding subset:"),
      p("\\[\\begin{aligned}j & = 1,\\ldots,J \\\\ M_j & = |P_j| \\\\ C_j & = \\text{coding-coordinate subset of }P_j \\\\ K & = \\sum_{j=1}^{J}K_j \\quad\\text{insertion-positive picks}\\end{aligned}\\]"),
      p("Complexity-normalized allocation:"),
      p("\\[\\begin{aligned}w_j & = M_j/\\sum_{l=1}^{J}M_l \\\\ K_j & \\approx K w_j \\\\ \\lambda_j & = K_j/M_j\\end{aligned}\\]"),
      p("Largest-remainder apportionment floors quotas, then assigns remaining picks by descending fractional remainder; pool order breaks exact ties. These are allocation weights, not physical mixture or abundance estimates."))),
    math_section("C. Insertion-site and within-pool gene recovery",tagList(
      p("Equal-Insertion draws each coordinate in P_j with probability 1/M_j. Coding coordinates form subset C_j; a noncoding/intergenic/unassigned outcome consumes one pick without recovering a protein-coding gene."),
      p("If gene g has m_gj coding coordinates in pool j:"),
      p("\\[P(g\\text{ recovered in }j)=1-(1-m_{gj}/M_j)^{K_j}\\]"),
      p("For abundance sensitivity, let a_ij be the mean depth-normalized support for coordinate i across parent/control observations in pool j:"),
      p("\\[\\begin{aligned}p_{ij} & = a_{ij}/\\sum_{u\\in P_j}a_{uj} \\\\ q_{gj} & = \\sum_{i\\mapsto g}p_{ij} \\\\ P(g\\text{ recovered in }j) & = 1-(1-q_{gj})^{K_j}\\end{aligned}\\]"),
      p("For abundance sensitivity, probabilities are normalized over every coordinate in P_j. CDS length, m_gj, coordinates per kb, and the number of pools representing each gene are exposed as diagnostics. Gene length and Hermes +2/+7 target-site preference are not fitted again because realized effects are reflected in the observed catalogue. Bottleneck estimates remain QC context rather than an added probability layer; current estimates are around three or fewer reads in most pools, with some slightly higher."))),
    math_section("D. Global recovery after separate pool sampling",tagList(
      p("For a gene accessible in multiple pools, the explanatory independent-draw probability is:"),
      p("\\[\\begin{aligned}P(g\\text{ missed globally}) & = \\prod_{j=1}^{J}P(g\\text{ missed in }j) \\\\ P(g\\text{ recovered globally}) & = 1-\\prod_{j=1}^{J}P(g\\text{ missed in }j)\\end{aligned}\\]"),
      p("The Monte Carlo implementation samples each pool's full P_j and unions recovered coordinates/genes only after pool-specific sampling."))),
    math_section("E. Balanced saturation-aware allocation",tagList(
      p("Pool j is evaluated on its own represented-gene denominator N_j:"),
      p("\\[S_j(k)=E[G_j(k)]/N_j\\]"),
      p("The implemented objective maximizes the minimum pool saturation fraction under an integer total budget:"),
      p("\\[\\max_{K_1,\\ldots,K_J}\\min_j S_j(K_j)\\quad\\text{subject to }K_j\\ge1,\\;\\sum_jK_j=K\\]"),
      p("Binary search finds the largest feasible common saturation threshold; residual picks are apportioned by largest remainder in proportion to M_j. The objective retains every selected pool and is not a global-gene-count optimizer."))),
    math_section("F. Equal allocation and recovery denominators",tagList(
      p("Equal allocation uses deterministic largest remainders around:"),p("\\[K_j\\approx K/J\\]"),
      p("Observed representation is represented genes divided by reference protein-coding genes. Parent-union recovery and whole-reference recovery use their respective active-catalogue denominators; the latter cannot exceed the observed representation ceiling. Recovery of a target fraction of the represented union is not the same as recovery of that fraction of the whole reference."))),
    math_section("G. Advanced immediate-gene optimizer",p("The retained comparison maximizes expected immediate global ≥1-gene recovery. It can concentrate picks in the currently richest pool, so it is not the default allocation for an exhaustive independent-pool screen.")),
    math_section("Monte Carlo sampling",tagList(
      p("Each replicate draws K_j insertion-positive colony/mutant picks with replacement from pool j's full insertion-coordinate distribution. A gene is hit if a sampled coding coordinate maps to it; noncoding/intergenic/unassigned coordinates consume picks. Pool-specific outcomes are combined only after sampling. The target is expected mean recovery, not a confidence guarantee."),
      p("Repeated recovery of a coordinate counts once toward distinct-coordinate alleles. The same coordinate in two independent pools is two pool-derived occurrences but one unique genomic coordinate. Read abundance is a sequencing proxy, not a measured colony frequency. Mahmutovic et al.'s growth/bottleneck model is not implemented.")))
  )))
}

screen_coverage_ui<-function(){
  if(!requireNamespace("plotly",quietly=TRUE))return(panel_card("Screen Coverage",p("Install plotly to display Screen Coverage.")))
  panel_card("Screen Coverage",ytab_two_column_layout(
    controls=ytab_control_panel("Screen planning",
      selectInput("sc_model_mode","Coverage model",c(
        "Theoretical Poisson / Equal-Gene Benchmark"="poisson",
        "Hermes Monte Carlo — Equal Insertion"="equal",
        "Hermes Monte Carlo — Observed-Abundance Sensitivity"="abundance"),selected="poisson"),
      selectInput("sc_target_preset","Target coverage",c("90%"=.9,"95%"=.95,"99%"=.99,"99.9%"=.999,"Custom"="custom"),selected="0.95"),
      conditionalPanel("input.sc_target_preset === 'custom'",numericInput("sc_target","Custom target gene coverage (%)",95,min=.001,max=99.9999,step=.1)),
      uiOutput("sc_model_controls"),
  p(class="help-block","Pools are independently constructed insertion libraries on the same parental genetic background. Empirical models sample each pool's full observed insertion-coordinate catalogue; coding coordinates determine gene recovery, and complexity weights use total genomic coordinate counts. Abundance is a separate sensitivity model."),
      ytab_plot_customization_controls("sc",include_labels=FALSE,default_height="medium",default_width="standard"),
      tags$details(tags$summary("Downloads"),downloadButton("sc_download_pool_json","All pool-model results & provenance"),downloadButton("sc_download_pool_allocations","Pool allocations CSV"),downloadButton("sc_download_pool_predictions","Predictions CSV"),downloadButton("sc_download_pool_alleles","Independent alleles CSV"))),
    main=navset_tab(id="screen_coverage_tabs",selected="planner",
      nav_panel("Screen Planner",value="planner",
        uiOutput("sc_primary_summary"),
        conditionalPanel("input.sc_model_mode === 'poisson'",panel_card("The theoretical benchmark: Poisson / 3N",withMathJax(tags$div(style="font-size:1.05em;line-height:1.7;margin:6px 0 12px;",p("\\[\\begin{aligned} C(K) & = 1-e^{-K/N} \\\\ C & = \\text{target fraction of represented genes} \\\\ K & = \\text{colonies screened} \\\\ N & = \\text{represented genes in the theoretical benchmark} \\\\ K_{95} & = -N\\ln(0.05) \\approx 2.9957N \\approx 3N \\end{aligned}\\]"))),
        p("Equal gene probabilities are assumed. The model does not account for insertion bias, gene length, abundance skew, growth or bottlenecks. It gives a represented-union benchmark and does not determine pool allocation."))),
        conditionalPanel("input.sc_model_mode !== 'poisson'",panel_card("Hermes empirical Monte Carlo",p("Each independent pool is sampled from its full observed insertion-coordinate catalogue. K is total insertion-positive colony/mutant picks. Coding coordinates can recover protein-coding genes; noncoding/intergenic/unassigned coordinates consume picks but do not recover those genes."),
          conditionalPanel("input.sc_model_mode === 'equal'",withMathJax(tags$div(style="font-size:1.02em;line-height:1.65;margin-top:8px;",p("\\[\\begin{aligned} K & = \\sum_{j=1}^{J} K_j \\\\ M_j & = |P_j| \\quad\\text{full observed coordinates in pool }j \\\\ C_j & \\subseteq P_j \\quad\\text{coding/gene-assigned coordinates} \\\\ m_{gj} & = \\text{coding coordinates assigned to gene }g \\\\ p_{ij} & = 1/M_j \\quad(i\\in P_j) \\\\ P(g\\text{ recovered in pool }j) & = 1-\\left(1-\\frac{m_{gj}}{M_j}\\right)^{K_j} \\\\ P(g\\text{ recovered globally}) & = 1-\\prod_{j=1}^{J} \\left(1-\\frac{m_{gj}}{M_j}\\right)^{K_j} \\end{aligned}\\]")))),
          conditionalPanel("input.sc_model_mode === 'abundance'",p("Observed sequencing abundance is a sensitivity proxy, not a direct measurement of viable colony frequency.")))),
        conditionalPanel("input.sc_model_mode === 'poisson'",tags$h3("Expected recovery versus total colonies")),
        conditionalPanel("input.sc_model_mode !== 'poisson'",tags$h3("Expected recovery versus total insertion-positive picks")),uiOutput("sc_frame_primary"),
        conditionalPanel("input.sc_model_mode === 'poisson'",uiOutput("sc_poisson_target")),
        conditionalPanel("input.sc_model_mode !== 'poisson'",
          tags$h3("Pool allocation at the selected target"),DT::DTOutput("sc_primary_allocation"),
          tags$h3("Expected recovery at the selected target"),DT::DTOutput("sc_primary_prediction"),
          tags$h3("Insertion-positive picks required for the selected target"),DT::DTOutput("sc_primary_targets"),
          tags$details(class="ytab-technical-details",tags$summary("Target requirements and other planner diagnostics"),
        p("Empirical thresholds are expected-coverage crossings under the selected model and allocation. Targets refer to the active parent-union represented genes, not the whole reference. Empirical K includes insertion-positive picks at noncoding/intergenic/unassigned coordinates."),
          selectInput("sc_planner_plot","Explore planner graphics",c(
            "Pool allocation at the evaluated K"="allocation",
            "Allocation across total colony budgets"="allocations_by_k",
            "Sampling depth relative to pool complexity"="sampling_depth",
            "Represented-union vs whole-reference recovery"="union_reference",
            "Distinct insertion-allele recovery"="allele_curves"),selected="allocation"),
          uiOutput("sc_selected_planner_plot")),
          tags$details(class="ytab-technical-details",tags$summary("Distinct insertion-allele recovery"),
            DT::DTOutput("sc_primary_alleles"),sc_method("alleles")),
        tags$details(class="ytab-technical-details",tags$summary("Pool-specific contributions and assumptions"),
          DT::DTOutput("sc_pool_marginal"),sc_pool_method()))),
      nav_panel("Pool Contribution / QC",value="pool_qc",
        p("Pools are independently constructed transposon libraries on the same genetic background. Their insertion complexity and represented-gene repertoires differ. Overlap describes insertion-site and gene content; the same coordinate across pools is one global distinct coordinate."),
        DT::DTOutput("sc_pool_inventory"),tags$h3("Overlap between physical pools"),DT::DTOutput("sc_pool_overlap"),
        selectInput("sc_pool_qc_plot","Explore pool diagnostics",c(
          "Gene saturation by pool"="pool_gene_saturation",
          "Distinct-allele recovery by pool"="pool_allele_saturation",
          "Insertion-coordinate saturation by pool"="pool_site_saturation",
          "Marginal gene recovery by pool"="marginal_curves",
          "Allocation across total budgets"="allocation_curve"),selected="pool_gene_saturation"),
        uiOutput("sc_selected_pool_qc_plot"),
        tags$h3("Observed complexity allocation weights"),DT::DTOutput("sc_complexity_weights"),
        p("Balanced allocation maximizes the minimum pool-specific gene-saturation fraction; complexity allocation approximately equalizes colonies per observed insertion coordinate."),
        tags$h3("Pool-specific sequencing abundance"),DT::DTOutput("sc_pool_abundance"),
        tags$h3("Which pools provide access to each gene?"),DT::DTOutput("sc_gene_access"),DT::DTOutput("sc_pool_probabilities"),downloadButton("sc_download_pool_probabilities","Pool gene probabilities CSV")),
      nav_panel("Library QC",value="representation",uiOutput("sc_inventory"),sc_method("ceiling"),
        selectInput("sc_diagnostic","Library diagnostic",c("Coding length vs insertion count"="length","Gene probabilities (technical selected model)"="probability","Cumulative abundance (technical model)"="lorenz","Read support vs insertion count"="abundance")),
        sc_plot_ui("diagnostic","alleles"),DT::DTOutput("sc_genes"),downloadButton("sc_download_genes_inline","Gene diagnostics CSV")),
      nav_panel("Parent Sequencing",value="parent_sequencing",uiOutput("sc_library_controls"),
        textInput("sc_audit_dir","Audited catalogue directory (blank uses the accepted parent audit)",""),
        p("Repeated sequencing observations characterize each physical pool. Additional detections are newly observed evidence, not newly generated insertion mutants."),
        DT::DTOutput("sc_parent_observations"),uiOutput("sc_frame_discovery"),
        p("Discovery curves average 1,000 random observation orderings. Their intervals describe order dependence, not colony recovery or unseen-library estimates.")),
      nav_panel("Statistical Details",value="statistics",
        tags$h3("Advanced allocation comparisons"),
        p("Maximize immediate >=1-gene recovery optimizes a narrow expected-gene-count objective, not an exhaustive experimental design or independent-allele recovery. Compare exact expected values under identical probability assumptions; no strategy is declared biologically optimal."),DT::DTOutput("sc_pool_baselines"),
        p("These technical models do not supply the primary pool-stratified predictions. Global-union calculations below are hypothetical/conditional diagnostics, not a homogeneous physical-pool screen plan."),
        selectInput("sc_universe","Technical gene denominator",c("Represented protein-coding genes"="informative","All protein-coding genes"="all","Observed ORF-span insertion"="observed",">=2 coding sites"="alleles2",">=3 coding sites"="alleles3")),
        selectInput("sc_basis","Technical diagnostic probability model",c("Hypothetical global equal-site catalogue"="sites","H2O2 observation sensitivity"="h2o2","Zn observation sensitivity"="zn")),
        selectInput("sc_conditioning","Technical sampling denominator",c("Whole catalogue (hypothetical global sampling)"="physical","Conditional on informative coding insertion"="conditional")),
        selectInput("sc_r_preset","Colony redundancy preset",c(">=1"=1,">=2"=2,">=3"=3,">=5"=5,"Custom"="custom")),numericInput("sc_r","Minimum colony representatives r",1,min=1,max=10000),
        numericInput("sc_max_k","Technical curve maximum",100000,min=1),
        checkboxInput("sc_technical_simulations","Also run legacy uniform/global-catalogue and finite-site diagnostics",FALSE),
        sc_plot_ui("theory","poisson"),DT::DTOutput("sc_targets"),sc_plot_ui("redundancy","redundancy"),
        selectInput("sc_pmf_mode","Poisson PMF input",c("Planned colonies"="planned","Set k/N"="depth")),numericInput("sc_depth","k/N",1,min=0,max=1000),uiOutput("sc_pmf_note"),sc_plot_ui("pmf","poisson"),
        sc_plot_ui("comparison","empirical"),DT::DTOutput("sc_comparison_targets"),
        sc_plot_ui("simulation","simulation"),DT::DTOutput("sc_simulation_table"),sc_plot_ui("rarefaction","rarefaction"),DT::DTOutput("sc_rarefaction_targets"),
        downloadButton("sc_download_curves","Technical curve data"),downloadButton("sc_download_rarefaction","Finite-site rarefaction")),
      nav_panel("Math & Equations",value="math",sc_math_equations()),
      nav_panel("Methods & References",value="methods",sc_pool_method(),uiOutput("sc_pool_limitations"),sc_method("poisson"),sc_method("redundancy"),sc_method("alleles"),sc_method("rarefaction"),
        tags$h3("References"),
        p(tags$a(href="https://doi.org/10.1038/ncomms13270","Baym et al. (2016). Knockout Sudoku; Poisson and Monte Carlo precedent.")),
        p(tags$a(href="https://doi.org/10.3389/fmicb.2019.00280","Conte et al. (2019). Poisson coverage and the approximately 3N rule.")),
        p(tags$a(href="https://doi.org/10.1016/j.csbj.2020.03.021","Mahmutovic et al. (2020). Sampling heterogeneity and bottlenecks.")),
        p(tags$a(href="https://doi.org/10.1534/g3.120.401595","Gale et al. (2020). Candida glabrata Hermes library."))))))
}

screen_coverage_legacy_ui <- function() {
  if (!requireNamespace("plotly",quietly=TRUE)) return(panel_card("Screen Coverage",p("Install the plotly R package to use Screen Coverage.")))
  panel_card("Screen Coverage",
    p("Library representation records observed insertion-bearing genes. Screen coverage estimates their recovery after colony picking. These are different quantities; no essentiality classifier defines the default universe."),
    ytab_two_column_layout(
      controls=ytab_control_panel("Screen planning",
        uiOutput("sc_library_controls"),
        selectInput("sc_universe","Gene denominator",c("Represented protein-coding genes (>=1 CDS insertion)"="informative",
          "All annotated protein-coding genes"="all","Genes with an observed ORF-span insertion"="observed",">=2 distinct coding sites observed"="alleles2",">=3 distinct coding sites observed"="alleles3")),
        selectInput("sc_basis","Displayed empirical model",c("Equal-insertion union model"="sites","H2O2 abundance sensitivity (primary proxy)"="h2o2","Zn abundance sensitivity (outgrowth comparison)"="zn")),
        selectInput("sc_conditioning","Sampling denominator",c("All observed sites in the parent union"="physical","Only informative insertions in selected genes (conditional)"="conditional")),
        p(class="help-block","Shared coordinates count once. Treated samples are excluded. A union describes available insertion sites, not a measured physical mixture or viable-colony abundance."),
        numericInput("sc_k","Planned colonies to screen",10000,min=0,step=100),
        selectInput("sc_target_preset","Target preset",c("90%"=.9,"95%"=.95,"99%"=.99,"99.9%"=.999,"Custom"="custom"),selected="0.95"),
        numericInput("sc_target","Target gene coverage (%)",95,min=.001,max=99.9999,step=.1),
        selectInput("sc_r_preset","Minimum representatives",c(">=1"=1,">=2"=2,">=3"=3,">=5"=5,"Custom"="custom")),
        numericInput("sc_r","Minimum colony representatives per gene r",1,min=1,max=10000),
        numericInput("sc_max_k","Curve maximum (colonies)",100000,min=1),
        checkboxInput("sc_simulate","Include Monte Carlo and site rarefaction",TRUE),
        tags$details(tags$summary("Simulation settings"),numericInput("sc_iterations","Iterations",1000,min=2,max=10000),numericInput("sc_seed","Random seed",2026,min=0)),
        actionButton("sc_run","Calculate coverage",class="btn-primary"),actionButton("sc_cancel","Cancel"),
        p(class="help-block","After the first calculation, scientific settings update analytical curves automatically. Use Calculate coverage to refresh simulations. Plot styling never reruns sampling."),
        uiOutput("sc_status"),
        ytab_plot_customization_controls("sc",include_labels=FALSE),
        tags$details(tags$summary("Downloads"),downloadButton("sc_download_settings","Settings and provenance JSON"),downloadButton("sc_download_curves","Saturation curves CSV"),downloadButton("sc_download_targets","Coverage targets CSV"),downloadButton("sc_download_genes","Gene probabilities CSV"),downloadButton("sc_download_simulation","Simulation summary CSV"))),
      main=tagList(uiOutput("sc_inventory"),uiOutput("sc_summary"),navset_tab(id="screen_coverage_tabs",
        nav_panel("Theoretical",value="theoretical",
          tags$h3("Equal-probability Poisson coverage"),sc_plot_ui("theory","poisson"),DT::DTOutput("sc_targets"),
          panel_card("The 3N rule",withMathJax(p("\\(k=-N\\ln(0.05)=2.9957N\\approx3N\\)")),
            p("About 3N picks gives 95% expected gene representation. This is a shortcut derived from Poisson, not a separate model."),
            p("Assumptions: independent sampling, equal gene probabilities, no abundance skew, insertion bias, viability/fitness bias or sampling bottleneck; one representative suffices. These assumptions are not fully satisfied by Hermes libraries."))),
        nav_panel("Redundancy",value="redundancy",tags$h3("Colony representation redundancy"),sc_plot_ui("redundancy","redundancy"),
          selectInput("sc_pmf_mode","Poisson distribution input",c("Use planned colonies"="planned","Set k/N"="depth")),numericInput("sc_depth","k/N",1,min=0,max=1000),
          uiOutput("sc_pmf_note"),sc_plot_ui("pmf","poisson"),tags$h3("Distinct insertion-allele recovery"),sc_plot_ui("alleles","alleles"),DT::DTOutput("sc_alleles_table"),downloadButton("sc_download_alleles","Download distinct-allele simulations")),
        nav_panel("Hermes & comparison",value="comparison",
          p("Click legend entries to toggle curves; hover for numerical values. Both insertion-site and read-abundance estimates are shown."),
          sc_plot_ui("comparison","empirical"),sc_method("simulation"),DT::DTOutput("sc_comparison_targets"),
          downloadButton("sc_download_targets_inline","Download coverage targets")),
        nav_panel("Simulation & rarefaction",value="simulation",
          sc_plot_ui("simulation","simulation"),DT::DTOutput("sc_simulation_table"),
          tags$h3("Distinct insertion-site rarefaction"),sc_plot_ui("rarefaction","rarefaction"),DT::DTOutput("sc_rarefaction_targets"),
          downloadButton("sc_download_rarefaction","Download site rarefaction")),
        nav_panel("Library representation",value="representation",
          sc_method("ceiling"),DT::DTOutput("sc_pools"),
          selectInput("sc_diagnostic","Why the Hermes estimate differs from 3N",c("Cumulative abundance (Lorenz)"="lorenz","Gene sampling probabilities"="probability","Coding length vs informative sites"="length","Read abundance vs informative sites"="abundance")),
          sc_plot_ui("diagnostic","empirical"),tags$h3("Rare and hard-to-sample genes"),
          p("Zero-probability genes are retained with no finite waiting time. Present genes are ordered from lowest estimated recovery probability."),DT::DTOutput("sc_genes"),downloadButton("sc_download_genes_inline","Download gene probabilities")),
        nav_panel("Practical planner",value="planner",tags$h3("Expected coverage at the planned screen size"),DT::DTOutput("sc_planned"),uiOutput("sc_planned_interval"),DT::DTOutput("sc_planned_alleles"),sc_method("alleles"),
          p("Additional colonies below target the selected redundancy r. Predictions are expectations, not guarantees that the stated fraction will be recovered."),
          numericInput("sc_per_plate","Colonies per plate",96,min=1),numericInput("sc_plates","Number of plates",100,min=0),numericInput("sc_desired","Desired colonies",10000,min=0),
          uiOutput("sc_plate_result"),sc_method("plate")),
        nav_panel("Methods & references",value="methods",uiOutput("sc_limitations"),
          sc_method("poisson"),sc_method("redundancy"),sc_method("empirical"),sc_method("simulation"),sc_method("rarefaction"),sc_method("ceiling"),
          tags$h3("References"),
          p(tags$a(href="https://doi.org/10.1038/ncomms13270",target="_blank","Baym et al. (2016). Nature Communications 7:13270.")),
          p(tags$a(href="https://doi.org/10.3389/fmicb.2019.00280",target="_blank","Conte et al. (2019). Frontiers in Microbiology 10:280.")),
          p(tags$a(href="https://doi.org/10.1016/j.csbj.2020.03.021",target="_blank","Mahmutovic et al. (2020). Computational and Structural Biotechnology Journal 18:791–804.")),
          p(tags$a(href="https://doi.org/10.1534/g3.120.401595",target="_blank","Gale et al. (2020). G3 10:3859–3870.")),
          p("Literature context only: Gale reported 513,123 unique sites across three pools (164,000–335,000 per pool) and 1,278 likely essential genes after filtering genes with >50% unmappable segments. These are not the computational inputs for this project."))))))
}

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

sc_curve_plot <- function(data,n,title,xlabel="Colonies screened",references=TRUE,reference_n=NULL) {
  data$tooltip <- sprintf("%s<br>%s: %s<br>Expected gene coverage: %.4f%%<br>Expected represented: %.2f / %s<br>Minimum representatives: %s",
    data$method,xlabel,format(data$k,big.mark=",",trim=TRUE),100*data$mean,n*data$mean,n,data$r)
  if(!is.null(reference_n))data$tooltip<-sprintf("%s<br>Colonies: %s<br>Parent-union genes recovered: %.2f / %s (%.4f%%)<br>Reference genes recovered: %.2f / %s (%.4f%%)",
    data$method,format(data$k,big.mark=",",trim=TRUE),n*data$mean,n,100*data$mean,n*data$mean,reference_n,100*n*data$mean/reference_n)
  p <- ggplot2::ggplot(data,ggplot2::aes(x=k,y=100*mean,colour=method,group=method,text=tooltip))
  if(all(c("lower","upper")%in%names(data))) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=100*lower,ymax=100*upper,fill=method),alpha=.13,colour=NA,show.legend=FALSE)
  p <- p + ggplot2::geom_line(linewidth=1) + ggplot2::scale_y_continuous(limits=c(0,100),labels=function(x)paste0(x,"%")) +
    ggplot2::scale_x_continuous(labels=function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE)) +
    ggplot2::labs(title=title,x=xlabel,y=sprintf("Expected coverage of %s selected genes (%%)",n))
  if(references) p <- p + ggplot2::geom_hline(yintercept=c(90,95,99,99.9),linetype="dashed",colour="grey65",linewidth=.3)
  if(!is.null(reference_n))p<-p+ggplot2::aes(linetype=method)+ggplot2::scale_linetype_manual(values=setNames(
    ifelse(unique(data$method)=="Theoretical Poisson / Equal-Gene Benchmark","dashed","solid"),unique(data$method)))
  p
}

sc_is_parent <- function(row) {
  values<-tolower(trimws(as.character(unlist(row[intersect(names(row),c("control_or_treated","library_role","condition","guessed_condition","fitness_role","treatment"))]))))
  values<-values[!is.na(values)]
  if(any(grepl("\\btreated\\b|\\btreatment\\b",values,perl=TRUE)))return(FALSE)
  if(any(values%in%c("parent","control","mock","untreated","parent/control","parental")))return(TRUE)
  name<-tolower(row$sample)
  !grepl("(^|[-_ ])treated($|[-_ ])",name)&&grepl("(^|[-_ ])(parent|control|mock)($|[-_ ])",name)
}

sc_parent_catalogue <- function(repo_root, current_project, audit_dir=NULL) {
  path<-file.path(audit_dir%||%file.path(repo_root,"docs/codex/parent_union_audit"),"parent_union_audit.json")
  if(!file.exists(path))return(data.frame())
  samples<-jsonlite::fromJSON(path)$samples
  audit<-jsonlite::fromJSON(path)
  namespace<-substr(digest::digest(file=audit$reference,algo="sha256"),1,12)
  identities<-labels<-character(nrow(samples))
  for(j in seq_len(nrow(samples))){
    cfg<-yaml::read_yaml(file.path(repo_root,"output/projects",samples$project[j],"config/project.yaml"))
    sheet<-cfg$sample_sheet;if(!grepl("^/",sheet))sheet<-file.path(repo_root,sheet)
    rows<-read.csv(sheet,stringsAsFactors=FALSE,check.names=FALSE);row<-rows[rows$sample==samples$sample[j],,drop=FALSE]
    stopifnot(nrow(row)==1L,sc_is_parent(row))
    field<-function(name,default=""){x<-row[[name]];if(is.null(x)||!length(x)||is.na(x[1])||!nzchar(as.character(x[1])))default else as.character(x[1])}
    pool<-field("pool_id",field("pool"))
    identities[j]<-paste(field("library_id",namespace),field("species"),field("physical_pool_id",paste0("pool",pool)),sep="::")
    labels[j]<-paste0("Pool ",pool)
  }
  data.frame(key=paste(samples$project,samples$sample,sep="::"),sample=samples$sample,
    project_id=samples$project,label=paste(samples$project,samples$sample,sep=" / "),
    physical_pool_id=identities,pool_label=labels)
}

sc_model_label<-function(basis)switch(basis,sites="Equal-insertion union",h2o2="H2O2 abundance sensitivity",zn="Zn abundance sensitivity")

screen_coverage_server <- function(input,output,session,repo_root,active,active_project_path,python_bin) {
  if(!requireNamespace("plotly",quietly=TRUE))return(invisible(NULL))
  result <- reactiveVal(NULL); message <- reactiveVal("Theoretical preview ready. Select a Hermes model to calculate pool-aware recovery.")
  submitted <- reactiveVal(NULL)
  job <- new.env(parent=emptyenv());job$process<-NULL;job$output<-NULL;job$project<-NULL;job$started<-NULL;job$request<-NULL
  revision <- reactiveVal(0L)
  planner_target <- reactive({
    preset<-input$sc_target_preset%||%"0.95"
    value<-if(identical(preset,"custom"))(input$sc_target%||%95)/100 else as.numeric(preset)
    max(.000001,min(.999999,value))
  })
  stop_worker <- function(){if(!is.null(job$process)&&job$process$is_alive())job$process$kill();job$process<-NULL}
  session$onSessionEnded(stop_worker)
  observeEvent(active_project_path(),{stop_worker();result(NULL);message("Theoretical preview ready. Select a Hermes model to calculate pool-aware recovery.")},ignoreNULL=TRUE)
  audit_directory<-reactive({if(!is.null(input$sc_audit_dir)&&nzchar(trimws(input$sc_audit_dir)))trimws(input$sc_audit_dir)else file.path(repo_root,"docs/codex/parent_union_audit")})
  parent_catalogue<-reactive({req(active_project_path());sc_parent_catalogue(repo_root,active_project_path(),audit_directory())})
  output$sc_library_controls <- renderUI({
    samples<-parent_catalogue()
    if(!is.data.frame(samples)||!nrow(samples))return(p("No included parent/control samples available."))
    tagList(selectizeInput("sc_samples","Parent union: projects / pools",choices=setNames(samples$key,samples$label),selected=samples$key,multiple=TRUE),
      p(class="help-block","Accepted parent-union audit. H2O2 parents were thawed and sequenced; Zn parents were grown to log phase. Their union combines evidence of presence, not physical populations. Essentiality predictions are not used."))
  })
  observeEvent(input$sc_target_preset,{if(input$sc_target_preset!="custom")updateNumericInput(session,"sc_target",value=100*as.numeric(input$sc_target_preset))})
  observeEvent(input$sc_r_preset,{if(input$sc_r_preset!="custom")updateNumericInput(session,"sc_r",value=as.numeric(input$sc_r_preset))})
  output$sc_model_controls<-renderUI({
    mode<-input$sc_model_mode%||%"poisson"
    if(identical(mode,"poisson"))return(NULL)
    tagList(
      p(class="help-block","The selected target determines the expected total insertion-positive colony/mutant picks and their pool allocation under this empirical model. Change the target or allocation strategy and the prediction updates."),
      selectInput("sc_allocation_mode","Pool allocation strategy",
        c("Complexity-normalized"="complexity","Balanced saturation-aware"="balanced",
          "Equal allocation baseline"="equal","Manual per-pool counts"="manual"),selected="complexity"),
      uiOutput("sc_manual_controls"),
      tags$details(tags$summary("Optional allocation constraint"),numericInput("sc_min_pool","Minimum insertion-positive picks per selected pool (0 = none)",0,min=0,step=100)),
      tags$details(tags$summary("Monte Carlo settings"),numericInput("sc_iterations","Simulation iterations",1000,min=100,max=10000,step=100),numericInput("sc_seed","Reproducible random seed",2026,min=0,step=1)),
      actionButton("sc_run","Calculate Hermes coverage",class="btn-primary"),actionButton("sc_cancel","Cancel"),uiOutput("sc_status"))
  })
  output$sc_poisson_target<-renderUI({
    path<-file.path(audit_directory(),"parent_union_audit.json")
    validate(need(file.exists(path),"Accepted parent-union audit is required for Screen Coverage."))
    audit<-jsonlite::fromJSON(path)
    n<-as.numeric(audit$union$protein_genes_with_coding_insertion)
    reference<-as.numeric(audit$gene_annotation$protein_coding_genes)
    q<-planner_target();k<-ceiling(-n*log1p(-q))
    recovered<-n*(1-exp(-k/n))
      panel_card("Poisson target requirement",
      tags$div(class="sc-target-result",style="display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;",tags$div(tags$strong(fmt(k)),tags$span(" colonies required")),
        tags$div(tags$strong(fmt(round(recovered))),tags$span(sprintf(" of %s represented genes expected",fmt(n)))),
        tags$div(tags$strong(sprintf("%.2f%%",100*recovered/reference)),tags$span(sprintf(" of %s reference genes",fmt(reference))))),
      withMathJax(tags$div(style="font-size:1.02em;line-height:1.65;margin-top:8px;",
        p(sprintf("\\[\\begin{aligned} K_q & = \\left\\lceil -N\\ln(1-q) \\right\\rceil \\\\ N & = %s \\\\ q & = %.4f \\\\ K_q & = %s\\text{ colonies} \\end{aligned}\\]",fmt(n),q,fmt(k))),
        p("The blue marker shows the required colony count; the red line shows the selected target. This is an equal-gene theoretical benchmark, not a Hermes-specific prediction."))))
  })
  # Technical conditioning never changes the units of the primary pool planner.
  output$sc_manual_controls<-renderUI({
    samples<-parent_catalogue();req(nrow(samples)>0L)
    if(length(input$sc_samples))samples<-samples[samples$key%in%input$sc_samples,,drop=FALSE]
    # Identity is verified independently by the Python worker. These metadata
    # keys let controls adapt to any number of audited physical pools.
    pools<-unique(samples[,c("physical_pool_id","pool_label")])
    if(!identical(input$sc_allocation_mode,"manual"))return(NULL)
    audit<-jsonlite::fromJSON(file.path(audit_directory(),"parent_union_audit.json"));n<-as.numeric(audit$union$protein_genes_with_coding_insertion);total<-ceiling(-n*log1p(-planner_target()));base<-total%/%nrow(pools);rem<-total%%nrow(pools)
    tagList(lapply(seq_len(nrow(pools)),function(j)numericInput(paste0("sc_manual_",j),pools$pool_label[j],input[[paste0("sc_manual_",j)]]%||%(base+as.integer(j<=rem)),min=0,step=100)))
  })
  request <- reactive({
    catalogue<-parent_catalogue();req(nrow(catalogue)>0L)
    selected_sources<-input$sc_samples
    if(is.null(selected_sources)||!length(selected_sources))selected_sources<-catalogue$key
    chosen<-catalogue[match(selected_sources,catalogue$key),,drop=FALSE]
    req(!anyNA(chosen$key))
    pools<-unique(chosen[,c("physical_pool_id","pool_label")])
    manual<-if(identical(input$sc_allocation_mode,"manual"))setNames(lapply(seq_len(nrow(pools)),function(j)input[[paste0("sc_manual_",j)]]%||%2500),pools$physical_pool_id)else NULL
    mode<-input$sc_model_mode%||%"poisson"
    audit<-jsonlite::fromJSON(file.path(audit_directory(),"parent_union_audit.json"))
    n<-as.numeric(audit$union$protein_genes_with_coding_insertion)
    q<-planner_target()
    # The worker computes the empirical target crossing for the selected
    # probability model/allocation. Until that result is available, K is only
    # the theoretical benchmark used to seed practical-grid output.
    k<-if(identical(mode,"poisson"))ceiling(-n*log1p(-q))else if(identical(input$sc_allocation_mode,"manual")&&!is.null(manual))sum(unlist(manual))else ceiling(-n*log1p(-q))
    list(repo_root=repo_root,audit_dir=audit_directory(),selected_sources=as.list(chosen$key),
      universe=input$sc_universe,basis=input$sc_basis,conditioning=input$sc_conditioning,
      k=k,r=input$sc_r,target=q,max_k=input$sc_max_k,model_mode=mode,
      allocation_mode=input$sc_allocation_mode%||%"complexity",manual_allocation=manual,
      minimum_per_pool=if(!is.null(input$sc_min_pool)&&input$sc_min_pool>0)setNames(as.list(rep(input$sc_min_pool,nrow(pools))),pools$physical_pool_id)else NULL,
      technical_simulations=isTRUE(input$sc_technical_simulations),
      simulate=!identical(mode,"poisson"),iterations=input$sc_iterations%||%1000,seed=input$sc_seed%||%2026)
  })
  launch <- function(include_simulations=TRUE,automatic=FALSE,include_discovery=TRUE,planner_only=FALSE){
    if(!is.null(job$process)&&job$process$is_alive()){
      if(automatic)stop_worker()else{showNotification("Coverage calculation is running.");return()}
    }
    req_data<-request()
    if(!include_simulations)req_data$simulate<-FALSE
    job$request<-req_data
    req_data$include_discovery<-include_discovery
    req_data$planner_only<-planner_only
    work<-tempfile("ytab-screen-coverage-");dir.create(work)
    req_file<-file.path(work,"request.json");job$output<-file.path(work,"result.json")
    jsonlite::write_json(req_data,req_file,auto_unbox=TRUE,pretty=TRUE,null="null")
    job$project<-active_project_path();job$started<-Sys.time()
    if(!identical(req_data$model_mode,"poisson"))message(sprintf("Running %s with %s simulations…",if(req_data$model_mode=="equal")"equal-insertion Hermes Monte Carlo"else"observed-abundance sensitivity Monte Carlo",format(req_data$iterations,big.mark=",")))
    stop_worker()
    tryCatch({runner<-sc_python_command(python_bin,repo_root);job$process<-processx::process$new(runner$command,c(runner$prefix,file.path(repo_root,"scripts/local/ytab_screen_coverage.py"),"--request",req_file,"--output",job$output,"--cache-dir",file.path(tempdir(),"ytab-screen-coverage-cache")),
      stdout=file.path(work,"stdout.log"),stderr=file.path(work,"stderr.log"),wd=repo_root,cleanup=TRUE)
      message("Calculating coverage. The app remains available.");revision(revision()+1L)},error=function(e)message(conditionMessage(e)))
  }
  observeEvent(input$sc_run,launch(include_discovery=FALSE,planner_only=TRUE))
  # The Poisson preview is analytical and needs no catalogue worker. Selecting
  # the Equal-Insertion model starts its pool-aware calculation immediately
  # instead of leaving a theoretical-only curve on screen while the user waits.
  observeEvent(input$sc_model_mode,{
    if(identical(input$sc_model_mode,"equal"))launch(include_simulations=TRUE,automatic=TRUE,include_discovery=FALSE,planner_only=TRUE)
  },ignoreInit=TRUE)
  scientific_settings<-debounce(request,750)
  observeEvent(scientific_settings(),{
    if(!identical(input$sc_model_mode%||%"poisson","poisson")&&
       !is.null(result())&&!identical(scientific_settings(),submitted())&&!identical(scientific_settings(),job$request))
      launch(include_simulations=TRUE,automatic=TRUE,include_discovery=FALSE,planner_only=TRUE)
  },ignoreInit=TRUE)
  observeEvent(result(),{
    if(!is.null(result())&&!is.null(submitted())&&!identical(request(),submitted())&&is.null(job$process))
      launch(include_simulations=TRUE,automatic=TRUE,include_discovery=FALSE,planner_only=TRUE)
  },ignoreInit=TRUE)
  observeEvent(input$sc_cancel,{stop_worker();message("Coverage calculation cancelled; previous results retained.")})
  observe({revision();if(is.null(job$process))return();invalidateLater(500,session)
    if(job$process$is_alive())return()
    status<-job$process$get_exit_status();job$process<-NULL
    stderr_path<-file.path(dirname(job$output),"stderr.log")
    stderr_text<-if(file.exists(stderr_path))paste(readLines(stderr_path,warn=FALSE),collapse="\n")else ""
    payload<-if(file.exists(job$output))tryCatch(jsonlite::fromJSON(job$output),error=function(e)list(error=conditionMessage(e)))else list(error=if(nzchar(stderr_text))paste("Worker produced no result:",stderr_text)else "Worker exited without producing a result. Check the Screen Coverage worker log.")
    if(!is.null(payload$error)){message(paste("Calculation failed:",payload$error));return()}
    if(status!=0){message("Calculation failed; previous results retained.");return()}
    if(!identical(job$project,isolate(active_project_path())))return()
    result(payload);submitted(job$request);message(sprintf("Calculation complete in %.1f seconds. Seed %s; %s iterations when simulations enabled.%s",as.numeric(difftime(Sys.time(),job$started,units="secs")),payload$settings$seed,payload$settings$iterations,if(!isTRUE(payload$settings$simulate))" Analytical results updated; Calculate coverage refreshes simulations."else""))
  })
  output$sc_status<-renderUI({tagList(p(message()),if(!is.null(result())&&!identical(request(),submitted()))p(class="text-warning","Controls differ from displayed results. Calculate coverage to apply scientific settings."))})
  current<-reactive({d<-result();validate(need(!is.null(d),"Calculate coverage to display results."));d})
  fmt<-function(x)format(round(x),big.mark=",",scientific=FALSE,trim=TRUE)
  output$sc_summary<-renderUI({d<-current();emp<-d$planned[d$planned$method==paste0(sc_model_label(d$settings$basis)," — exact"),]
    vals<-list("Screenable genes"=fmt(d$n),"Observed represented genes"=fmt(d$n*d$ceiling),"Planned colonies"=fmt(d$settings$k),
      "Poisson coverage (>=1)"=sprintf("%.2f%%",100*d$planned$coverage[1]),"Selected model coverage (>=1)"=sprintf("%.2f%%",100*emp$coverage),"Selected sampling-model ceiling"=sprintf("%.2f%%",100*d$ceiling))
    tagList(p(sprintf("Sampling coverage of %s selected genes; selected model: %s. Library representation above is a separate quantity.",fmt(d$n),sc_model_label(d$settings$basis))),
      if(d$settings$conditioning=="conditional")p(class="text-warning","Conditional estimates count useful picks only; they do not count every colony picked from the physical pool."),
      tags$div(style="display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;",
        lapply(names(vals),function(name)tags$div(class="ytab-release-card",style="padding:12px;",tags$strong(vals[[name]]),tags$p(name)))))})
  make_plot<-function(id,fn){
    output[[paste0("sc_frame_",id)]]<-renderUI({v<-ytab_plot_customization_values(input,"sc");ytab_plot_frame(plotly::plotlyOutput(paste0("sc_plot_",id),height=v$height),width=v$width)})
    output[[paste0("sc_plot_",id)]]<-plotly::renderPlotly({d<-result();v<-ytab_plot_customization_values(input,"sc")
      mode<-input$sc_model_mode%||%"poisson"
      if(id=="primary"&&is.null(d)){
        validate(need(file.exists(file.path(audit_directory(),"parent_union_audit.json")),
          "Accepted parent-union audit is required for Screen Coverage."))
      }
      if(id=="primary"&&identical(mode,"poisson")){
        n<-if(is.null(d))as.integer(jsonlite::fromJSON(file.path(audit_directory(),"parent_union_audit.json"))$union$protein_genes_with_coding_insertion)else pool_data(d)$summary$represented
        target<-planner_target();planned<-ceiling(-n*log1p(-target));target_k<-planned
        xmax<-max(25000,planned,target_k)*1.04;ks<-seq(0,xmax,length.out=101)
        preview<-data.frame(k=ks,mean=1-exp(-ks/n),r=1,method="Theoretical Poisson / Equal-Gene Benchmark")
        point<-data.frame(k=planned,mean=1-exp(-planned/n),text=sprintf("Planned K: %s<br>Expected recovery: %.2f%%",fmt(planned),100*(1-exp(-planned/n))))
        p<-ggplot2::ggplot(preview,ggplot2::aes(x=k,y=100*mean,group=1))+ggplot2::geom_line(colour="#E76F64",linewidth=1.15)+
          ggplot2::geom_hline(yintercept=100*target,colour="#B34A36",linetype="dashed",linewidth=.5)+
          ggplot2::geom_vline(xintercept=planned,colour="#2C7FB8",linetype="dotted",linewidth=.5)+
          ggplot2::geom_point(data=point,ggplot2::aes(x=k,y=100*mean),inherit.aes=FALSE,colour="#2C7FB8",size=2.5)+
          ggplot2::scale_y_continuous(limits=c(0,100),labels=function(x)paste0(x,"%"))+
          ggplot2::scale_x_continuous(limits=c(0,xmax),labels=function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE))+
          ggplot2::labs(title="Theoretical Poisson / Equal-Gene Benchmark",x="Total colonies screened",y="Represented genes recovered (%)",
            subtitle=sprintf("Target: %.1f%% of %s represented genes  |  required K = %s colonies",100*target,fmt(n),fmt(target_k)))+
          sc_theme(v)+ggplot2::theme(legend.position="none")
      }else if(id=="primary"&&(is.null(d)||!identical(request(),submitted()))){
        selected_model<-if(identical(mode,"equal"))"Equal-Insertion"else"Observed-Abundance Sensitivity"
        audit<-jsonlite::fromJSON(file.path(audit_directory(),"parent_union_audit.json"));n<-as.numeric(audit$union$protein_genes_with_coding_insertion)
        xmax<-max(25000,ceiling(-n*log1p(-planner_target())))*1.04
        calculating<-!is.null(job$process)&&job$process$is_alive()
        p<-ggplot2::ggplot(data.frame(k=c(0,xmax),mean=c(NA_real_,NA_real_)),ggplot2::aes(x=k,y=mean))+ggplot2::geom_blank()+
          ggplot2::scale_x_continuous(limits=c(0,xmax),labels=function(x)format(x,big.mark=",",scientific=FALSE,trim=TRUE))+
          ggplot2::scale_y_continuous(limits=c(0,100),labels=function(x)paste0(x,"%"))+
    ggplot2::labs(title=paste(if(calculating)"Calculating Hermes Monte Carlo —"else"Hermes Monte Carlo —",selected_model),x="Gene-relevant mutant picks K",y="Represented genes recovered (%)",
            subtitle=if(calculating)"Sampling each selected independent pool separately; the curve will update when the calculation finishes."else"Select Calculate Hermes coverage to generate the pool-aware result.")+sc_theme(v)+ggplot2::theme(legend.position="none")
      }else{
        d<-current();p<-fn(d)+sc_theme(v)
        if(id%in%c("theory","redundancy","comparison","simulation","alleles")&&identical(d$settings$conditioning,"conditional"))p<-p+ggplot2::labs(x="Conditional coding picks (selected genes only; not all colonies)")
        if(id%in%c("theory","redundancy","comparison","simulation","rarefaction"))p<-p+ggplot2::geom_hline(yintercept=100*d$settings$target,colour="#B34A36",linetype="dotdash",linewidth=.4)
      }
      plotly::config(plotly::layout(plotly::ggplotly(p,tooltip="text"),legend=list(orientation="h",y=-0.24,x=0,font=list(size=10)),margin=list(t=65,b=95)),displaylogo=FALSE,toImageButtonOptions=list(format="svg"))})
  }
  make_plot("theory",function(d){df<-subset(d$curves,method=="Equal-probability Poisson"&r==d$settings$r);marks<-subset(d$inverse,method=="Equal-probability Poisson"&r==d$settings$r)
    sc_curve_plot(df,d$n,"Equal-probability gene coverage")+ggplot2::geom_vline(data=marks,ggplot2::aes(xintercept=colonies),linetype="dotted",colour="#2C7FB8",inherit.aes=FALSE)})
  make_plot("redundancy",function(d){df<-subset(d$curves,method=="Equal-probability Poisson");df$method<-paste0(">=",df$r," representatives");sc_curve_plot(df,d$n,"Expected redundancy of gene representation")})
  make_plot("pmf",function(d){lam<-if(identical(input$sc_pmf_mode,"depth"))input$sc_depth else d$lambda_value;req(is.finite(lam),lam>=0)
    low<-max(0,qpois(.000001,lam));high<-qpois(.999999,lam);xs<-unique(round(seq(low,high,length.out=min(500,high-low+1))))
    df<-data.frame(x=xs,y=dpois(xs,lam));df$text<-sprintf("X = %s<br>P(X) = %.8f",df$x,df$y)
    ggplot2::ggplot(df,ggplot2::aes(x=x,y=y,text=text))+ggplot2::geom_col(width=.6,fill="#2C7FB8",colour="black",linewidth=1)+ggplot2::labs(x="Sampled representatives of a gene, X",y="P(X = x)",title=sprintf("Poisson distribution: lambda = %.3f",lam))})
  output$sc_pmf_note<-renderUI({d<-current();p("At k=N, lambda=1: about 36.79% of genes remain missed. N picks does not imply every gene is represented. The displayed PMF omits extreme tails; large ranges are thinned to 500 discrete x values.")})
  make_plot("comparison",function(d){df<-subset(d$curves,r==d$settings$r);df$lower<-df$upper<-df$sd<-NA_real_
    if(is.data.frame(d$simulations)&&nrow(d$simulations)){ss<-subset(d$simulations,r==d$settings$r);df<-rbind(df[,names(ss)],ss)}
    sc_curve_plot(df,d$n,"Model comparison: gene sampling coverage")})
  make_plot("simulation",function(d){validate(need(is.data.frame(d$simulations)&&nrow(d$simulations),"Enable Monte Carlo, then calculate coverage."));sc_curve_plot(subset(d$simulations,r==d$settings$r),d$n,"Colony sampling: mean and 95% simulation interval")})
  make_plot("alleles",function(d){validate(need(is.data.frame(d$alleles)&&nrow(d$alleles),"Enable Monte Carlo to calculate distinct-allele recovery."));df<-d$alleles;df$method<-paste(df$method,paste0(">=",df$r," distinct alleles"));sc_curve_plot(df,d$n,"Distinct insertion-allele recovery")})
  make_plot("rarefaction",function(d){validate(need(is.data.frame(d$rarefaction)&&nrow(d$rarefaction),"Enable simulations to compute site rarefaction."));df<-subset(d$rarefaction,r%in%c(1,2,3,d$settings$r));df$method<-paste0(">=",df$r," distinct informative sites");sc_curve_plot(df,d$n,"Finite insertion catalogue (without replacement)","Distinct insertion sites subsampled")})
  make_plot("diagnostic",function(d){df<-d$library$genes;abundance_label<-if(d$settings$basis=="zn")"Zn"else"H2O2";df$read_abundance<-if(d$settings$basis=="zn")df$read_abundance_zn else df$read_abundance_h2o2;df$text<-sprintf("%s %s<br>CDS length: %s<br>Informative sites: %s<br>Read support (%s only): %s<br>Estimated p: %.8g",df$gene_id,df$common_name,df$length,df$informative_sites,abundance_label,df$read_abundance,df$p)
    switch(input$sc_diagnostic%||%"lorenz",
      lorenz={vals<-sort(df$p);validate(need(sum(vals)>0,"No positive sampling mass."));z<-data.frame(x=c(0,seq_along(vals)/length(vals)),y=c(0,cumsum(vals)/sum(vals)));z$text<-sprintf("Genes %.2f%%<br>Cumulative target abundance %.2f%%",100*z$x,100*z$y)
        ggplot2::ggplot(z,ggplot2::aes(x=100*x,y=100*y,text=text))+ggplot2::geom_line(colour="#087E8B",linewidth=1)+ggplot2::geom_abline(slope=1,intercept=0,linetype="dashed")+ggplot2::labs(x="Selected genes, sorted by abundance (%)",y="Cumulative abundance within selected genes (%)",title="Library non-uniformity")},
      probability={df<-df[order(df$p),];df$rank<-seq_len(nrow(df));ggplot2::ggplot(df,ggplot2::aes(x=rank,y=p,text=text))+ggplot2::geom_point(colour="#087E8B",size=1)+ggplot2::labs(x="Gene rank by sampling probability",y="Estimated gene sampling probability",title="Rare and abundant genes")},
      length=ggplot2::ggplot(df,ggplot2::aes(x=length,y=informative_sites,text=text))+ggplot2::geom_point(colour="#087E8B",alpha=.8)+ggplot2::labs(x="Coding length (bp)",y="Distinct informative insertion sites",title="Gene length and insertion representation"),
      abundance=ggplot2::ggplot(df,ggplot2::aes(x=informative_sites,y=log10(1+read_abundance),text=text))+ggplot2::geom_point(colour="#087E8B",alpha=.8)+ggplot2::labs(x="Distinct informative insertion sites",y=paste0("log10(1 + ",abundance_label," read support)"),title="Read support and insertion representation"))})
  humanize_columns<-function(df){
    names(df)<-vapply(names(df),function(name){
      label<-gsub("[._]+"," ",name)
      label<-gsub("([a-z])([A-Z])","\\1 \\2",label)
      tools::toTitleCase(label)
    },character(1))
    df
  }
  table<-function(id,fn,placeholder=NULL){output[[id]]<-DT::renderDT({d<-result();if(is.null(d))shiny::req(FALSE);df<-humanize_columns(fn(d));DT::datatable(df,rownames=FALSE,filter="top",options=list(pageLength=10,scrollX=TRUE),escape=TRUE)})}
  target_display<-function(df){data.frame(Method=df$method,Target=paste0(100*df$target,"%"),Minimum.representatives=df$r,Required=ifelse(is.na(df$colonies),df$status,format(df$colonies,big.mark=",",trim=TRUE)),Status=df$status,check.names=FALSE)}
  table("sc_targets",function(d){x<-subset(d$inverse,method=="Equal-probability Poisson"&r==d$settings$r);data.frame(Target=paste0(100*x$target,"%"),Depth=round(x$depth,4),Colonies=x$colonies,Expected.represented=round(x$represented,2),Expected.missed=round(x$missed,2))})
  table("sc_comparison_targets",function(d)target_display(subset(d$inverse,method!="Site catalogue rarefaction" & ((r==1&target%in%c(.9,.95,.99,.999))|(r%in%c(2,3,d$settings$r)&target==.95)))))
  table("sc_rarefaction_targets",function(d){x<-subset(d$inverse,method=="Site catalogue rarefaction");if(!nrow(x))return(data.frame(Status="No site rarefaction calculated"));y<-target_display(x);names(y)[4]<-"Distinct sites required";y})
  table("sc_genes",function(d){df<-d$library$genes;df<-df[order(df$p==0,df$p),];data.frame(Gene=df$gene_id,CDS.length=df$length,Informative.sites=df$informative_sites,H2O2.read.support=df$read_abundance_h2o2,Zn.read.support=df$read_abundance_zn,Estimated.p=signif(df$p,5),Expected.picks.until.first=ifelse(is.na(df$expected_first),"Not recoverable",format(round(df$expected_first),big.mark=",",trim=TRUE)),Recovery.at.k=round(100*df$recovery,3),Colony.redundancy.ge.2=round(100*df$recovery_2,3))})
  table("sc_pools",function(d)d$library$pools)
  table("sc_planned",function(d){x<-subset(d$planned,method=="Equal-probability Poisson"|grepl(" — exact$",method));interval<-vapply(x$method,function(method){if(!is.data.frame(d$simulations))return("Not simulated");name<-if(method=="Equal-probability Poisson")"Uniform Monte Carlo"else sub(" — exact$"," — Monte Carlo",method);ss<-subset(d$simulations,method==name&k==d$settings$k&r==1);if(!nrow(ss))return("Not simulated");sprintf("%.2f–%.2f%%",100*ss$lower,100*ss$upper)},character(1));data.frame(Method=x$method,Colonies=x$k,Expected.genes=round(x$represented,2),Expected.missed=round(x$missed,2),Coverage.percent=round(100*x$coverage,3),Simulation.interval.95=interval,Colony.redundancy.ge2.percent=round(100*x$at_least_2,3),Colony.redundancy.ge3.percent=round(100*x$at_least_3,3),Additional.to.95=x$additional_95,Additional.to.99=x$additional_99)})
  allele_table<-function(d,planned=FALSE){if(!is.data.frame(d$alleles)||!nrow(d$alleles))return(data.frame(Status="Distinct-allele simulation not available"));x<-d$alleles;if(planned)x<-subset(x,k==d$settings$k);if(!nrow(x))return(data.frame(Status=d$allele_simulation_note));data.frame(Method=x$method,Colonies=x$k,Minimum.distinct.alleles=x$r,Mean.genes=round(x$represented,2),Mean.percent=round(100*x$mean,3),SD.genes=round(d$n*x$sd,2),Lower.95.genes=round(d$n*x$lower,2),Upper.95.genes=round(d$n*x$upper,2))}
  table("sc_alleles_table",allele_table)
  table("sc_planned_alleles",function(d)allele_table(d,TRUE))
  table("sc_simulation_table",function(d){if(!is.data.frame(d$simulations))return(data.frame(Status="Simulations not run"));x<-subset(d$simulations,k==d$settings$k);data.frame(Method=x$method,Colonies=x$k,Minimum.representatives=x$r,Mean.percent=round(100*x$mean,3),SD.percentage.points=round(100*x$sd,3),Lower.95.percent=round(100*x$lower,3),Upper.95.percent=round(100*x$upper,3),Expected.represented=round(x$represented,2),Expected.missed=round(x$missed,2))})
  pool_data<-function(d){validate(need(!is.null(d$pool_screen),"Calculate the pool-aware models."));d$pool_screen}
  selected_target_row<-function(d){
    z<-pool_data(d);q<-planner_target()
    method<-if(identical(input$sc_model_mode,"equal"))"Hermes Empirical Monte Carlo — Equal Insertion"else"Hermes Empirical Monte Carlo — Observed-Abundance Sensitivity"
    strategy<-input$sc_allocation_mode%||%"complexity"
    rows<-z$targets[z$targets$method==method&z$targets$strategy==strategy&abs(z$targets$target-q)<1e-7,,drop=FALSE]
    if(nrow(rows)&&is.finite(as.numeric(rows$k[1])))rows[1,,drop=FALSE]else NULL
  }
  selected_target_k<-function(d){
    z<-pool_data(d);row<-selected_target_row(d)
    if(identical(input$sc_allocation_mode,"manual"))return(sum(vapply(seq_along(z$pools$pool_id),function(j)as.numeric(input[[paste0("sc_manual_",j)]]%||%0),numeric(1))))
    if(!is.null(row)&&is.finite(as.numeric(row$k[1])))return(as.numeric(row$k[1]))
    NA_real_
  }
  selected_target_allocation<-function(d){
    z<-pool_data(d);row<-selected_target_row(d)
    if(!is.null(row)&&"allocation"%in%names(row)){
      value<-as.numeric(unlist(row$allocation[[1]]))
      if(length(value)==nrow(z$pools)&&all(is.finite(value)))return(value)
    }
    if(identical(input$sc_allocation_mode,"manual"))return(vapply(seq_along(z$pools$pool_id),function(j)as.numeric(input[[paste0("sc_manual_",j)]]%||%0),numeric(1)))
    NULL
  }
  simple_table<-function(x){
    for(name in names(x))if(is.list(x[[name]]))x[[name]]<-vapply(x[[name]],function(v)paste(unlist(v),collapse=", "),character(1))
    x
  }
  table("sc_pool_inventory",function(d)simple_table(pool_data(d)$pools))
  table("sc_pool_overlap",function(d)pool_data(d)$overlaps)
  table("sc_complexity_weights",function(d){z<-pool_data(d);x<-z$complexity_weights;coding<-x$eligible_coding_sites%||%z$pools$coding_sites;fraction<-x$coding_fraction%||%(coding/x$observed_sites);weight<-x$weight%||%(x$observed_sites/sum(x$observed_sites));data.frame(Pool=x$label,Observed.genomic.coordinates.M_j=x$observed_sites,Coding.coordinates.C_j=coding,Coding.fraction=round(fraction,4),Allocation.weight.percent=round(100*weight,3))})
  table("sc_gene_access",function(d)pool_data(d)$gene_access)
  table("sc_parent_observations",function(d)pool_data(d)$observation_qc)
  table("sc_pool_baselines",function(d){x<-pool_data(d)$baselines
    data.frame(Method=x$method,Strategy=x$strategy_label,K=x$k,
      Exact.expected.genes=round(x$exact_expected_genes,2),
      Percent.parent.union=round(x$exact_expected_genes_union_percent,3),
      Percent.reference=round(x$exact_expected_genes_reference_percent,3),
      Allocation=vapply(x$allocation,function(v)paste(unlist(v),collapse=" / "),character(1)))})
  table("sc_pool_probabilities",function(d){x<-pool_data(d)$probabilities;n<-nrow(x);length_bp<-x$coding_length_bp%||%rep(NA_real_,n);m_gj<-x$coding_coordinates_in_pool%||%rep(NA_real_,n);global_m<-x$coding_coordinates_global%||%rep(NA_real_,n);pool_count<-x$pools_represented%||%rep(NA_real_,n);per_kb<-x$insertion_coordinates_per_kb%||%(m_gj/(length_bp/1000));data.frame(Gene=x$gene_id,Pool=x$pool_id,
    `CDS length (bp)`=length_bp,`Coding coordinates in pool (m_gj)`=m_gj,
    `Coding coordinates across selected pools`=global_m,`Pools represented`=pool_count,
    `Observed coding coordinates per kb`=round(per_kb,4),
    `Equal-insertion per-pick gene probability`=signif(x$p_equal,5),
    `Abundance-sensitivity per-pick gene probability`=signif(x$p_abundance,5))})
  table("sc_pool_abundance",function(d)pool_data(d)$pool_abundance)
  table("sc_pool_marginal",function(d){z<-pool_data(d);x<-subset(z$marginal,strategy==z$settings$primary_strategy)
    pool_name<-if("label"%in%names(x))x$label else z$pools$label[match(x$pool_id,z$pools$pool_id)]
    noncoding<-x$pool_noncoding_sites%||%(x$pool_genomic_sites-x$pool_coding_sites);coding_fraction<-x$pool_coding_fraction%||%(x$pool_coding_sites/x$pool_genomic_sites)
    data.frame(Model=x$method,Strategy=x$strategy,Pool=pool_name,`Insertion-positive picks`=x$colonies,
      Observed.genomic.coordinates=x$pool_genomic_sites,Coding.coordinates=x$pool_coding_sites,
      Noncoding.or.unassigned.coordinates=noncoding,Coding.fraction=round(coding_fraction,4),
      Represented.genes=x$pool_represented_genes,Genes.unique.in.observed.catalogue=x$pool_unique_genes,
      Unique.coding.coordinates=x$pool_unique_coding_alleles,Expected.unique.coordinates=round(x$expected_unique_sites,2),
      Expected.genes=round(x$pool_expected_genes,2),Pool.saturation.percent=round(100*x$pool_saturation_fraction,3),
      Expected.genes.ge1.distinct.coordinate=round(x$pool_genes_ge1_distinct,2),
      Expected.genes.ge2.distinct.coordinates=round(x$pool_genes_ge2_distinct,2),
      Expected.genes.ge3.distinct.coordinates=round(x$pool_genes_ge3_distinct,2),
      K.over.M=x$normalized_sampling_depth,Leave.out.global.genes=round(x$leave_out_additional_genes,2),
      Marginal.genes.per.100=x$expected_gain_next_100)})
  table("sc_primary_targets",function(d){
    z<-pool_data(d);q<-planner_target()
    selected_method<-if(identical(input$sc_model_mode,"equal"))"Hermes Empirical Monte Carlo — Equal Insertion"else"Hermes Empirical Monte Carlo — Observed-Abundance Sensitivity"
    all_targets<-z$targets;x<-all_targets[all_targets$method==selected_method & all_targets$strategy==z$settings$primary_strategy & abs(all_targets$target-q)<1e-7,,drop=FALSE]
    empirical<-data.frame(Method=x$method,Allocation=x$strategy_label,Target.union.percent=paste0(100*x$target,"%"),Target.reference.percent=round(100*x$target*z$summary$represented/z$summary$annotated,3),`Insertion-positive picks K`=x$k,
      Pool.allocation=if("allocation"%in%names(x))vapply(x$allocation,function(v)paste(unlist(v),collapse=" / "),character(1))else "",
      Status=x$status)
    t<-z$theory[abs(z$theory$target-q)<1e-7,,drop=FALSE]
    rbind(data.frame(Method=t$method,Allocation="Not defined",Target.union.percent=paste0(100*t$target,"%"),Target.reference.percent=round(100*t$target*z$summary$represented/z$summary$annotated,3),`Insertion-positive picks K`=t$k,Pool.allocation="Not defined",Status="Analytical Poisson benchmark"),empirical)
  })
  table("sc_primary_prediction",function(d){z<-pool_data(d);target_row<-selected_target_row(d)
    if(is.null(target_row))return(data.frame(Status="Target crossing is unavailable for this allocation; use the manual allocation diagnostics."))
    method<-target_row$method[1];strategy<-target_row$strategy[1];K<-as.numeric(target_row$k[1]);genes<-as.numeric(target_row$exact_expected_genes[1]);
    data.frame(Method=method,Strategy=target_row$strategy_label[1],`Gene-relevant picks K`=K,Expected.genes=round(genes,1),Expected.missed=round(z$summary$represented-genes,1),
      Percent.parent.union=round(100*genes/z$summary$represented,3),Percent.reference=round(100*genes/z$summary$annotated,3),
      Note="Exact expected-coverage target crossing; simulation intervals are reported at evaluated budgets.")
  })
  selected_allocations<-function(d){z<-pool_data(d);allocation<-selected_target_allocation(d);K<-selected_target_k(d)
    if(is.null(allocation)){
      return(data.frame(Method="Selected Hermes model",Pool="Target unavailable",Colonies=NA_real_,Percent.of.total=NA_real_))
    }
    data.frame(method=if(identical(input$sc_model_mode,"equal"))"Hermes Empirical Monte Carlo — Equal Insertion"else"Hermes Empirical Monte Carlo — Observed-Abundance Sensitivity",
      label=z$pools$label,pool_id=z$pools$pool_id,colonies=allocation,fraction=if(K>0)allocation/K else 0,k=K,strategy=input$sc_allocation_mode%||%"complexity")
  }
  table("sc_primary_allocation",function(d){x<-selected_allocations(d);pools<-pool_data(d)$pools;m<-x$observed_sites%||%pools$sites[match(x$pool_id,pools$pool_id)];c<-x$observed_coding_sites%||%pools$coding_sites[match(x$pool_id,pools$pool_id)];data.frame(Method=x$method,Pool=x$label,`Insertion-positive picks K_j`=x$colonies,`Percent of total`=round(100*x$fraction,2),`Observed genomic coordinates M_j`=m,`Coding coordinates C_j`=c,`Sampling depth K_j/M_j`=round(x$colonies/m,5))})
  table("sc_primary_alleles",function(d){z<-pool_data(d);x<-z$alleles
    if(!is.data.frame(x)||!nrow(x))return(data.frame(Status="Run simulations to calculate independent-allele recovery (planned k up to 1,000,000)."))
    x<-subset(x,k==z$settings$k & strategy==z$settings$primary_strategy)
    data.frame(Method=x$method,Distinct.alleles.at.least=x$r,Expected.genes=round(x$represented,1),
      Percent.parent.union=round(x$represented_union_percent,3),Percent.reference=round(x$represented_reference_percent,3),
      Lower.95.genes=x$lower_genes,Upper.95.genes=x$upper_genes)})
  output$sc_primary_summary<-renderUI({
    if(is.null(result())){
      path<-file.path(audit_directory(),"parent_union_audit.json")
      validate(need(file.exists(path),"Select an accepted parent-union audit."))
      audit<-jsonlite::fromJSON(path);n<-audit$gene_annotation$protein_coding_genes;u<-audit$union
      return(panel_card("Library representation — before colony sampling",tagList(
        tags$div(style="display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:8px;",
          lapply(list(c(fmt(n),"Reference protein-coding genes"),c(fmt(u$protein_genes_with_coding_insertion),"Parent-union represented genes"),c(sprintf("%.2f%%",100*u$protein_genes_with_coding_insertion/n),"Observed representation ceiling"),c(fmt(u$unique_sites),"Observed genomic sites"),c(fmt(u$unambiguous_protein_coding_sites),"Coding insertion sites")),function(v)tags$div(class="ytab-release-card",style="padding:9px;",tags$strong(v[1]),tags$div(v[2])))),
        tags$details(class="ytab-technical-details",tags$summary("How to interpret these values"),p("Pools are independently constructed insertion libraries on the same genetic background. Complexity-normalized allocation uses observed site counts, not abundance or a physical mixing ratio. Library representation is separate from colony-screen recovery.")))))
    }
    z<-pool_data(current());s<-z$summary
    panel_card("Library representation — before colony sampling",tagList(
      tags$div(style="display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:8px;",
        lapply(list(c(fmt(s$annotated),"Reference protein-coding genes"),c(fmt(s$represented),"Parent-union represented genes"),c(sprintf("%.2f%%",100*s$reference_representation),"Observed representation ceiling"),c(fmt(s$sites),"Observed genomic sites"),c(fmt(s$coding_sites),"Coding insertion sites")),function(v)tags$div(class="ytab-release-card",style="padding:9px;",tags$strong(v[1]),tags$div(v[2])))),
      tags$details(class="ytab-technical-details",tags$summary("Pool and model context"),
        p(sprintf("%s independently constructed insertion pools are selected. At the selected target, K = %s under %s. Theoretical Poisson is a benchmark; empirical models sample each pool separately and union recovered coordinates afterward.",nrow(z$pools),if(is.finite(selected_target_k(current())))fmt(selected_target_k(current()))else"not reached",z$settings$primary_strategy_label)),
        p("Observed library representation is relative to the active reference annotation, not colony-screen recovery. Observed-abundance results use sequencing reads as an imperfect proxy for viable-colony frequency.")),
      if(!isTRUE(z$settings$simulate))tags$details(class="ytab-technical-details",tags$summary("Simulation details"),p("Displayed expectations are analytical. Enable Monte Carlo and calculate to add sampling intervals and distinct insertion-allele recovery."))))})
  output$sc_selected_planner_plot<-renderUI({id<-input$sc_planner_plot%||%"allocation";allowed<-c("allocation","allocations_by_k","sampling_depth","union_reference","allele_curves");if(!id%in%allowed)id<-"allocation";uiOutput(paste0("sc_frame_",id))})
  output$sc_selected_pool_qc_plot<-renderUI({id<-input$sc_pool_qc_plot%||%"pool_gene_saturation";allowed<-c("pool_gene_saturation","pool_allele_saturation","pool_site_saturation","marginal_curves","allocation_curve");if(!id%in%allowed)id<-"pool_gene_saturation";uiOutput(paste0("sc_frame_",id))})
  make_plot("primary",function(d){
    z<-pool_data(d);target<-planner_target();target_k<-selected_target_k(d)
    allocation<-selected_target_allocation(d);ids<-z$pools$pool_id;labels<-z$pools$label;n_global<-z$summary$represented
    prob<-as.data.frame(z$probabilities);genes<-unique(as.character(prob$gene_id));pcol<-if(identical(input$sc_model_mode,"abundance"))"p_abundance"else"p_equal"
    P<-matrix(0,nrow=length(genes),ncol=length(ids),dimnames=list(genes,ids))
    P[cbind(match(prob$gene_id,genes),match(prob$pool_id,ids))]<-prob[[pcol]]
    # Pool records are the canonical source for coding-coordinate complexity;
    # tolerate older cached fixtures whose weight table only records genomic
    # site counts.
    site_counts<-as.numeric(z$pools$coding_sites%||%z$complexity_weights$eligible_coding_sites)
    apportion<-function(total,weights){total<-as.integer(round(total));weights<-as.numeric(weights);if(!length(weights)||any(!is.finite(weights))||any(weights<0)||sum(weights)<=0)stop("Pool allocation requires finite, nonnegative weights with positive total.");quota<-total*weights/sum(weights);out<-floor(quota);left<-total-sum(out);if(left>0){ix<-order(-(quota-out),seq_along(out))[seq_len(left)];out[ix]<-out[ix]+1};as.numeric(out)}
    pool_denominators<-colSums(P>0)
    pool_sat<-function(j,k){present<-P[,j]>0;if(!any(present))return(0);sum(-expm1(k*log1p(-pmin(P[present,j],1))))/sum(present)}
    balanced_allocation<-function(total){
      if(total<length(ids))return(apportion(total,site_counts))
      depth_for<-function(j,s){lo<-1;hi<-total;while(lo<hi){mid<-floor((lo+hi)/2);if(pool_sat(j,mid)>=s)hi<-mid else lo<-mid+1};lo}
      lo<-0;hi<-1
      for(iter in seq_len(32)){mid<-(lo+hi)/2;need<-vapply(seq_along(ids),function(j)depth_for(j,mid),numeric(1));if(sum(need)<=total)lo<-mid else hi<-mid}
      out<-vapply(seq_along(ids),function(j)depth_for(j,lo),numeric(1));left<-as.integer(total-sum(out));if(left>0)out<-out+apportion(left,site_counts);out
    }
    allocation_at<-function(total){
      total<-as.integer(round(total))
      if(is.finite(target_k)&&total==target_k&&!is.null(allocation))return(allocation)
      switch(input$sc_allocation_mode%||%"complexity",equal=apportion(total,rep(1,length(ids))),balanced=balanced_allocation(total),
        manual=apportion(total,allocation%||%site_counts),apportion(total,site_counts))
    }
    total_budget_for_pool_target<-function(j){
      low<-0;high<-max(1000,if(is.finite(target_k))target_k else 1000)
      while(pool_sat(j,allocation_at(high)[j])<target&&high<1e10){low<-high;high<-min(1e10,high*2)}
      if(pool_sat(j,allocation_at(high)[j])<target)return(NA_real_)
      while(low<high){mid<-floor((low+high)/2);if(pool_sat(j,allocation_at(mid)[j])>=target)high<-mid else low<-mid+1}
      low
    }
    pool_crossings<-vapply(seq_along(ids),total_budget_for_pool_target,numeric(1))
    poisson_k<-ceiling(-n_global*log1p(-target))
    xmax<-max(c(25000,poisson_k,target_k,pool_crossings),na.rm=TRUE);if(!is.finite(xmax))xmax<-25000
    ks<-unique(sort(c(seq(0,min(xmax,25000),by=1000),if(xmax>25000)exp(seq(log(25000),log(xmax),length.out=45)),poisson_k,pool_crossings[is.finite(pool_crossings)])))
    allocations_by_k<-vapply(ks,allocation_at,numeric(length(ids)))
    if(is.null(dim(allocations_by_k)))allocations_by_k<-matrix(allocations_by_k,nrow=length(ids))
    model_name<-if(identical(input$sc_model_mode,"equal"))"Equal-insertion"else"Observed-abundance sensitivity"
    pool_rows<-do.call(rbind,lapply(seq_along(ids),function(j){
      kj<-allocations_by_k[j,]
      frac<-vapply(seq_along(ks),function(i)pool_sat(j,kj[i]),numeric(1))
      data.frame(k=ks,mean=frac,series=paste(labels[j],"(% of own represented genes)"),pool=labels[j],text=sprintf("%s<br>Total screen budget K=%s<br>Pool allocation K_j=%s<br>Expected genes recovered: %.2f%% of %s represented genes",labels[j],fmt(ks),fmt(kj),100*frac,fmt(pool_denominators[j])),stringsAsFactors=FALSE)
    }))
    theory<-data.frame(k=ks,mean=1-exp(-ks/n_global),series="Theoretical Poisson / Equal-Gene Benchmark",pool="",text=sprintf("Theoretical Poisson benchmark<br>Total K=%s<br>Expected recovery: %.2f%% of %s represented-union genes",fmt(ks),100*(1-exp(-ks/n_global)),fmt(n_global)))
    df<-rbind(theory[,c("k","mean","series","pool","text")],pool_rows)
    markers<-data.frame(k=c(poisson_k,pool_crossings),y=rep(100*target,length(ids)+1),series=c("Theoretical Poisson / Equal-Gene Benchmark",vapply(labels,function(x)paste(x,"(% of own represented genes)"),character(1))),pool=c("Poisson",labels),text=c(sprintf("Poisson target K=%s",fmt(poisson_k)),sprintf("%s reaches %.2f%% of its own represented genes at total K=%s (pool picks=%s)",labels,100*target,fmt(pool_crossings),vapply(seq_along(pool_crossings),function(j)if(is.finite(pool_crossings[j]))fmt(allocation_at(pool_crossings[j])[j])else"not reached",character(1)))))
    markers<-markers[is.finite(markers$k),,drop=FALSE]
    pool_series<-unique(pool_rows$series);pool_colours<-setNames(grDevices::hcl.colors(length(pool_series),"Dark 3"),pool_series)
    ggplot2::ggplot(df,ggplot2::aes(x=k,y=100*mean,colour=series,group=series,text=text))+
      ggplot2::geom_line(data=subset(df,pool!=""),linetype="dotted",linewidth=1.05,alpha=.95)+
      ggplot2::geom_line(data=subset(df,pool==""),linewidth=1.2)+
      ggplot2::geom_hline(yintercept=100*target,linetype="dashed",colour="grey35",linewidth=.5)+
      ggplot2::geom_segment(data=markers,ggplot2::aes(x=k,xend=k,y=0,yend=y,colour=series,text=text),inherit.aes=FALSE,linetype="dotted",linewidth=.65,show.legend=FALSE)+
      ggplot2::geom_point(data=markers,ggplot2::aes(x=k,y=y,colour=series,text=text),inherit.aes=FALSE,size=2.2,show.legend=FALSE)+
      ggplot2::scale_colour_manual(values=c("Theoretical Poisson / Equal-Gene Benchmark"="#E76F64",pool_colours))+
      ggplot2::scale_y_continuous(limits=c(0,100),labels=function(v)paste0(v,"%"))+
      ggplot2::scale_x_continuous(limits=c(0,xmax*1.02),labels=function(v)format(v,big.mark=",",scientific=FALSE,trim=TRUE))+
      ggplot2::labs(x="Total insertion-positive pick budget K",y="Expected recovery (%)",colour=NULL,title=paste("Hermes Monte Carlo —",model_name),
        subtitle=sprintf("Dotted pool curves use each pool’s own gene denominator; Poisson uses the %s-gene global represented union. Vertical markers show total K needed to reach %.2f%%.",fmt(n_global),100*target))
  })
  make_plot("allocation",function(d){x<-selected_allocations(d);x$text<-sprintf("%s<br>%s: %s insertion-positive picks (%.2f%%)",x$method,x$label,fmt(x$colonies),100*x$fraction)
    ggplot2::ggplot(x,ggplot2::aes(x=label,y=colonies,fill=method,text=text))+ggplot2::geom_col(position="dodge",width=.6,colour="black",linewidth=1)+
    ggplot2::labs(x="Independent insertion pool",y="Insertion-positive picks (K_j)",fill=NULL,title="Allocation for the selected total budget")})
  make_plot("allocation_curve",function(d){z<-pool_data(d);x<-subset(z$allocations,strategy%in%c("equal","complexity","balanced","optimized"));x$text<-sprintf("%s<br>%s<br>Total K: %s<br>Pool insertion-positive picks: %s",x$strategy_label,x$label,fmt(x$k),fmt(x$colonies))
    ggplot2::ggplot(x,ggplot2::aes(x=k,y=colonies,colour=strategy_label,group=interaction(strategy_label,label),text=text))+ggplot2::geom_line(linewidth=1)+ggplot2::facet_wrap(~label)+ggplot2::labs(x="Total insertion-positive pick budget K",y="Insertion-positive picks allocated to pool",colour="Allocation",title="Allocation by total insertion-positive pick budget")})
  make_plot("allocations_by_k",function(d){z<-pool_data(d);x<-subset(z$allocations,strategy%in%c("equal","complexity","balanced"));x$text<-sprintf("%s<br>%s<br>K=%s; pool allocation=%s",x$strategy_label,x$label,fmt(x$k),fmt(x$colonies))
    ggplot2::ggplot(x,ggplot2::aes(x=factor(k),y=colonies,fill=label,text=text))+ggplot2::geom_col(colour="black",linewidth=.25)+ggplot2::facet_grid(model~strategy_label)+ggplot2::labs(x="Total insertion-positive pick budget K",y="Insertion-positive picks assigned (stacked by pool)",fill="Independent pool",title="Pool allocation across budgets")})
  make_plot("sampling_depth",function(d){z<-pool_data(d);x<-subset(z$allocations,strategy%in%c("equal","complexity","balanced"));x$text<-sprintf("%s<br>%s<br>K=%s<br>K_j/M_j=%.5f",x$strategy_label,x$label,fmt(x$k),x$normalized_depth)
    ggplot2::ggplot(x,ggplot2::aes(x=factor(k),y=normalized_depth,fill=label,text=text))+ggplot2::geom_col(position="dodge",colour="black",linewidth=.25)+ggplot2::facet_grid(model~strategy_label)+ggplot2::labs(x="Total insertion-positive pick budget K",y="Full-catalogue sampling depth K_j / M_j",fill="Independent pool",title="Sampling depth by pool")})
  make_plot("union_reference",function(d){z<-pool_data(d);x<-subset(z$curves,method!="Theoretical Poisson / Equal-Gene Benchmark"&strategy%in%c("complexity","balanced","equal"));x$series<-paste(x$method,x$strategy_label,sep=" — ");theory<-subset(z$curves,method=="Theoretical Poisson / Equal-Gene Benchmark");theory$strategy_label<-"Theoretical benchmark";theory$series<-theory$strategy_label;x<-rbind(x,theory[,names(x)]);x$percent<-100*x$mean*z$summary$represented/z$summary$annotated;x$lo<-100*x$lower*z$summary$represented/z$summary$annotated;x$hi<-100*x$upper*z$summary$represented/z$summary$annotated;x$text<-sprintf("%s<br>K=%s<br>Represented union: %.2f%%<br>Reference: %.2f%%",x$series,fmt(x$k),100*x$mean,x$percent)
    ggplot2::ggplot(x,ggplot2::aes(x=k,y=percent,colour=series,group=series,text=text))+ggplot2::geom_line(linewidth=.9)+ggplot2::geom_hline(yintercept=100*z$summary$reference_representation,linetype="dashed",colour="grey40")+ggplot2::scale_y_continuous(labels=function(v)paste0(v,"%"))+ggplot2::labs(x="Insertion-positive mutant picks K",y=sprintf("Recovered genes (%% of %s reference genes)",fmt(z$summary$annotated)),title="Recovery relative to the reference",subtitle=sprintf("Observed representation ceiling: %.2f%%",100*z$summary$reference_representation))})
  make_plot("allele_curves",function(d){z<-pool_data(d);model<-if(identical(input$sc_model_mode,"abundance"))"abundance"else"equal";strategy<-input$sc_allocation_mode%||%"complexity"
    # Avoid showing every model/allocation in the planner's focused result. The
    # previous faceting duplicated strategy labels and collided with the legend.
    x<-subset(z$alleles,z$alleles$model==model&z$alleles$strategy==strategy);validate(need(nrow(x)>0,"Run Monte Carlo to display distinct-coordinate allele recovery and intervals."));x$series<-paste0(">=",x$r," distinct coordinates");x$text<-sprintf("%s<br>K=%s<br>Expected genes: %.1f (%.2f%% of represented union)",x$series,fmt(x$k),x$represented,100*x$mean)
    ggplot2::ggplot(x,ggplot2::aes(x=k,y=100*mean,colour=series,group=series,text=text))+ggplot2::geom_line(linewidth=.9)+ggplot2::scale_y_continuous(limits=c(0,100),labels=function(v)paste0(v,"%"))+ggplot2::labs(x="Insertion-positive mutant picks K",y="Represented-gene recovery (%)",colour="Distinct insertion alleles per gene",title="Distinct insertion-coordinate allele recovery",subtitle=paste(unique(x$method),"|",unique(x$strategy_label)))})
  make_plot("pool_gene_saturation",function(d){x<-pool_data(d)$pool_saturation;x$text<-sprintf("%s<br>%s<br>Insertion-positive picks=%s<br>Expected genes: %.1f / %s<br>Distinct alleles: >=1 %.1f; >=2 %.1f; >=3 %.1f",x$method,x$pool_label,fmt(x$k),x$expected_genes,fmt(x$gene_denominator),x$expected_genes_ge1_allele,x$expected_genes_ge2_alleles,x$expected_genes_ge3_alleles)
    ggplot2::ggplot(x,ggplot2::aes(x=k,y=100*gene_fraction,colour=pool_label,group=pool_label,text=text))+ggplot2::geom_ribbon(ggplot2::aes(ymin=100*lower_gene_fraction,ymax=100*upper_gene_fraction,fill=pool_label),alpha=.1,colour=NA,show.legend=FALSE,na.rm=TRUE)+ggplot2::geom_line(linewidth=.9)+ggplot2::facet_wrap(~method)+ggplot2::scale_y_continuous(limits=c(0,100),labels=function(v)paste0(v,"%"))+ggplot2::labs(x="Insertion-positive picks sampled from this pool",y="Pool's represented genes recovered (%)",colour="Independent pool",title="Pool-specific gene saturation (denominator is N_j)")})
  make_plot("pool_allele_saturation",function(d){x<-pool_data(d)$pool_saturation
    vals<-rbind(data.frame(x,threshold="≥1 distinct coordinate",recovered=x$expected_genes_ge1_allele),data.frame(x,threshold="≥2 distinct coordinates",recovered=x$expected_genes_ge2_alleles),data.frame(x,threshold="≥3 distinct coordinates",recovered=x$expected_genes_ge3_alleles))
    vals$text<-sprintf("%s<br>%s<br>k=%s<br>%s genes: %.1f / %s pool-represented genes",vals$method,vals$pool_label,fmt(vals$k),vals$threshold,vals$recovered,fmt(vals$gene_denominator))
    ggplot2::ggplot(vals,ggplot2::aes(x=k,y=100*recovered/gene_denominator,colour=threshold,group=threshold,text=text))+ggplot2::geom_line(linewidth=.8)+ggplot2::facet_grid(method~pool_label)+ggplot2::scale_y_continuous(limits=c(0,100),labels=function(v)paste0(v,"%"))+ggplot2::labs(x="Gene-relevant picks sampled from this pool",y="Pool genes recovered (%)",colour="Distinct alleles per gene",title="Pool-specific recovery of genes with independent insertion coordinates")})
  make_plot("pool_site_saturation",function(d){x<-pool_data(d)$pool_saturation;x$text<-sprintf("%s<br>%s<br>Insertion-positive picks=%s<br>Expected unique genomic coordinates: %.1f / %s (%.2f%%)",x$method,x$pool_label,fmt(x$k),x$expected_unique_sites,fmt(x$site_denominator),100*x$insertion_sites_fraction)
    ggplot2::ggplot(x,ggplot2::aes(x=k,y=100*insertion_sites_fraction,colour=pool_label,group=pool_label,text=text))+ggplot2::geom_line(linewidth=.9)+ggplot2::facet_wrap(~method)+ggplot2::scale_y_continuous(limits=c(0,100),labels=function(v)paste0(v,"%"))+ggplot2::labs(x="Insertion-positive picks from this pool",y="Observed genomic coordinates recovered (%)",colour="Independent pool",title="Pool-specific genomic-coordinate saturation")})
  make_plot("marginal_curves",function(d){x<-pool_data(d)$pool_saturation;x$text<-sprintf("%s<br>%s<br>k=%s<br>Marginal genes per next 100 picks: %.4f",x$method,x$pool_label,fmt(x$k),x$marginal_genes_per_100)
    ggplot2::ggplot(x,ggplot2::aes(x=k,y=marginal_genes_per_100,colour=pool_label,group=pool_label,text=text))+ggplot2::geom_line(linewidth=.9)+ggplot2::facet_wrap(~method,scales="free_y")+ggplot2::labs(x="Gene-relevant picks sampled from this pool",y="Expected additional genes per 100 picks",colour="Independent pool",title="Diminishing marginal gene recovery")})
  make_plot("discovery",function(d){x<-pool_data(d)$discovery;validate(need(nrow(x)>0,"Parent-sequencing discovery is calculated by the full Calculate coverage action; the default planner skips this secondary diagnostic."));x$text<-sprintf("%s<br>Observations: %s<br>Mean: %.1f<br>95%% random-order range: %.1f–%.1f",x$metric,x$observations,x$mean,x$lower,x$upper)
    ggplot2::ggplot(x,ggplot2::aes(x=observations,y=mean,colour=metric,group=metric,text=text))+ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper,fill=metric),alpha=.15,colour=NA)+ggplot2::geom_line(linewidth=1)+ggplot2::facet_wrap(~metric,scales="free_y")+ggplot2::labs(x="Parent sequencing observations incorporated",y="Distinct sites / represented genes",title="Sequencing discovery, not colony-screen coverage")})
  output$sc_pool_limitations<-renderUI({tags$details(class="ytab-technical-details",tags$summary("Scientific limitations and assumptions"),tags$ul(lapply(pool_data(current())$limitations,tags$li)))})
  output$sc_inventory<-renderUI({d<-result()
    if(is.null(d)){
      path<-file.path(audit_directory(),"parent_union_audit.json")
      validate(need(file.exists(path),"Accepted parent-union audit is required for Screen Coverage."))
      audit<-jsonlite::fromJSON(path);u<-audit$union;overlap<-audit$overlap;n<-audit$gene_annotation$protein_coding_genes
      return(tagList(tags$h3("Accepted library representation — before colony sampling"),
        p(sprintf("%s annotated protein-coding genes; %s with >=1 coding insertion; %s with >=2 distinct coding sites; %s with >=3; %s without coding insertions. Observed reference representation: %.2f%%.",fmt(n),fmt(u$protein_genes_with_coding_insertion),fmt(u$protein_genes_with_ge2_sites),fmt(u$protein_genes_with_ge3_sites),fmt(u$protein_genes_without_coding_insertion),100*u$protein_genes_with_coding_insertion/n)),
        p(sprintf("Union: %s genomic sites; %s unambiguous protein-coding CDS sites. Parent sequencing provenance is listed in its own subtab.",fmt(u$unique_sites),fmt(u$unambiguous_protein_coding_sites))),
        p("This is the full accepted parent-union inventory, not predicted colony recovery. Calculate coverage to apply the selected parent subset and sampling model.")))
    }
    l<-d$library;a<-l$representation;tagList(tags$h3("Library representation — before colony sampling"),
    p(sprintf("%s annotated protein-coding genes; %s with >=1 coding insertion; %s with >=2 distinct coding sites; %s with >=3; %s without an observed coding insertion. Reference-level observed representation: %.2f%%.",fmt(a$annotated),fmt(a$represented),fmt(a$ge2),fmt(a$ge3),fmt(a$missing),100*a$reference_fraction)),
    p(sprintf("Coding-represented genes: H2O2 %s; Zn %s; union %s. H2O2-only %s; Zn-only %s; shared %s.",fmt(a$h2o2_genes),fmt(a$zn_genes),fmt(a$represented),fmt(a$h2o2_only_genes),fmt(a$zn_only_genes),fmt(a$shared_genes))),
    p(sprintf("%s distinct genomic union sites; %s unambiguous protein-coding CDS sites. Shared coordinates count once. Pseudogenes and essentiality-classifier filtering are excluded from the primary gene denominator.",fmt(l$unique_sites),fmt(l$informative_sites))),
    p("Representation is evidence of an observed insertion, not the percentage recovered after colony picking. Sampling can approach 100% within the union-represented universe under the equal-insertion model; abundance-model ceilings can be lower when that dataset did not observe a gene."),
    p(l$informative_definition),
    if(!isTRUE(l$audit_matches))p(class="text-warning","A subset of audited parents is selected; representation counts reflect that subset."),
    if(d$settings$conditioning=="conditional")p(class="text-warning","Conditional calculation: k counts coding picks in selected genes only, not all colonies."))})
  output$sc_limitations<-renderUI({d<-current();tagList(p(d$library$informative_definition),tags$ul(lapply(d$library$limitations,tags$li)),p(d$allele_simulation_note),p("Representation derives from the accepted parent-union audit. Provenance downloads include audit hashes and the model-specific seed map."))})
  output$sc_planned_interval<-renderUI({d<-current();if(!is.data.frame(d$simulations))return(p("Monte Carlo interval unavailable: calculate with simulations enabled."));name<-paste0(sc_model_label(d$settings$basis)," — Monte Carlo");x<-subset(d$simulations,method==name&k==d$settings$k&r==d$settings$r);req(nrow(x)==1);p(sprintf("%s at %s picks, >=%s colony representatives: mean %.2f%%; 95%% simulation interval %.2f–%.2f%%. This describes sampling variability conditional on assumed probabilities.",name,fmt(x$k),x$r,100*x$mean,100*x$lower,100*x$upper))})
  output$sc_plate_result<-renderUI({d<-current();req(input$sc_per_plate>0,input$sc_plates>=0,input$sc_desired>=0)
    per<-ceiling(input$sc_per_plate);capacity<-per*floor(input$sc_plates);p<-d$library$genes$p
    theory<-ppois(d$settings$r-1,capacity/d$n,lower.tail=FALSE);emp<-mean(pbinom(d$settings$r-1,capacity,p,lower.tail=FALSE))
    tagList(p(sprintf("Capacity: %s colonies. Plates required for %s desired colonies: %s.",fmt(capacity),fmt(input$sc_desired),fmt(ceiling(input$sc_desired/per)))),p(sprintf("At capacity: Poisson %.2f%%; exact Hermes %.2f%% of the selected %s genes with >=%s representatives.",100*theory,100*emp,fmt(d$n),d$settings$r)))})
  export<-function(id,fn,name){output[[id]]<-downloadHandler(filename=function()paste0("screen_coverage_",name,".csv"),content=function(file)write.csv(fn(current()),file,row.names=FALSE,na=""))}
  export("sc_download_curves",function(d)d$curves,"curves")
  for(id in c("sc_download_targets","sc_download_targets_inline"))export(id,function(d)d$inverse,"targets")
  for(id in c("sc_download_genes","sc_download_genes_inline"))export(id,function(d)d$library$genes,"gene_probabilities")
  export("sc_download_simulation",function(d)d$simulations,"simulations")
  export("sc_download_rarefaction",function(d)d$rarefaction,"site_rarefaction")
  export("sc_download_alleles",function(d)d$alleles,"distinct_alleles")
  export("sc_download_pool_allocations",function(d)d$pool_screen$allocations,"pool_allocations")
  export("sc_download_pool_predictions",function(d)simple_table(d$pool_screen$practical),"pool_predictions")
  export("sc_download_pool_alleles",function(d)d$pool_screen$alleles,"pool_alleles")
  export("sc_download_pool_probabilities",function(d)d$pool_screen$probabilities,"pool_gene_probabilities")
  output$sc_download_pool_json<-downloadHandler(filename=function()"pool_screen_coverage.json",content=function(file)jsonlite::write_json(current()$pool_screen,file,pretty=TRUE,auto_unbox=TRUE,na="null"))
  output$sc_download_settings<-downloadHandler(filename=function()"screen_coverage_provenance.json",content=function(file){d<-current();jsonlite::write_json(list(settings=d$settings,software=d$software,seeds=d$seed_map,sources=d$library$sources,definition=d$library$informative_definition,universes=d$library$universes,representation=d$library$representation,limitations=d$library$limitations),file,pretty=TRUE,auto_unbox=TRUE)})
}
