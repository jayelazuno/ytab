mutant_screen_controls <- function(){
  tagList(tags$h4("Identify individual mutants"),
    radioButtons("ms_input_mode","Sequence input",c("Local file or folder"="path","Upload files"="upload"),inline=TRUE),
    conditionalPanel("input.ms_input_mode !== 'upload'",textAreaInput("ms_paths","Paths (one per line)","",rows=3,placeholder=".seq/.ab1 files, folders or .tar archives")),
    conditionalPanel("input.ms_input_mode === 'upload'",fileInput("ms_files","Sanger files or archive",multiple=TRUE,
      accept=c(".seq",".fa",".fasta",".fna",".fastq",".fq",".ab1",".abi",".tar",".tgz",".gz",".zip"))),
    textInput("ms_condition","Screen condition (optional)",""),textInput("ms_pool","Pool / clone group (optional)",""),
    uiOutput("ms_reference_control"),uiOutput("ms_link_project_control"),uiOutput("ms_fitness_control"),
    tags$details(class="ytab-technical-details",tags$summary("Alignment and clipping settings"),
      numericInput("ms_min_length","Minimum aligned genomic bases",80,min=20,max=5000),
      numericInput("ms_identity","Minimum identity (%)",90,min=50,max=100),
      numericInput("ms_mapq","Minimum MAPQ when available",20,min=0,max=60),
      numericInput("ms_quality","End-trimming Phred threshold (ABI / FASTQ)",20,min=0,max=60),
      numericInput("ms_overlap","Minimum long Hermes overlap (bases)",12,min=8,max=200),
      numericInput("ms_junction_overlap","Hermes overlap for junction-review flag (bp)",20,min=12,max=200),
      numericInput("ms_gene_window","Multiple-CDS resolution window (bp)",100,min=50,max=500),
      numericInput("ms_gene_min_bases","Minimum proximal aligned bases",50,min=20,max=500),
      selectInput("ms_junction_end","Junction end if Hermes is undetected",c("Read start (HR2 protocol)"="start","Read end (reverse-complemented read)"="end","Require Hermes orientation evidence"="require")),
      numericInput("ms_error","Allowed Hermes mismatch fraction",.15,min=0,max=.3,step=.01),
      textInput("ms_hermes","Complete Hermes / construct FASTA (blank uses bundled reference)","")),
    actionButton("ms_run","Identify mutants",class="btn-primary"),actionButton("ms_cancel","Cancel"),uiOutput("ms_status"),
    tags$details(class="ytab-more-options",tags$summary("Current saved identification"),uiOutput("ms_saved_control")),
    tags$details(class="ytab-more-options",tags$summary("Downloads"),downloadButton("ms_download_table","Identification CSV"),
      downloadButton("ms_download_json","Results and provenance JSON"),downloadButton("ms_download_flanks","Trimmed genomic FASTA")))
}

mutant_screen_ui <- function(){
  tagList(panel_card("Hermes insertion-site identification by nested arbitrary PCR",
    p("Two rounds of nested arbitrary PCR recover genomic DNA adjacent to a Hermes insertion, followed by Sanger sequencing with HR-2."),
    tags$ol(tags$li("PCR1: genomic DNA or cell lysate is amplified with Hermes-specific HR-1 and degenerate ARB6. Low-stringency priming establishes a genomic flank."),
      tags$li("PCR2: HR-2 and ARB2 enrich the junction-containing product from PCR1 using OneTaq conditions."),
      tags$li("Assess PCR2 products on an agarose gel, purify suitable products, and sequence with HR-2.")),
    tags$figure(tags$img(src="ytab-mutant-method/gel_image.png",alt="Supplied PCR product gel: 100 bp ladder, M1–M9 and NC",style="max-width:100%;width:760px;height:auto;"),
      tags$figcaption("Supplied gel: 100 bp ladder, M1–M9 and negative control (NC). Gel lanes are not automatically matched to sequencing filenames.")),
    tags$details(class="ytab-technical-details",tags$summary("Primer sequences and full protocol"),
      tags$table(class="table table-condensed",tags$thead(tags$tr(tags$th("Primer"),tags$th("Sequence (5′–3′)"))),tags$tbody(
        tags$tr(tags$td("HR-1"),tags$td("AATACTTGCACTCAAAAGGCTTGACACC")),
        tags$tr(tags$td("ARB6"),tags$td("GGCACGCGTCGACTAGTACNNNNNNNNNNACGCC")),
        tags$tr(tags$td("HR-2"),tags$td("ACACTTGTGCTTATCTATGTGGC")),
        tags$tr(tags$td("ARB2"),tags$td("GGCACGCGTCGACTAGTAC")))),
      tags$a(href="ytab-mutant-method/mutant_screen_method.pdf",target="_blank","Open the supplied protocol PDF")),
    tags$details(class="ytab-technical-details",tags$summary("References"),
      p(tags$a(href="https://doi.org/10.1016/j.mimet.2005.02.011","Das et al. (2005). An improved arbitrary primed PCR method.")),
      p(tags$a(href="https://doi.org/10.3389/fbioe.2014.00046","Martínez-García et al. (2014). New transposon tools.")),
      p(tags$a(href="https://doi.org/10.1186/s44330-025-00054-3","Alejaldre et al. (2026). Automated workflow for genotyping individual transposon library variants.")),
      p(tags$a(href="https://www.ncbi.nlm.nih.gov/nuccore/L34807.1","Hermes reference: GenBank L34807.1")))),
    panel_card("Mutant identifications",p("Select input and click Identify mutants. Genomic flanks are aligned to the selected reference genome; multiple-CDS overlaps are resolved using the junction-facing read orientation. Identification is universal unless a project is explicitly selected for screen links."),
      p(class="help-block","One row per mutant. Matching .seq/.ab1 exports are combined using ABI quality scores. Multiple candidate genes are listed together; inspect the mutant below for their context. Aligned bases are query bases aligned after quality trimming and adapter clipping, not the raw read length."),
      DT::DTOutput("ms_results"),
      tags$details(id="ms_inspect_mutant",class="ytab-technical-details",
        ontoggle="if (this.open) window.dispatchEvent(new Event('resize'));",
        tags$summary("Inspect mutant"),uiOutput("ms_mutant_control"),uiOutput("ms_details"))),
    tags$details(id="ms_gene_information",class="ytab-card ytab-technical-details",
      ontoggle="if (this.open) window.dispatchEvent(new Event('resize'));",
      tags$summary("Gene and optional screen information"),DT::DTOutput("ms_genes"),uiOutput("ms_link_context"),tags$h5("Fitness summary"),DT::DTOutput("ms_screen_results"),
      tags$details(class="ytab-technical-details",tags$summary("Per-sample insertion support"),DT::DTOutput("ms_sample_support")),uiOutput("ms_explorer_control")),
    tags$details(class="ytab-technical-details",tags$summary("Run provenance and limitations"),verbatimTextOutput("ms_provenance")))
}

ms_python_command <- function(python_bin,repo_root){
  candidates<-unique(c(python_bin(),file.path(R.home(),"../..","bin","python"),file.path(repo_root,".venv/bin/python")))
  for(candidate in candidates){
    if(!file.exists(candidate))next
    check<-tryCatch(processx::run(candidate,c("-c","import Bio, pysam, yaml, cutadapt"),error_on_status=FALSE),error=function(e)NULL)
    if(!is.null(check)&&check$status==0L)return(normalizePath(candidate))
  }
  stop("Select a Python environment containing Biopython, pysam, PyYAML and cutadapt.",call.=FALSE)
}

ms_rows <- function(rows){
  if(!length(rows))return(data.frame())
  keys<-unique(unlist(lapply(rows,names)))
  do.call(rbind,lapply(rows,function(r)as.data.frame(setNames(lapply(keys,function(k){v<-r[[k]];if(is.null(v)||!length(v))"" else paste(as.character(v),collapse="; ")}),keys),stringsAsFactors=FALSE,check.names=FALSE)))
}

ms_natural_order <- function(values){
  keys<-tolower(as.character(values));matches<-gregexpr("[0-9]+",keys,perl=TRUE)
  regmatches(keys,matches)<-lapply(regmatches(keys,matches),function(x)sprintf("%020.0f",as.numeric(x)))
  order(keys,seq_along(keys))
}

ms_reference_label <- function(project){
  filename<-basename(project$reference$fasta%||%"")
  accession<-regmatches(filename,regexpr("GC[AF]_[0-9]+(\\.[0-9]+)?",filename,perl=TRUE))
  short<-if(length(accession)&&nzchar(accession))accession else if(nchar(filename)>24L)paste0(substr(filename,1,21),"…")else filename
  paste(project$species%||%"Reference",short,sep=" — ")
}

# Keep native IDs in saved results and explorer queries; reuse Fitness display annotations.
ms_gene_display <- function(ids,repo_root){
  x<-ytab_glabrata_regular_display_columns(data.frame(gene_id=ids),repo_root,"gene_id",include_details=TRUE)
  x[,setdiff(names(x),"S. cerevisiae orthology"),drop=FALSE]
}

ms_identification_table <- function(reads,repo_root,project){
  reads<-reads[ms_natural_order(vapply(reads,function(r)r$mutant,""))]
  rows<-lapply(reads,function(r){
    ids<-if(length(r$candidates))vapply(r$candidates,function(g)g$gene_id,"")else""
    display<-ms_gene_display(ids,repo_root)
    row<-list(Mutant=r$mutant,`Gene ID`=paste(display[["Gene ID"]],collapse="; "),
      `Gene name`=paste(display[["Gene name"]],collapse="; "),
      Chromosome=qc_plot_chromosome_display(r$chromosome%||%"",project),
      `Aligned bases after QC (bp)`=r$aligned_bases%||%NA_real_,
      `Hermes overlap (bp)`=r$hermes_overlap_before_qc_bases%||%NA_real_,
      `Identity (%)`=if(is.null(r$identity))NA_real_ else round(100*r$identity,2))
    if(!length(r$candidates)){
      row[["Gene name"]]<-if(identical(r$gene_assignment$rule,"no_cds_overlap"))"Intergenic" else
        if(identical(r$gene_assignment$rule,"noncoding_gene_overlap"))"Noncoding gene region" else
        if(identical(r$status,"No alignment passes thresholds"))"No passing alignment" else
        if(grepl("Ambiguous mapping",r$status,fixed=TRUE))"Ambiguous mapping" else"Unassigned"
    }
    for(col in setdiff(names(display),c("Gene ID","Gene name"))){
      row[[col]]<-paste(vapply(seq_len(nrow(display)),function(i){
        value<-display[[col]][i]
        if(nrow(display)>1L)paste0(display[["Gene ID"]][i],": ",value)else value
      },""),collapse="\n\n")
    }
    row
  })
  x<-ms_rows(rows)
  for(col in intersect(c("Aligned bases after QC (bp)","Hermes overlap (bp)","Identity (%)"),names(x))){
    values<-x[[col]];values[values%in%c("NA","")]<-NA_character_
    x[[col]]<-as.numeric(values)
  }
  x
}

ms_identification_widget <- function(data){
  if(!nrow(data))return(ytab_datatable(data))
  # Sort by an explicit numeric rank: browser natural/date parsing of labels such
  # as 1-HR2 is unreliable. Keep displayed mutant labels unchanged.
  data$.ms_mutant_order<-match(seq_len(nrow(data)),ms_natural_order(data$Mutant))
  rank_col<-match(".ms_mutant_order",names(data))-1L
  ytab_gene_details_datatable(data,options=list(order=list(list(0,"asc")),
    columnDefs=list(list(orderData=rank_col,targets=0),
      list(visible=FALSE,searchable=FALSE,targets=rank_col))))
}

ms_candidate_table <- function(r,repo_root,project){
  x<-ms_rows(r$candidates);if(!nrow(x))return(x)
  display<-ms_gene_display(x$gene_id,repo_root)
  location<-vapply(r$candidates,function(g){
    if(!is.null(g$gene_start)&&!is.null(g$gene_end))paste0(format(g$gene_start,big.mark=","),"–",format(g$gene_end,big.mark=","))
    else if(!is.null(g$distance_bp))paste(format(g$distance_bp,big.mark=","),"bp from mapped flank")else""
  },"")
  cbind(display[c("Gene ID","Gene name")],data.frame(
    Chromosome=rep(qc_plot_chromosome_display(r$chromosome%||%"",project),nrow(x)),
    Location=location,Interpretation=gsub("; (exact )?junction unconfirmed","",x$context),check.names=FALSE),
    display[setdiff(names(display),c("Gene ID","Gene name"))])
}

ms_linked_table <- function(r,repo_root,kind="Fitness screen"){
  rows<-Filter(function(x)identical(x$kind,kind),r$screen_results)
  if(!length(rows))return(data.frame())
  x<-ms_rows(rows);display<-ms_gene_display(x$gene_id,repo_root)
  metrics<-if(identical(kind,"Fitness screen"))data.frame(
    Analysis=x$analysis_or_sample,
    `Fitness call`=vapply(rows,function(z)z$final_call%||%z$fitness_call%||%"Not available",""),
    `Mean log2 fold change`=vapply(rows,function(z)as.character(z$mean_log2FC%||%z$mean_log2fc%||%"Not available"),""),check.names=FALSE)
  else data.frame(Sample=x$analysis_or_sample,
    `Insertion hits`=vapply(rows,function(z)as.character(z$hits%||%"Not available"),""),
    Reads=vapply(rows,function(z)as.character(z$reads%||%"Not available"),""),check.names=FALSE)
  cbind(display[c("Gene ID","Gene name")],metrics,display[setdiff(names(display),c("Gene ID","Gene name"))])
}

mutant_screen_server <- function(input,output,session,repo_root,active,active_project_path,python_bin,go_to,
                                storage_root=file.path(repo_root,"output/mutant_screen")){
  assets<-file.path(repo_root,"data/mutant_screen")
  if(dir.exists(assets))addResourcePath("ytab-mutant-method",assets)
  results<-reactiveVal(NULL);status<-reactiveVal("Ready to identify mutants.");tick<-reactiveVal(0L)
  job<-new.env(parent=emptyenv());job$process<-NULL;job$project<-NULL;job$output<-NULL;job$work<-NULL
  stop_worker<-function(){if(!is.null(job$process)&&job$process$is_alive())job$process$kill_tree();job$process<-NULL}
  session$onSessionEnded(stop_worker)
  projects<-reactive({discover_ytab_projects(repo_root)})
  output$ms_reference_control<-renderUI({available<-projects();if(!nrow(available))return(NULL)
    refs<-lapply(available$project_config,function(path)read_project_summary(path,repo_root))
    keys<-vapply(refs,function(p)paste(p$species,p$reference$fasta%||%""),"")
    refs<-refs[!duplicated(keys)]
    choices<-setNames(vapply(refs,function(p)p$project_config,""),make.unique(vapply(refs,ms_reference_label,"")))
    # default: the reference of the open project (same species and FASTA), else the first
    open<-tryCatch(active(),error=function(e)NULL);key<-if(is.null(open))"" else paste(open$species,open$reference$fasta%||%"")
    selected<-if(key%in%keys[!duplicated(keys)])choices[[match(key,keys[!duplicated(keys)])]] else choices[[1]]
    selectInput("ms_reference_config","Alignment reference genome",choices,selected=isolate(input$ms_reference_config)%||%selected)})
  reference_config<-reactive({input$ms_reference_config%||%active_project_path()})
  display_reference<-reactive({d<-results()
    if(!is.null(d$reference_species))return(list(species=d$reference_species,repo_root=repo_root))
    path<-d$project_config%||%reference_config()
    read_project_summary(path,repo_root)})
  output$ms_link_project_control<-renderUI({p<-projects()
    choices<-if(nrow(p))setNames(p$project_config,p$display_name)else character()
    selectInput("ms_link_project","Link screen results to project (optional)",c("None — universal identification"="",choices),selected="")})
  output$ms_fitness_control<-renderUI({path<-input$ms_link_project%||%"";if(!nzchar(path))return(NULL)
    p<-read_project_summary(path,repo_root);runs<-discover_fitness_results(p)
    choices<-c("All available gene-level Fitness analyses"="",setNames(vapply(runs,function(r)r$table,""),vapply(runs,function(r)paste(r$analysis_id,if(r$legacy)"(legacy result)"else"(completed result)"),"")))
    selectInput("ms_fitness_result","Screen result to link",choices)})
  output$ms_saved_control<-renderUI({tick();file<-file.path(storage_root,"current/results.json")
    choices<-if(file.exists(file))c("Current identification"=file)else character()
    selectInput("ms_saved_run","Saved identification",c("Load current result…"="",choices))})
  observeEvent(TRUE,{file<-file.path(storage_root,"current/results.json")
    if(file.exists(file))tryCatch({results(jsonlite::fromJSON(file,simplifyVector=FALSE));status("Current identification loaded.")},error=function(e)status("Current result could not be loaded."))},once=TRUE)
  observeEvent(input$ms_saved_run,{if(!nzchar(input$ms_saved_run%||%""))return()
    tryCatch({
      file<-normalizePath(input$ms_saved_run,mustWork=TRUE)
      expected<-normalizePath(file.path(storage_root,"current/results.json"),mustWork=FALSE)
      if(!identical(file,expected))stop("Select the current universal identification.")
      d<-jsonlite::fromJSON(file,simplifyVector=FALSE)
      if(!identical(d$project_id,"universal")||!is.list(d$results))stop("Saved result is not a universal identification.")
      results(d);status("Saved identification run loaded.")
    },error=function(e)status(paste("Could not load saved run:",conditionMessage(e))))},ignoreInit=TRUE)
  observeEvent(input$ms_run,{
    if(!is.null(job$process)&&job$process$is_alive())return(showNotification("Identification is already running."))
    tryCatch({
      p<-active();req(p);work<-tempfile("ytab-mutant-input-");dir.create(work)
      if(identical(input$ms_input_mode,"upload")){
        uploaded<-input$ms_files;validate(need(!is.null(uploaded)&&nrow(uploaded)>0,"Upload sequence files or an archive."))
        folders<-file.path(work,as.character(seq_len(nrow(uploaded))))
        invisible(lapply(folders,dir.create))
        paths<-file.path(folders,basename(uploaded$name))
        if(!all(file.copy(uploaded$datapath,paths)))stop("Could not copy uploaded inputs.")
      }else{
        paths<-trimws(strsplit(input$ms_paths%||%"","\n",fixed=TRUE)[[1]]);paths<-paths[nzchar(paths)]
        if(!length(paths))stop("Select input files or a directory.")
        paths<-vapply(paths,function(x)if(grepl("^/|^~",x))path.expand(x)else file.path(repo_root,x),"")
      }
      runner<-ms_python_command(python_bin,repo_root)
      out<-file.path(storage_root,"current")
      args<-c(file.path(repo_root,"scripts/local/ytab_identify_mutants.py"),"--project-config",reference_config(),"--input",paths,"--output",out,"--replace-current",
        "--min-length",as.character(input$ms_min_length%||%80),"--min-identity",as.character((input$ms_identity%||%90)/100),
        "--min-mapq",as.character(input$ms_mapq%||%20),"--quality",as.character(input$ms_quality%||%20),
        "--overlap",as.character(input$ms_overlap%||%12),"--junction-min-overlap",as.character(input$ms_junction_overlap%||%20),"--error-rate",as.character(input$ms_error%||%.15),
        "--junction-gene-window",as.character(input$ms_gene_window%||%100),"--junction-gene-min-bases",as.character(input$ms_gene_min_bases%||%50),"--junction-read-end",input$ms_junction_end%||%"start",
        "--condition",input$ms_condition%||%"","--pool",input$ms_pool%||%"")
      if(nzchar(trimws(input$ms_hermes%||%""))){h<-trimws(input$ms_hermes);if(!grepl("^/|^~",h))h<-file.path(repo_root,h);args<-c(args,"--hermes",path.expand(h))}
      if(nzchar(input$ms_link_project%||%"")){
        args<-c(args,"--link-project-config",input$ms_link_project)
        if(nzchar(input$ms_fitness_result%||%""))args<-c(args,"--fitness-result",input$ms_fitness_result)
      }
      job$project<-active_project_path();job$output<-out;job$work<-work
      job$process<-processx::process$new(runner,args,wd=repo_root,stdout=file.path(work,"stdout.log"),stderr=file.path(work,"stderr.log"),cleanup=TRUE,cleanup_tree=TRUE)
      status("Clipping Hermes sequence and aligning genomic flanks…");tick(tick()+1L)
    },error=function(e){status(paste("Identification could not start:",conditionMessage(e)));showNotification(conditionMessage(e),type="error")})
  })
  observeEvent(input$ms_cancel,{stop_worker();status("Identification cancelled; previous results retained.")})
  observe({tick();if(is.null(job$process))return();invalidateLater(500,session);if(job$process$is_alive())return()
    exit<-job$process$get_exit_status();job$process<-NULL
    info<-tryCatch(jsonlite::fromJSON(tail(readLines(file.path(job$work,"stdout.log"),warn=FALSE),1)),error=function(e)NULL)
    if(exit==0L&&!is.null(info$output_dir)){
      returned<-normalizePath(info$output_dir,mustWork=FALSE)
      expected<-normalizePath(file.path(storage_root,"current"),mustWork=FALSE)
      if(identical(returned,expected))job$output<-returned
    }
    file<-file.path(job$output,"results.json")
    if(exit==0L&&file.exists(file)){results(jsonlite::fromJSON(file,simplifyVector=FALSE));status(if(isTRUE(info$reused))"Previous matching identification reused."else"Identification complete.");tick(isolate(tick())+1L)}
    else{err<-file.path(job$work,"stderr.log");status(paste("Identification failed:",if(file.exists(err))paste(tail(readLines(err,warn=FALSE),15),collapse="\n")else"Worker returned no result."))}})
  output$ms_status<-renderUI(p(status()))
  reads<-reactive({d<-results();if(is.null(d))return(list());r<-d$results
    # Older saved runs remain immutable; clarify sufficient-overlap labels at display time.
    threshold<-d$parameters$junction_min_overlap%||%20
    r<-lapply(r,function(x){
      if(identical(x$status,"Genomic flank mapped; junction unconfirmed")&&
         (x$hermes_overlap_bases%||%0)>=threshold){
        x$status<-if(!is.null(x$junction_adjacent_base))"Hermes junction supported"else"Hermes overlap supported; exact junction position unresolved"
        x$junction_review_reason<-gsub("leaves the junction unconfirmed","prevents an exact junction-position call",x$junction_review_reason%||%"")
        x$candidates<-lapply(x$candidates,function(g){g$context<-gsub("junction unconfirmed","exact junction position unresolved",g$context%||%"");g})
      };x
    })
    r[ms_natural_order(vapply(r,function(x)x$mutant,""))]})
  output$ms_results<-DT::renderDT({ms_identification_widget(ms_identification_table(reads(),repo_root,display_reference()))})
  output$ms_mutant_control<-renderUI({r<-reads();if(!length(r))return(NULL)
    selectInput("ms_selected","Inspect mutant",setNames(vapply(r,function(x)x$read_id,""),vapply(r,function(x)paste(x$mutant,"—",x$source),"")))})
  selected<-reactive({r<-reads();if(!length(r))return(NULL);key<-input$ms_selected%||%r[[1]]$read_id;hit<-Filter(function(x)x$read_id==key,r);if(length(hit))hit[[1]]else NULL})
  output$ms_details<-renderUI({r<-selected();if(is.null(r))return(NULL)
    tagList(if(grepl("No alignment|Ambiguous mapping|Insufficient genomic",r$status))p(class="text-warning",r$status),
      p(paste0("Chromosome: ",qc_plot_chromosome_display(r$chromosome%||%"Not mapped",display_reference()),
      "; identity: ",if(is.null(r$identity))"Unavailable"else paste0(round(100*r$identity,2),"%"),
      "; MAPQ: ",r$mapq%||%"Unavailable","; Hermes overlap before QC/clipping: ",r$hermes_overlap_before_qc_bases%||%"Not stored"," bp",
      "; Hermes overlap after QC: ",r$hermes_overlap_bases%||%"Not stored"," bp",
      "; trimmed bases: ",r$clipped_bases,"; genomic read length: ",r$clean_length," bp.")),
      p(if(isTRUE(r$quality_available))"Per-base quality scores were available."else"No per-base quality scores in this file; review the Sanger trace for mixed or low-quality sequence."),
      p(class=if(identical(r$sequence_export_agrees,FALSE))"text-warning"else"help-block",r$input_merge%||%"Single input in this older run."),
      if(!is.null(r$gene_assignment))p(class="help-block",if(identical(r$gene_assignment$rule,"multiple_cds_junction_proximal"))
        paste(if(length(r$candidates))"Multiple-CDS overlap resolved from the junction-facing end."else"Multiple-CDS overlap remains unresolved.",
          r$gene_assignment$orientation_basis%||%"",paste0(r$gene_assignment$aligned_bases%||%0," aligned bases in the proximal window."))
        else sub("; exact junction unconfirmed","",r$gene_assignment$reason,fixed=TRUE)),
      tags$details(tags$summary("Clipping evidence and alternate alignments"),tags$pre(jsonlite::toJSON(list(
        junction_status=r$status,junction_review=r$junction_review_reason,junction_adjacent_base=r$junction_adjacent_base,
        gene_assignment=r$gene_assignment,hermes_before_qc=r$hermes_clipping_before_qc,clipping_after_qc=r$clipping,alignments=r$alignments),auto_unbox=TRUE,pretty=TRUE))))})
  output$ms_genes<-DT::renderDT({r<-selected();if(is.null(r))return(ytab_datatable(data.frame()))
    ytab_gene_details_datatable(ms_candidate_table(r,repo_root,display_reference()))})
  output$ms_link_context<-renderUI({d<-results();p(class="help-block",if(is.null(d$linked_project_id))
    "Universal identification: no project screen results linked. Select a project explicitly to add screen associations."
    else paste("Screen associations explicitly linked to:",d$linked_project_id))})
  output$ms_screen_results<-DT::renderDT({r<-selected();if(is.null(r))return(ytab_datatable(data.frame()))
    ytab_gene_details_datatable(ms_linked_table(r,repo_root))})
  output$ms_sample_support<-DT::renderDT({r<-selected();if(is.null(r))return(ytab_datatable(data.frame()))
    ytab_gene_details_datatable(ms_linked_table(r,repo_root,"Raw sample summary"))})
  output$ms_explorer_control<-renderUI({r<-selected();if(is.null(r)||!length(r$candidates))return(NULL)
    tagList(p("Screen associations are gene-level results; nearby intergenic genes and uncertain junctions are not confirmed disrupted genes."),
      selectInput("ms_explorer_gene","Gene to inspect",{
        ids<-vapply(r$candidates,function(g)g$gene_id,"");d<-ms_gene_display(ids,repo_root)
        setNames(ids,paste(d[["Gene ID"]],d[["Gene name"]],sep=" — "))
      }),if(!is.null(results()$linked_project_id)&&identical(results()$linked_project_id,active()$project_id))
        actionButton("ms_open_explorer","View gene & domain insertions")
      else p(class="help-block","Insertion tracks require an explicitly linked project; open that project to inspect its tracks."))})
  observeEvent(input$ms_open_explorer,{req(input$ms_explorer_gene,results()$linked_project_id)
    req(identical(results()$linked_project_id,active()$project_id))
    go_to("gene_explorer","gene_explorer_tabs","gene_domain");updateTextInput(session,"gene_domain_query",value=input$ms_explorer_gene)},ignoreInit=TRUE)
  output$ms_provenance<-renderText({d<-results();if(is.null(d))return("No identification run selected.")
    jsonlite::toJSON(d[c("created","project_id","linked_project_id","reference","reference_sha256","annotation_metadata","hermes_reference","hermes_sha256","parameters","versions","limitations")],auto_unbox=TRUE,pretty=TRUE)})
  download<-function(id,name){output[[id]]<-downloadHandler(filename=function()name,content=function(file){d<-results();req(d);path<-file.path(d$output_dir,name);req(file.exists(path));file.copy(path,file)})}
  download("ms_download_table","mutant_identifications.csv");download("ms_download_json","results.json");download("ms_download_flanks","genomic_flanks.fasta")
}
