sample_display_value <- function(data, columns, fallback = "") {
  column <- intersect(columns, names(data))
  if (length(column)) as.character(data[[column[[1L]]]]) else rep(fallback, nrow(data))
}

sample_display_included <- function(data) {
  if (!"include" %in% names(data)) return(rep("Yes", nrow(data)))
  ifelse(tolower(as.character(data$include)) %in% c("true", "1", "yes"), "Yes", "No")
}

compact_sample_table <- function(samples, status = NULL, selected = NULL,
                                 include_status = TRUE,
                                 include_selected = FALSE) {
  if (!is.data.frame(samples) || !nrow(samples)) return(data.frame())
  out <- data.frame(
    Sample = as.character(samples$sample),
    Role = sample_display_value(samples, c("guessed_condition", "condition", "treatment")),
    Background = sample_display_value(samples, c("guessed_background", "background")),
    Pool = sample_display_value(samples, c("guessed_pool", "pool")),
    Layout = sample_display_value(samples, c("layout")),
    Included = sample_display_included(samples),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (isTRUE(include_selected))
    out <- cbind(Selected = ifelse(out$Sample %in% selected, "Yes", "No"),
                 out, stringsAsFactors = FALSE)
  if (isTRUE(include_status) && is.data.frame(status) && nrow(status)) {
    index <- match(out$Sample, as.character(status$sample))
    status_value <- function(column) {
      if (!column %in% names(status)) return(rep("", nrow(out)))
      value <- as.character(status[[column]][index])
      value[is.na(value)] <- ""
      value
    }
    out$Mapping <- status_value("mapping")
    out$`Hit file` <- status_value("hit_file")
    out$`Summary table` <- status_value("summary")
  }
  out
}

sample_selector_ui <- function(id,title="Select samples",show_presets=TRUE,show_status=TRUE) {
  ns<-NS(id);tagList(tags$h3(title),selectizeInput(ns("samples"),"Stage selection",choices=character(),multiple=TRUE),
    if(show_presets)tags$div(class="ytab-actions",actionButton(ns("all"),"Select all included"),actionButton(ns("none"),"Select none"),uiOutput(ns("condition_buttons")),actionButton(ns("incomplete"),"Select incomplete"),actionButton(ns("failed"),"Select failed")),
    textOutput(ns("count")),if(show_status)DT::DTOutput(ns("table")))
}

sample_selector_server <- function(id,sample_data,default_samples=NULL,status_data=NULL,selected_state=NULL) {
  moduleServer(id,function(input,output,session){selected<-selected_state%||%reactiveVal(character());initialized<-reactiveVal(FALSE)
    included_names<-reactive({d<-sample_data();if(!nrow(d))return(character());inc<-if("include"%in%names(d))tolower(as.character(d$include))%in%c("true","1","yes")else rep(TRUE,nrow(d));as.character(d$sample[inc])})
    observeEvent(sample_data(),{d<-sample_data();choices<-if(nrow(d))as.character(d$sample)else character();defaults<-default_samples;if(is.null(defaults))defaults<-included_names();if(!initialized()&&!length(selected()))selected(defaults);initialized(TRUE);selected(intersect(selected(),choices));updateSelectizeInput(session,"samples",choices=choices,selected=selected(),server=TRUE)},ignoreInit=FALSE)
    observeEvent(selected(),updateSelectizeInput(session,"samples",selected=selected(),server=TRUE),ignoreInit=TRUE)
    observeEvent(input$samples,selected(input$samples %||% character()),ignoreInit=TRUE)
    observeEvent(input$all,{selected(included_names());updateSelectizeInput(session,"samples",selected=selected(),server=TRUE)},ignoreInit=TRUE)
    observeEvent(input$none,{selected(character());updateSelectizeInput(session,"samples",selected=character(),server=TRUE)},ignoreInit=TRUE)
    condition_col<-reactive({d<-sample_data();if("guessed_condition"%in%names(d))"guessed_condition"else if("condition"%in%names(d))"condition"else NULL})
    choose_condition<-function(value){d<-sample_data();col<-condition_col();req(col);wanted<-intersect(as.character(d$sample[tolower(as.character(d[[col]]))==value]),included_names());selected(wanted);updateSelectizeInput(session,"samples",selected=wanted,server=TRUE)}
    observeEvent(input$parents,choose_condition("parent"),ignoreInit=TRUE);observeEvent(input$treated,choose_condition("treated"),ignoreInit=TRUE)
    observeEvent(input$incomplete,{s<-status_data();wanted<-s$sample[s$included&(!s$mapping%in%c("success","skipped")|!s$hit_file%in%c("success","skipped")|!s$summary%in%c("success","skipped"))];selected(wanted);updateSelectizeInput(session,"samples",selected=wanted,server=TRUE)},ignoreInit=TRUE)
    observeEvent(input$failed,{s<-status_data();wanted<-s$sample[s$included&(s$mapping=="failed"|s$hit_file=="failed"|s$summary=="failed")];selected(wanted);updateSelectizeInput(session,"samples",selected=wanted,server=TRUE)},ignoreInit=TRUE)
    output$condition_buttons<-renderUI({d<-sample_data();col<-condition_col();if(is.null(col))return(NULL);values<-tolower(as.character(d[[col]]));tagList(if("parent"%in%values)actionButton(session$ns("parents"),"Select parents"),if("treated"%in%values)actionButton(session$ns("treated"),"Select treated"))})
    output$count<-renderText(sprintf("Selected samples: %d",length(selected())))
    output$table<-DT::renderDT({d<-sample_data();s<-status_data();if(!nrow(d))return(DT::datatable(data.frame(),selection="none"));out<-compact_sample_table(d,s,selected(),include_status=TRUE,include_selected=TRUE);DT::datatable(out,rownames=FALSE,filter="top",selection="none",options=list(scrollX=TRUE,pageLength=12))})
    selected
  })
}
