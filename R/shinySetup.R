#' Shiny setup for DOQTL selection
#'
#' Shiny module for phenotype selection.
#'
#' @param input,output,session standard shiny arguments
#' @param pheno_type,peaks_tbl,pmap_obj,analyses_tbl,cov_mx reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinySetup <- function(input, output, session,
                        pheno_type, peaks_tbl, pmap_obj, analyses_tbl, cov_mx) {
  ns <- session$ns
  
  # Select phenotype dataset
  output$dataset <- renderUI({
    req(choices <- pheno_type())
    if(is.null(selected <- input$dataset))
      selected <- NULL
    selectInput(ns("dataset"), "Phenotype Group",
                choices = as.list(choices),
                selected = selected,
                multiple = TRUE)
  })
  
  ## Locate Peak.
  win_par <- callModule(shinyPeaks, "shinypeaks",
                         input, pheno_type, peaks_tbl, pmap_obj)
  
  chr_pos <- reactive({
    make_chr_pos(win_par$chr_id, 
                 win_par$peak_Mbp, win_par$window_Mbp)
  })
  output$chr_pos <- renderText({
    paste0("Region: ", chr_pos(), "Mbp")
  })
  output$num_pheno <- renderText({
    num_pheno(character(), analyses_tbl())
  })
  observeEvent(phe_par$pheno_names, {
    output$num_pheno <- renderText({
      num_pheno(phe_par$pheno_names, analyses_tbl())
    })
  })
  
  ## Use window as input to shinyPhenos.
  phe_par <- callModule(shinyPhenos, "phenos",
             input, peaks_tbl, analyses_tbl, win_par)
  
  ## Density or scatter plot of phenotypes.
  analyses_df <- reactive({
    analyses_tbl() %>%
      filter(pheno %in% req(phe_par$pheno_names))
  })
  phe_df <- reactive({
    get_pheno(pheno_data,
              analyses_df() %>%
                distinct(pheno, .keep_all=TRUE))
  })
  raw_phe_df <- reactive({
    get_pheno(pheno_data,
              analyses_df() %>%
                distinct(pheno, .keep_all=TRUE),
              FALSE)
  })
  callModule(shinyPhenoPlot, "PhenoPlotRaw", raw_phe_df, cov_mx)
  callModule(shinyPhenoPlot, "PhenoPlotTrans", phe_df, cov_mx)
  
  # Output the analyses table
  output$analyses_tbl <- renderDataTable({
    collapse_covar(analyses_df())
  }, options = list(scrollX = TRUE, pageLength = 10))
  
  # Output the peaks table
  
  ## Set up peaks data frame.
  peaks_df <- reactive({
    phename <- req(phe_par$pheno_names) 
    peaks %>%
      filter(pheno %in% phename)
  })
  output$peaks_tbl <- renderDataTable({
    peaks_df() %>%
      select(pheno,chr,pos,lod)
  }, options = list(scrollX = TRUE, pageLength = 10))
  
  ## Setup input logic.
  output$title <- renderUI({
    strong(req(input$setup))
  })
  output$sidebar_setup <- renderUI({
    switch(req(input$setup),
           Phenotypes = tagList(
             shinyPhenosUI(ns("phenos")),
             radioButtons(ns("show_data"), NULL,
                          c("LOD Peaks","Covariates",
                            "Trans Data","Raw Data"),
                          input$show_data)),
           Region = shinyPeaksInput(ns("shinypeaks")))
  })
  output$main_setup <- renderUI({
    switch(req(input$setup),
           Phenotypes = {
             switch(req(input$show_data),
                    "LOD Peaks"  = dataTableOutput(ns("peaks_tbl")),
                    "Raw Data"   = shinyPhenoPlotUI(ns("PhenoPlotRaw")),
                    "Trans Data" = shinyPhenoPlotUI(ns("PhenoPlotTrans")),
                    "Covariates" = dataTableOutput(ns("analyses_tbl")))
             },
           Region = shinyPeaksOutput(ns("shinypeaks")))
  })
  
  output$radio <- renderUI({
    radioButtons(ns("setup"), NULL,
                 c("Region", "Phenotypes"),
                 input$setup,
                 inline=TRUE)
  })
  
  ## Return.
  reactive({
    list(phe_par = phe_par,
         win_par = win_par)
  })
}
#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinySetup
#' @export
shinySetupOutput <- function(id) {
  ns <- NS(id)
  tagList(
    textOutput(ns("num_pheno")),
    uiOutput(ns("chr_pos"))
  )
}
#' @rdname shinySetup
#' @export
shinySetupUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarPanel(
      uiOutput(ns("title")),
      uiOutput(ns("radio")),
      uiOutput(ns("dataset")),
      uiOutput(ns("sidebar_setup"))
    ),
    mainPanel(uiOutput(ns("main_setup")))
  )
}
