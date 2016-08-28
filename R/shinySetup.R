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
  
  ## Window slider
  output$window_Mbp <- renderUI({
    if(is.null(pos <- input$window_Mbp))
      pos <- 3
    sliderInput(ns("window_Mbp"), "Window Half Width",
                0, 6, pos, step=0.5)
  })
  
  ## Peak counts.
  hot_peak <- callModule(shinyPeaks, "shinypeaks",
                         input, pheno_type, peaks_tbl, pmap_obj)

  ## Use peaks as input to shinyChrPeak.
  chr_peak <- callModule(shinyChrPeak, "chr_peak",
                        hot_peak, pmap_obj)
  
  chr_pos <- reactive({
    make_chr_pos(chr_peak$chr_id, 
                 chr_peak$peak_Mbp, input$window_Mbp)
  })
  output$chr_pos <- renderText({
    paste0("Region: ", chr_pos(), "Mbp")
  })
  output$num_pheno <- renderText({
    num_pheno(character(), analyses_tbl())
  })
  observeEvent(pheno_anal(), {
    output$num_pheno <- renderText({
      pheno <- req(pheno_anal())
      num_pheno(pheno, analyses_tbl())
    })
  })
  
  ## Use window as input to shinyPhenos.
  pheno_anal <- callModule(shinyPhenos, "phenos",
                           input, peaks_tbl, analyses_tbl,
                           chr_peak)
  
  ## Density or scatter plot of phenotypes.
  analyses_df <- reactive({
    phenos <- req(pheno_anal())
    analyses_tbl() %>%
      filter(output %in% names(phenos))
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
    phenos <- req(pheno_anal())
    peaks %>%
      filter(output %in% names(phenos))
  })
  output$peaks_tbl <- renderDataTable({
    peaks_df() %>%
      select(pheno,chr,pos,lod)
  }, options = list(scrollX = TRUE, pageLength = 10))
  
  output$sidebar_setup <- renderUI({
    side <- req(input$setup)
    if(side != "Show Data") {
      
    }
    switch(side,
           "Select Data" = {
             tagList(
               uiOutput(ns("dataset")),
               shinyPhenosUI(ns("phenos")),
               shinyChrPeakUI(ns("chr_peak")),
               uiOutput(ns("window_Mbp"))
             )
           },
           "Explore Hotspots" = {
             tagList(
               uiOutput(ns("dataset")),
               shinyPeaksInput(ns("shinypeaks")),
               uiOutput(ns("window_Mbp"))
             )
           },
           "Show Data" = {
             tagList(
               shinyPhenosUI(ns("phenos")),
               radioButtons(ns("show_data"), NULL,
                            c("Raw Data","Trans Data",
                              "LOD Peaks","Covariates"))
             )
           })
  })
  output$main_setup <- renderUI({
    switch(input$setup,
           "Select Data" = {
           },
           "Explore Hotspots" = {
             shinyPeaksOutput(ns("shinypeaks"))
           },
           "Show Data" = {
             switch(input$setup,
                    "LOD Peaks"  = dataTableOutput(ns("peaks_tbl")),
                    "Raw Data"   = shinyPhenoPlotUI(ns("PhenoPlotRaw")),
                    "Trans Data" = shinyPhenoPlotUI(ns("PhenoPlotTrans")),
                    "Covariates" = dataTableOutput(ns("analyses_tbl")))
           })
  })
  
  ## Return.
  reactive({
    list(pheno_anal = pheno_anal(),
         win_par = list(chr_id     = chr_peak$chr_id,
                        peak_Mbp   = chr_peak$peak_Mbp,
                        window_Mbp = input$window_Mbp))
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
    sidebarPanel(tagList(
      h4(strong("Phenotypes")),
      radioButtons(ns("setup"), NULL,
                   c("Select Data",
                     "Explore Hotspots",
                     "Show Data"))),
      uiOutput(ns("sidebar_setup"))
    ),
    mainPanel(uiOutput(ns("main_setup")))
  )
}
