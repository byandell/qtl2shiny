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
  chr_peak <- callModule(shinyPeaks, "shinypeaks",
                         input, pheno_type, peaks_tbl, pmap_obj)
  
  chr_pos <- reactive({
    make_chr_pos(chr_peak$chr_id, 
                 chr_peak$peak_Mbp, chr_peak$window_Mbp)
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
  
  ## Setup input logic.
  output$title <- renderUI({
    if(req(input$setup) == "Phenotypes") {
      mytitle <- "Phenotypes"
    } else {
      mytitle <- "Region"
    }
    h4(strong(mytitle))
  })
  output$sidebar_setup <- renderUI({
    if(req(input$setup) == "Phenotypes") {
      tagList(
        shinyPhenosUI(ns("phenos")),
        radioButtons(ns("show_data"), NULL,
                     c("Raw Data","Trans Data",
                       "LOD Peaks","Covariates"))
      )
    } else {
      shinyPeaksInput(ns("shinypeaks"))
    }
  })
  output$main_setup <- renderUI({
    if(req(input$setup) == "Phenotypes") {
      switch(req(input$show_data),
             "LOD Peaks"  = dataTableOutput(ns("peaks_tbl")),
             "Raw Data"   = shinyPhenoPlotUI(ns("PhenoPlotRaw")),
             "Trans Data" = shinyPhenoPlotUI(ns("PhenoPlotTrans")),
             "Covariates" = dataTableOutput(ns("analyses_tbl")))
    } else {
      shinyPeaksOutput(ns("shinypeaks"))
    }
  })
  
  ## Return.
  reactive({
    list(pheno_anal = pheno_anal(),
         win_par = list(chr_id     = chr_peak$chr_id,
                        peak_Mbp   = chr_peak$peak_Mbp,
                        window_Mbp = chr_peak$window_Mbp))
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
      uiOutput(ns("title")),
      radioButtons(ns("setup"), NULL,
                   c("Phenotypes",
                     "Region"),
                   inline=TRUE)),
      uiOutput(ns("dataset")),
      uiOutput(ns("sidebar_setup"))
    ),
    mainPanel(uiOutput(ns("main_setup")))
  )
}
