#' Shiny setup for DOQTL selection
#'
#' Shiny module for phenotype selection.
#'
#' @param input,output,session standard shiny arguments
#' @param pheno_type,peaks_tbl,analyses_tbl,hot_peak,win_par reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinySetup <- function(input, output, session,
                        pheno_type, peaks_tbl, pmap_obj, analyses_tbl) {
  ns <- session$ns
  
  ## Peak counts.
  hot_peak <- callModule(shinyPeaks, "shinypeaks",
                         pheno_type, peaks_tbl, pmap_obj)
  
  ## Use peaks as input to shinyWindow.
  win_par <- callModule(shinyWindow, "window",
                        pmap_obj, hot_peak)
  
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
  observeEvent(pheno_anal(), {
    output$num_pheno <- renderText({
      pheno <- req(pheno_anal())
      num_pheno(pheno, analyses_tbl())
    })
  })
  
  ## Use window as input to shinyPhenos.
  pheno_anal <- callModule(shinyPhenos, "phenos",
                           pheno_type, peaks_tbl, analyses_tbl,
                           hot_peak, win_par)
  
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
  
  output$setup <- renderUI({
    switch(input$setup,
           "Peak Count" = shinyPeaksOutput(ns("shinypeaks")),
           "Peaks"      = dataTableOutput(ns("peaks_tbl")),
           "Raw Data"   = shinyPhenoPlotUI(ns("PhenoPlotRaw")),
           "Trans Data" = shinyPhenoPlotUI(ns("PhenoPlotTrans")),
           "Covariates" = dataTableOutput(ns("analyses_tbl")))
  })
  
  ## Return.
  reactive({
    list(pheno_anal = pheno_anal(),
         win_par = list(chr_id     = win_par$chr_id,
                        peak_Mbp   = win_par$peak_Mbp,
                        window_Mbp = win_par$window_Mbp))
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
      shinyPhenosUI(ns("phenos")),
      shinyWindowUI(ns("window")),
      shinyPeaksInput(ns("shinypeaks")),
      radioButtons(ns("setup"), NULL,
                   c("Peak Count","Peaks","Raw Data","Trans Data","Covariates")))),
    mainPanel(uiOutput(ns("setup")))
  )
}
