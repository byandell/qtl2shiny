#' Shiny phenotype selection
#'
#' Shiny module for phenotype selection, with interfaces \code{shinyPhenosUI} and  \code{shinyPhenosOutput}.
#'
#' @param input,output,session standard shiny arguments
#' @param set_par,win_par,peaks_df,analyses_tbl,pheno_data,cov_df,project_info reactive arguments
#' @param id shiny identifier
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom dplyr filter
#' @importFrom shiny NS reactive req 
#'   checkboxInput selectInput
#'   uiOutput
#'   renderUI
#'   tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
shinyPhenos <- function(input, output, session,
                        set_par, win_par, peaks_df, analyses_tbl, pheno_data, cov_df,
                        project_info) {
  ns <- session$ns
  # Output the peaks table
  output$peaks <- shiny::renderDataTable({
    dplyr::arrange(
      dplyr::select(
        peaks_df(), pheno, chr, pos, lod),
      dplyr::desc(lod))
  }, options = list(scrollX = TRUE, pageLength = 5,
                    lengthMenu = c(5,10,25)))
  
  ## Density or scatter plot of phenotypes.
  analyses_plot <- shiny::reactive({
    shiny::req(analyses_tbl())
    phename <- shiny::req(set_par$pheno_names)
    dplyr::filter(analyses_tbl(), pheno %in% phename)
  })
  phe_mx <- shiny::reactive({
    pheno_read(pheno_data(), analyses_plot())
  })
  raw_phe_mx <- shiny::reactive({
    pheno_read(pheno_data(), analyses_plot(), FALSE)
  })
  shiny::callModule(shinyPhenoPlot, "PhenoPlotRaw", raw_phe_mx, cov_df)
  shiny::callModule(shinyPhenoPlot, "PhenoPlotTrans", phe_mx, cov_df)
  
  # Show data.
  output$radio <- renderUI({
    shiny::radioButtons(ns("radio"), NULL,
                        c("LOD Peaks","Covariates",
                          "Trans Data","Raw Data"),
                        input$radio)
  })
  output$show_data <- renderUI({
    switch(shiny::req(input$radio),
           "LOD Peaks"  = shiny::dataTableOutput(ns("peaks")),
           "Raw Data"   = shinyPhenoPlotUI(ns("PhenoPlotRaw")),
           "Trans Data" = shinyPhenoPlotUI(ns("PhenoPlotTrans")),
           "Covariates" = shiny::dataTableOutput(ns("analyses_tbl")))
  })
  
  # Output the analyses table
  output$analyses_tbl <- shiny::renderDataTable({
    collapse_covar(analyses_plot())
  }, options = list(scrollX = TRUE, pageLength = 10))
}

shinyPhenosUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("radio"))
}
shinyPhenosOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("show_data"))
}
