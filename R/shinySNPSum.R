#' Shiny SNP summary module
#'
#' Shiny module for SNP summary, with interfaces \code{shinySNPInput}, \code{shinySNPUI} and  \code{shinySNPOutput}.
#'
#' @param input,output,session standard shiny arguments
#' @param chr_pos,top_snps_tbl,project_info,snp_action reactive arguments
#' @param id shiny identifier
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom dplyr filter select
#' @importFrom shiny NS reactive req 
#'   selectInput
#'   dataTableOutput uiOutput
#'   renderDataTable renderUI
#'   tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
#' @importFrom ggplot2 autoplot
shinySNPSum <- function(input, output, session,
                     chr_pos, top_snps_tbl, project_info,
                     snp_action = shiny::reactive({"basic"})) {
  ns <- session$ns
  
  best_snps <- shiny::reactive({
    shiny::req(top_snps_tbl())
    summary(top_snps_tbl(),"best")
  })
  best_href <- shiny::reactive({
    best <- shiny::req(best_snps())
    ensembl_gene(best, project_info(), TRUE)
  })
  best_http <- shiny::reactive({
    best <- shiny::req(best_snps())
    ensembl_gene(best, project_info())
  })

  output$top_snps_tbl <- shiny::renderDataTable({
    shiny::req(top_snps_tbl())
    shiny::withProgress(message = "Top SNP Range ...", value = 0,
    {
      shiny::setProgress(1)
      summary(top_snps_tbl())
    })
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  output$top_snps_best <- shiny::renderDataTable({
    shiny::withProgress(message = "Top SNP Best ...", value = 0,
    {
      shiny::setProgress(1)
      shiny::req(best_href())
  })
  }, escape = FALSE,
    options = list(scrollX = TRUE, pageLength = 10))
  output$top_indels <- shiny::renderDataTable({
    shiny::withProgress(message = "Top InDels ...", value = 0,
                 {
                   shiny::setProgress(1)
                   dplyr::filter(shiny::req(best_href()), type != "SNP")
                 })
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  output$top_snps_peak <- shiny::renderDataTable({
    shiny::req(top_snps_tbl())
    shiny::withProgress(message = "Top SNP Peaks ...", value = 0,
    {
      shiny::setProgress(1)
      summary(top_snps_tbl(),"peak")
    })
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  output$snp_sum <- shiny::renderUI({
    switch(input$snp_sum,
           best   = shiny::dataTableOutput(ns("top_snps_best")),
           indels = shiny::dataTableOutput(ns("top_indels")),
           peaks  = shiny::dataTableOutput(ns("top_snps_peak")),
           range  = shiny::dataTableOutput(ns("top_snps_tbl")))
  })
  
  ## Downloads.
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("top_snps_", chr_pos(), "_", snp_action(), ".csv")) },
    content = function(file) {
      write.csv(best_http(), file)
    }
  )
}

shinySNPSumInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::selectInput(ns("snp_sum"), NULL, c("best","indels","peaks","range"))
  )
}
shinySNPSumUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::downloadButton(ns("downloadData"), "CSV")
}
shinySNPSumOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("snp_sum"))
}
