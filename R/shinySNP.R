#' Shiny SNP analysis and plot module
#'
#' Shiny module for scan1 analysis and plots, with interfaces \code{shinySNPInput}, \code{shinySNPUI} and  \code{shinySNPOutput}.
#'
#' @param input,output,session standard shiny arguments
#' @param chr_pos,top_snps_tbl reactive arguments
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
shinySNP <- function(input, output, session,
                     chr_pos, top_snps_tbl,
                     snp_action = shiny::reactive({"basic"})) {
  ns <- session$ns
  
  ensembl_species <- shiny::reactive({"Mus_musculus"})
  best_snps <- shiny::reactive({
    shiny::req(top_snps_tbl())
    summary(top_snps_tbl(),"best")
  })
  best_href <- shiny::reactive({
    best <- shiny::req(best_snps())
    if(!is.null(ensembl_species())) {
      ensembl_URL <- file.path("http://www.ensembl.org",
                               ensembl_species(),
                               "Gene/Summary?g=")
      best$ensembl_gene <-
        unlist(apply(dplyr::select(best, ensembl_gene),
                     1,
                     function(x,y) {
                       if(is.na(x))
                         return("")
                       if(nchar(x))
                         as.character(shiny::a(x, href=paste0(y,x)))
                       else
                         x
                     },
                     ensembl_URL))
    }
    best
  })
  best_http <- shiny::reactive({
    best <- shiny::req(best_snps())
    if(!is.null(ensembl_species())) {
      ensembl_URL <- file.path("http://www.ensembl.org",
                               ensembl_species(),
                               "Gene/Summary?g=")
      best$ensembl_gene <-
        unlist(apply(dplyr::select(best, ensembl_gene),
                     1,
                     function(x,y) {
                       if(is.na(x))
                         return("")
                       if(nchar(x))
                         paste0(y,x)
                       else
                         x
                     },
                     ensembl_URL))
    }
    best
  })

  output$top_snps_tbl <- shiny::renderDataTable({
    shiny::req(top_snps_tbl())
    shiny::withProgress(message = "Top SNP Range ...", value = 0,
    {
      shiny::setProgress(1)
      summary(top_snps_tbl())
    })
  })
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
  })
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

shinySNPInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::selectInput(ns("snp_sum"), NULL, c("best","indels","peaks","range"))
  )
}
shinySNPUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::downloadButton(ns("downloadData"), "CSV")
}
shinySNPOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("snp_sum"))
}
