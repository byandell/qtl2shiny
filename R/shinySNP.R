#' Shiny SNP analysis and plot module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param chr_pos,top_snps_tbl reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinySNP <- function(input, output, session,
                     chr_pos, top_snps_tbl) {
  ensembl_species <- reactive({"Mus_musculus"})
  best_snps <- reactive({
    req(top_snps_tbl())
    summary(top_snps_tbl(),"best")
  })
  best_href <- reactive({
    best <- req(best_snps())
    if(!is.null(ensembl_species())) {
      ensembl_URL <- file.path("http://www.ensembl.org",
                               ensembl_species(),
                               "Gene/Summary?g=")
      best$ensembl_gene <-
        unlist(apply(best %>% select(ensembl_gene),
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
  best_http <- reactive({
    best <- req(best_snps())
    if(!is.null(ensembl_species())) {
      ensembl_URL <- file.path("http://www.ensembl.org",
                               ensembl_species(),
                               "Gene/Summary?g=")
      best$ensembl_gene <-
        unlist(apply(best %>% select(ensembl_gene),
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
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("top_snps_", chr_pos(), ".csv")) },
    content = function(file) {
      write.csv(best_http(), file)
    }
  )

  output$top_snps_tbl <- renderDataTable({
    req(top_snps_tbl())
    withProgress(message = "Top SNP Range ...", value = 0,
    {
      setProgress(1)
      summary(top_snps_tbl())
    })
  })
  output$top_snps_best <- renderDataTable({
    withProgress(message = "Top SNP Best ...", value = 0,
    {
      setProgress(1)
      req(best_href())
  })
  }, escape = FALSE,
    options = list(scrollX = TRUE, pageLength = 10))
  output$top_indels <- renderDataTable({
    withProgress(message = "Top InDels ...", value = 0,
                 {
                   setProgress(1)
                   req(best_href()) %>%
                     filter(type != "SNP")
                 })
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  output$top_snps_peak <- renderDataTable({
    req(top_snps_tbl())
    withProgress(message = "Top SNP Peaks ...", value = 0,
    {
      setProgress(1)
      summary(top_snps_tbl(),"peak")
    })
  })
  output$best <- renderText({"Best SNPs"})
  output$peak <- renderText({"Peak SNPs"})
  output$range <- renderText({"SNP Range"})
}

#' UI for shinySNP Shiny Module
#'
#' UI for scan1 analyses and summaries to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinySNP
#' @export
shinySNPUI <- function(id) {
  ns <- NS(id)
  tagList(fluidRow(
    column(10, tabsetPanel(
      tabPanel("best",
        dataTableOutput(ns("top_snps_best"))),
      tabPanel("indels",
        dataTableOutput(ns("top_indels"))),
      tabPanel("peak",
        dataTableOutput(ns("top_snps_peak"))),
      tabPanel("range",
        dataTableOutput(ns("top_snps_tbl")))
    )),
    column(2, downloadButton(ns("downloadData"), "CSV"))
  ))
}
