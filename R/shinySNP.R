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
                     chr_pos, top_snps_tbl,
                     snp_action = reactive({"basic"})) {
  ns <- session$ns
  
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
  output$snp_sum <- renderUI({
    switch(input$snp_sum,
           best   = dataTableOutput(ns("top_snps_best")),
           indels = dataTableOutput(ns("top_indels")),
           peaks  = dataTableOutput(ns("top_snps_peak")),
           range  = dataTableOutput(ns("top_snps_tbl")))
  })
  
  ## Downloads.
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("top_snps_", chr_pos(), "_", snp_action(), ".csv")) },
    content = function(file) {
      write.csv(best_http(), file)
    }
  )
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinySNP
#' @export
shinySNPInput <- function(id) {
  ns <- NS(id)
  tagList(
    selectInput(ns("snp_sum"), NULL, c("best","indels","peaks","range"))
  )
}
#' @rdname shinySNP
#' @export
shinySNPUI <- function(id) {
  ns <- NS(id)
  downloadButton(ns("downloadData"), "CSV")
}
#' @rdname shinySNP
#' @export
shinySNPOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_sum"))
}
