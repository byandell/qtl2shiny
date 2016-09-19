#' Shiny Genes and Exons with nearby SNPs module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param chr_pos,top_snps_tbl,gene_exon_tbl reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyGeneExon <- function(input, output, session,
                     chr_pos, top_snps_tbl, gene_exon_tbl) {
  ns <- session$ns
  output$gene_sum <- renderDataTable({
    withProgress(message = 'Gene Exon Table ...', value = 0, {
      setProgress(1)
      summary(gene_exon_tbl())
    })
  }, options = list(scrollX = TRUE, pageLength = 5,
                    lengthMenu = list(c(5,10,20,-1),
                                      list("5","10","20","all"))))
  output$gene_name <- renderUI({
    gene_exon <- req(gene_exon_tbl())
    selectInput(ns("gene_name"), NULL,
                choices = unique(gene_exon$gene),
                selected = input$gene_name)
  })
  plot_gene_exon <- function(gene_name) {
    plot(gene_exon_tbl(), top_snps_tbl(), genes = gene_name)
  }
  output$gene_plot <- renderPlot({
    req(top_snps_tbl(),gene_exon_tbl())
    gene_name <- input$gene_name
    withProgress(message = 'Gene Exon Plot ...', value = 0, {
      setProgress(1)
      plot_gene_exon(gene_name)
    })
  })
  
  ## Outputs
  output$exon_input <- renderUI({
    switch(req(input$button),
           Plot    = uiOutput(ns("gene_name")))
  })
  output$exon_output <- renderUI({
    switch(req(input$button),
           Plot    = plotOutput(ns("gene_plot")),
           Summary = dataTableOutput(ns("gene_sum")))
  })
  
  ## Downloads.
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("gene_exon_", chr_pos(), ".csv")) },
    content = function(file) {
      write.csv(req(gene_exon_tbl()), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("gene_exon_", chr_pos(), ".pdf")) },
    content = function(file) {
      req(top_snps_tbl())
      gene_exon <- req(gene_exon_tbl())
      pdf(file)
      for(gene_name in unique(gene_exon$gene))
        plot_gene_exon(gene_name)
      dev.off()
    }
  )
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyGeneExon
#' @export
shinyGeneExonInput <- function(id) {
  ns <- NS(id)
  fluidRow(
    selectInput(ns("button"), NULL, c("Plot","Summary")),
    uiOutput(ns("exon_input")))
}
#' @rdname shinyGeneExon
#' @export
shinyGeneExonUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(6, downloadButton(ns("downloadData"), "CSV")),
    column(6, downloadButton(ns("downloadPlot"), "Plots")))
}
#' @rdname shinyGeneExon
#' @export
shinyGeneExonOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("exon_output"))
}
