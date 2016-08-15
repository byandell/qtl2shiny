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
      cat(file=stderr(),"gene_sum\n")
      summary(gene_exon_tbl())
    })
  }, options = list(scrollX = TRUE, pageLength = 5,
                    lengthMenu = list(c(5,10,20,-1),
                                      list("5","10","20","all"))))
  output$gene_name <- renderUI({
    gene_exon <- req(gene_exon_tbl())
    cat(file=stderr(),"gene_name\n")
    selectInput(ns("gene_name"), NULL,
                choices = unique(gene_exon$gene))
  })
  plot_gene_exon <- function(gene_name) {
    plot(gene_exon_tbl(), top_snps_tbl(), genes = gene_name)
  }
  output$gene_plot <- renderPlot({
    req(top_snps_tbl(),gene_exon_tbl())
    cat(file=stderr(),"gene_plot\n")
    gene_name <- input$gene_name
    withProgress(message = 'Gene Exon Plot ...', value = 0, {
      setProgress(1)
      plot_gene_exon(gene_name)
    })
  })
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

#' UI for shinyGeneExon Shiny Module
#'
#' UI for scan1 analyses and summaries to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyGeneExon
#' @export
shinyGeneExonUI <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      tabPanel("summary",
               downloadButton(ns("downloadData"), "Download CSV"),
               dataTableOutput(ns("gene_sum"))),
      tabPanel("plots",
               downloadButton(ns("downloadPlot"), "Download Plots"),
               uiOutput(ns("gene_name")),
               plotOutput(ns("gene_plot")))
    )
  )
}
