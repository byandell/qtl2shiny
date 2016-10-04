#' Shiny Genes in SNP Region module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param top_snps_tbl,feature_file reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyGeneRegion <- function(input, output, session,
                     snp_par, top_snps_tbl, feature_file) {
  ns <- session$ns

  rng <- reactive({range(req(top_snps_tbl())$pos_Mbp)})
  chr_id <- reactive({unique(top_snps_tbl()$chr)[1]})
  gene_region_tbl <- reactive({
    req(top_snps_tbl())
    wrng <- rng()
    withProgress(message = 'Extract gene features ...',
                 value = 0, {
      setProgress(1)
      get_mgi_features(chr_id(), wrng[1], wrng[2],
                       sql_file = feature_file())
    })
  })
  output$gene_sum <- renderTable({
    summary(gene_region_tbl())
  })
  chr_pos_all <- reactive({
    chr <- req(chr_id())
    scan_w <- round(rng(), 2)
    paste(chr, scan_w[1], scan_w[2], sep = "_")
  })
  chr_pos <- reactive({
    chr <- req(chr_id())
    scan_w <- req(snp_par$scan_window)
    scan_w <- round(scan_w, 2)
    paste(chr, scan_w[1], scan_w[2], sep = "_")
  })

  plot_gene_region <- function() {
    req(gene_region_tbl(),snp_par$scan_window)
    ## Plot pseudogene and gene locations along with SNPs.
    ## Ordered by pseudogenes first, then genes
    ## negative (blue) strand, then unknown (grey), then positive (red) strand.
    wrng <- snp_par$scan_window
    pheno <- top_snps_tbl()$pheno[1]
    ## Filtering removes feature_tbl class, so need to be explicit.
    plot(subset(gene_region_tbl(), wrng[1], wrng[2]),
         top_snps_tbl=subset(top_snps_tbl(), wrng[1], wrng[2])) +
      ggtitle(paste("SNPs for", pheno))
  }
  output$gene_plot <- renderPlot({
    plot_gene_region()
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("gene_region_", chr_pos_all(), ".csv")) },
    content = function(file) {
      write.csv(req(gene_region_tbl()), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    ## This only gets plot with selected pheno_name.
    ## Would have to move pheno_name subsetting in here to change.
    filename = function() {
      file.path(paste0("gene_region_", chr_pos(), ".pdf")) },
    content = function(file) {
      ggsave(file, plot = plot_gene_region(), device = "pdf")
    }
  )
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyGeneRegion
#' @export
shinyGeneRegionUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    fluidRow(
      column(6, downloadButton(ns("downloadData"), "CSV")),
      column(6, downloadButton(ns("downloadPlot"), "Plot"))))
}
#' @rdname shinyGeneRegion
#' @export
shinyGeneRegionOutput <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("gene_plot")),
    tableOutput(ns("gene_sum"))
  )
}
