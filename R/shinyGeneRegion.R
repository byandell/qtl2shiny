#' Shiny Genes in SNP Region module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param top_snps_tbl reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom ggplot2 autoplot ggtitle
#' @importFrom shiny NS reactive req isTruthy
#'   checkboxInput
#'   tableOutput plotOutput uiOutput
#'   renderPlot renderTable renderUI
#'   fluidRow column tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
shinyGeneRegion <- function(input, output, session,
                     snp_par, top_snps_tbl,
                     snp_action = shiny::reactive({"basic"})) {
  ns <- session$ns

  rng <- shiny::reactive({
    range(shiny::req(top_snps_tbl())$pos)
  })
  chr_id <- shiny::reactive({
    unique(top_snps_tbl()$chr)[1]
  })
  gene_region_tbl <- shiny::reactive({
    wrng <- rng()
    shiny::withProgress(message = 'Extract gene features ...',
                 value = 0, {
      shiny::setProgress(1)
      CCSanger::get_gene(chr_id(), wrng[1], wrng[2])
    })
  })
  output$gene_sum <- shiny::renderTable({
    summary(gene_region_tbl())
  })
  chr_pos_all <- shiny::reactive({
    chr <- shiny::req(chr_id())
    scan_w <- round(rng(), 2)
    paste(chr, scan_w[1], scan_w[2], sep = "_")
  })
  chr_pos <- shiny::reactive({
    chr <- shiny::req(chr_id())
    scan_w <- shiny::req(snp_par$scan_window)
    scan_w <- round(scan_w, 2)
    paste(chr, scan_w[1], scan_w[2], sep = "_")
  })

  output$SNP <- shiny::renderUI({
    shiny::checkboxInput(ns("SNP"), "Add SNPs?", input$SNP)
  })
  output$gene_plot <- shiny::renderPlot({
    shiny::req(gene_region_tbl(), top_snps_tbl())
    wrng <- shiny::req(snp_par$scan_window)
    phename <- shiny::req(snp_par$pheno_name)
    use_snp <- shiny::isTruthy(input$SNP)
    plot_gene_region(phename, gene_region_tbl(), top_snps_tbl(), 
                     wrng, use_snp, snp_action())
  })
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("gene_region_", chr_pos_all(), "_", snp_action(), ".csv")) },
    content = function(file) {
      write.csv(shiny::req(gene_region_tbl()), file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("gene_region_", chr_pos(), "_", snp_action(), ".pdf")) },
    content = function(file) {
      shiny::req(gene_region_tbl(), top_snps_tbl())
      wrng <- shiny::req(snp_par$scan_window)
      phename <- shiny::req(snp_par$pheno_name)
      use_snp <- shiny::isTruthy(input$SNP)
      pheno_names <- unique(shiny::req(top_snps_tbl())$pheno)
      pdf(file, width = 9)
      for(pheno in pheno_names)
        print(plot_gene_region(pheno, gene_region_tbl(), top_snps_tbl(), 
                       wrng, use_snp, snp_action()))
      dev.off()
    }
  )
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyGeneRegion
#' @export
shinyGeneRegionInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("SNP"))
}
#' @rdname shinyGeneRegion
#' @export
shinyGeneRegionUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::fluidRow(
      shiny::column(6, shiny::downloadButton(ns("downloadData"), "CSV")),
      shiny::column(6, shiny::downloadButton(ns("downloadPlot"), "Plot"))))
}
#' @rdname shinyGeneRegion
#' @export
shinyGeneRegionOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::plotOutput(ns("gene_plot")),
    shiny::tableOutput(ns("gene_sum"))
  )
}
