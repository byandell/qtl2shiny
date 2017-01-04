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
#' @importFrom doqtl2 get_mgi_features
#' @importFrom ggplot2 ggtitle
#' @importFrom shiny NS reactive req isTruthy
#'   checkboxInput
#'   tableOutput plotOutput uiOutput
#'   renderPlot renderTable renderUI
#'   fluidRow column tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
shinyGeneRegion <- function(input, output, session,
                     snp_par, top_snps_tbl, feature_file,
                     snp_action = shiny::reactive({"basic"})) {
  ns <- session$ns

  rng <- shiny::reactive({
    range(shiny::req(top_snps_tbl())$pos_Mbp)
  })
  chr_id <- shiny::reactive({
    unique(top_snps_tbl()$chr)[1]
  })
  gene_region_tbl <- shiny::reactive({
    wrng <- rng()
    shiny::withProgress(message = 'Extract gene features ...',
                 value = 0, {
      shiny::setProgress(1)
      doqtl2::get_mgi_features(chr_id(), wrng[1], wrng[2],
                       sql_file = feature_file())
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

  plot_gene_region <- function(pheno) {
    shiny::req(gene_region_tbl(),snp_par$scan_window)
    ## Plot pseudogene and gene locations along with SNPs.
    ## Ordered by pseudogenes first, then genes
    ## negative (blue) strand, then unknown (grey), then positive (red) strand.
    wrng <- snp_par$scan_window
    ## Filtering removes feature_tbl class, so need to be explicit.
    use_snp <- shiny::isTruthy(input$SNP)
    if(use_snp) {
      top_snps_rng <- subset(top_snps_tbl(), 
                             wrng[1], wrng[2],
                             pheno)
      if(!nrow(top_snps_rng))
        top_snps_rng <- NULL
    } else {
      top_snps_rng <- NULL
    }
    p <- plot(subset(gene_region_tbl(), wrng[1], wrng[2]),
              top_snps_tbl = top_snps_rng)
    if(use_snp)
      p <- p + 
        ggplot2::ggtitle(paste("SNPs for", pheno, snp_action()))
    p
  }
  output$SNP <- shiny::renderUI({
    shiny::checkboxInput(ns("SNP"), "Add SNPs?", input$SNP)
  })
  output$gene_plot <- shiny::renderPlot({
    phename <- shiny::req(snp_par$pheno_name)
    plot_gene_region(phename)
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
      shiny::req(gene_region_tbl(),snp_par$scan_window)
      pheno_names <- unique(shiny::req(top_snps_tbl())$pheno)
      pdf(file, width = 9)
      for(pheno in pheno_names)
        print(plot_gene_region(pheno))
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
