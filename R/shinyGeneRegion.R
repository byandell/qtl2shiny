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
                     snp_par, top_snps_tbl, feature_file,
                     snp_action = reactive({"basic"})) {
  ns <- session$ns

  rng <- reactive({
    range(req(top_snps_tbl())$pos_Mbp)
  })
  chr_id <- reactive({
    unique(top_snps_tbl()$chr)[1]
  })
  gene_region_tbl <- reactive({
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

  plot_gene_region <- function(pheno) {
    req(gene_region_tbl(),snp_par$scan_window)
    ## Plot pseudogene and gene locations along with SNPs.
    ## Ordered by pseudogenes first, then genes
    ## negative (blue) strand, then unknown (grey), then positive (red) strand.
    wrng <- snp_par$scan_window
    ## Filtering removes feature_tbl class, so need to be explicit.
    use_snp <- isTruthy(input$SNP)
    if(use_snp) {
      top_snps_rng <- subset(top_snps_tbl(), 
                             wrng[1], wrng[2],
                             pheno)
    } else {
      top_snps_rng <- NULL
    }
    p <- plot(subset(gene_region_tbl(), wrng[1], wrng[2]),
              top_snps_tbl = top_snps_rng)
    if(use_snp)
      p <- p + 
        ggtitle(paste("SNPs for", pheno, snp_action()))
    p
  }
  output$SNP <- renderUI({
    checkboxInput(ns("SNP"), "Add SNPs?", input$SNP)
  })
  output$gene_plot <- renderPlot({
    phename <- req(snp_par$pheno_name)
    plot_gene_region(phename)
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("gene_region_", chr_pos_all(), "_", snp_action(), ".csv")) },
    content = function(file) {
      write.csv(req(gene_region_tbl()), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("gene_region_", chr_pos(), "_", snp_action(), ".pdf")) },
    content = function(file) {
      req(gene_region_tbl(),snp_par$scan_window)
      pheno_names <- unique(req(top_snps_tbl())$pheno)
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
  ns <- NS(id)
  uiOutput(ns("SNP"))
}
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
