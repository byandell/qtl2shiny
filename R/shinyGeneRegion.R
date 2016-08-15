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
                     top_snps_tbl, feature_file) {
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

  # Scan Window slider
  output$choose_scan_window <- renderUI({
    wrng <- round(rng(), 2)
    sliderInput(ns("scan_window"), NULL, wrng[1], wrng[2],
                wrng, step=.5)
  })
  chr_pos_all <- reactive({
    chr <- req(chr_id())
    scan_w <- round(rng(), 2)
    paste(chr, scan_w[1], scan_w[2], sep = "_")
  })
  chr_pos <- reactive({
    chr <- req(chr_id())
    scan_w <- req(input$scan_window)
    scan_w <- round(scan_w, 2)
    paste(chr, scan_w[1], scan_w[2], sep = "_")
  })

  plot_gene_region <- function() {
    req(gene_region_tbl(),input$scan_window)
    ## Plot pseudogene and gene locations along with SNPs.
    ## Ordered by pseudogenes first, then genes
    ## negative (blue) strand, then unknown (grey), then positive (red) strand.
    wrng <- input$scan_window
    ## Filtering removes feature_tbl class, so need to be explicit.
    plot(subset(gene_region_tbl(), wrng[1], wrng[2]),
         top_snps_tbl=subset(top_snps_tbl(), wrng[1], wrng[2]))
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
    filename = function() {
      file.path(paste0("gene_region_", chr_pos(), ".pdf")) },
    content = function(file) {
      ggsave(file, plot = plot_gene_region(), device = "pdf")
    }
  )
  reactive({input$scan_window})
}

#' UI for shinyGeneRegion Shiny Module
#'
#' UI for scan1 analyses and summaries to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyGeneRegion
#' @export
shinyGeneRegionUI <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      tabPanel("summary",
               downloadButton(ns("downloadData"), "Download CSV"),
               tableOutput(ns("gene_sum"))),
      tabPanel("plot",
               uiOutput(ns("choose_scan_window")),
               downloadButton(ns("downloadPlot"), "Download Plot"),
               plotOutput(ns("gene_plot")))
    )
  )
}
