#' Shiny Top Features in SNP Region module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param scan_obj,top_snps_tbl,gene_exon_tbl reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyTopFeature <- function(input, output, session,
                     chr_pos, scan_obj, top_snps_tbl, gene_exon_tbl) {
  ns <- session$ns
  ## This is just repeat of GeneRegion for now
  ## See child/tests/top_feature.R for filler.
  top_feature <- reactive({
    req(top_snps_tbl(),scan_obj(),gene_exon_tbl())
    withProgress(message = 'Merging gene info ...', value = 0,
    {
      setProgress(1)
      merge_feature(top_snps_tbl(), scan_obj(), 1.5, 0, gene_exon_tbl())
    })
  })
  output$top_snp_type <- renderDataTable({
    tops <- req(top_feature(), "SNP type")
    summary(tops)
  }, options = list(scrollX = TRUE, paging = FALSE, searching=FALSE))
  output$top_pattern <- renderDataTable({
    summary(top_feature(), "pattern")
  }, options = list(scrollX = TRUE, paging = FALSE, searching=FALSE))
  phename <- reactive({dimnames(scan_obj()$lod)[[2]]})
  output$top_names <- renderUI({
    selectInput(ns("top_name"), NULL, phename())
  })
  output$top_gene_by_snp <- renderPlot({
    req(top_feature(), input$top_name)
    plot(top_feature(), input$top_name, "consequence")
  })
  output$top_gene_by_pattern <- renderPlot({
    req(top_feature(),input$top_name)
    plot(top_feature(), input$top_name, "pattern")
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("top_feature_", chr_pos(), ".csv")) },
    content = function(file) {
      write.csv(req(top_feature()), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("top_feature_", chr_pos(), ".pdf")) },
    content = function(file) {
      req(top_feature())
      pdf(file)
      for(phenoi in phename()) {
        print(plot_top_feat_csq(phenoi))
        print(plot_top_feat_pat(phenoi))
      }
      dev.off()
    }
  )
}

#' UI for shinyTopFeature Shiny Module
#'
#' UI for top features to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyTopFeature
#' @export
shinyTopFeatureUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      uiOutput(ns("top_names")),
      downloadButton(ns("downloadPlot"), "Plots"),
      downloadButton(ns("downloadData"), "CSV")),
    tabsetPanel(
      tabPanel("By Allele Pattern",
               plotOutput(ns("top_gene_by_snp")),
               dataTableOutput(ns("top_pattern"))),
      tabPanel("By Consequence",
               plotOutput(ns("top_gene_by_pattern")),
               dataTableOutput(ns("top_snp_type"))))
  )
}
