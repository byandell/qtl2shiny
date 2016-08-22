#' Shiny SNP Consequence
#'
#' @param input,output,session standard shiny arguments
#' @param snp_scan_obj,chr_pos reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#' 
#' @return tbl with top SNPs
#'
#' @export
shinySNPCsq <- function(input, output, session,
                        snp_scan_obj, chr_pos) {
  ns <- session$ns
  
  top_snps_tbl <- reactive({
    withProgress(message = 'Get Top SNPs ...', value = 0, {
      setProgress(1)
      get_top_snps_tbl(snp_scan_obj())
    })
  })
  callModule(shinySNP, "best_snp", chr_pos, top_snps_tbl)
  feature_file <- reactive({file.path(datapath, "mgi_db.sqlite")})
  callModule(shinyGeneRegion, "gene_region",
             top_snps_tbl, feature_file)
  gene_exon_tbl <- reactive({
    withProgress(message = 'Gene Exon Calc ...', value = 0, {
      setProgress(1)
      sql_filename <- req(feature_file())
      top_snps_set <- req(top_snps_tbl())
      get_gene_exon_snp(top_snps_set, sql_filename)
    })
  })
  callModule(shinyTopFeature, "top_feature",
             chr_pos, snp_scan_obj, top_snps_tbl, gene_exon_tbl)
  callModule(shinyGeneExon, "gene_exon",
             chr_pos, top_snps_tbl, gene_exon_tbl)
  
  output$snp_csq <- renderUI({
    switch(req(input$snp_csq),
           "Top Features" = {
             shinyTopFeatureUI(ns("top_feature"))
             },
           Genes = {
             shinyGeneRegionUI(ns("gene_region"))
             },
           Exons = {
             shinyGeneExonUI(ns("gene_exon"))
             },
          Summary = {
            shinySNPUI(ns("best_snp"))
            })
  })
  top_snps_tbl
}
#' @rdname shinySNPCsq
#' @export
shinySNPCsqUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(3, 
             h4(strong("SNP Consequence")),
             radioButtons(ns("snp_csq"), "",
                          c("Top Features","Genes",
                            "Exons", "Summary"))),
      column(9,
             uiOutput(ns("snp_csq"))))
  )
}
