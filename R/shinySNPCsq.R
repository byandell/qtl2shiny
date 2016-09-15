#' Shiny SNP Consequence
#'
#' @param input,output,session standard shiny arguments
#' @param chr_pos,top_snps_tbl,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#' 
#' @return tbl with top SNPs
#'
#' @export
shinySNPCsq <- function(input, output, session,
                        chr_pos, top_snps_tbl,
                        snp_action = reactive({"basic"})) {
  ns <- session$ns

  ## Shiny Modules
  
  ## Gene Region
  feature_file <- reactive({file.path(datapath, "mgi_db.sqlite")})
  callModule(shinyGeneRegion, "gene_region",
             top_snps_tbl, feature_file)
  ## Gene Exon
  gene_exon_tbl <- reactive({
    withProgress(message = 'Gene Exon Calc ...', value = 0, {
      setProgress(1)
      sql_filename <- req(feature_file())
      top_snps_set <- req(top_snps_tbl())
      get_gene_exon_snp(top_snps_set, sql_filename)
    })
  })
  callModule(shinyGeneExon, "gene_exon",
             chr_pos, top_snps_tbl, gene_exon_tbl)
  
  output$snp_csq <- renderUI({
    switch(req(input$snp_csq),
           Scan = {
             shinyScan1SNPOutput(ns("snp_scan"))
           },
           Genes = {
             shinyGeneRegionOutput(ns("gene_region"))
             },
           Exons = {
             shinyGeneExonOutput(ns("gene_exon"))
             },
           Summary = {
             shinySNPOutput(ns("best_snp"))
             })
  })
  output$snp_choice <- renderUI({
    switch(req(input$snp_csq),
           Scan = {
             shinyScan1SNPUI(ns("snp_scan"))
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
  output$title <- renderUI({
    if(snp_action() == "basic")
      h4(strong("SNP Association"))
  })
  
  gene_exon_tbl
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinySNPCsq
#' @export
shinySNPCsqUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("title")),
    radioButtons(ns("snp_csq"), "",
                 c("Scan", "Genes", "Exons", "Summary")),
    uiOutput(ns("snp_choice")))
}
#' @rdname shinySNPCsq
#' @export
shinySNPCsqOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_csq"))
}
