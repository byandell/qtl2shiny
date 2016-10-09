#' Shiny SNP Consequence
#'
#' @param input,output,session standard shiny arguments
#' @param snp_par,chr_pos,pheno_names,snp_scan_obj,top_snps_tbl,gene_exon_tbl,data_path,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#' 
#' @return tbl with top SNPs
#'
#' @export
shinySNPAssoc <- function(input, output, session,
                        snp_par, chr_pos, pheno_names,
                        snp_scan_obj, top_snps_tbl, 
                        gene_exon_tbl, data_path,
                        snp_action = reactive({"basic"})) {
  ns <- session$ns
  feature_file <- reactive({file.path(data_path(), 
                                      "mgi_db.sqlite")})
  
  ## Shiny Modules
  ## SNP Association Scan
  callModule(shinyScan1SNP, "snp_scan",
             snp_par, chr_pos, pheno_names,
             snp_scan_obj, snp_action)
  ## SNP Summary
  callModule(shinySNP, "best_snp", 
             chr_pos, top_snps_tbl, snp_action)
  ## Gene Region
  callModule(shinyGeneRegion, "gene_region",
             snp_par, 
             top_snps_tbl, feature_file,
             snp_action)
  ## Genes and Exons
  callModule(shinyGeneExon, "gene_exon",
             snp_par, chr_pos, 
             top_snps_tbl, gene_exon_tbl,
             snp_action)
  
  output$snp_input <- renderUI({
    switch(req(input$button),
           Exons   = shinyGeneExonInput(ns("gene_exon")),
           Summary = shinySNPInput(ns("best_snp")))
  })
  output$snp_output <- renderUI({
    switch(req(input$button),
           Scan    = shinyScan1SNPOutput(ns("snp_scan")),
           Genes   = shinyGeneRegionOutput(ns("gene_region")),
           Exons   = shinyGeneExonOutput(ns("gene_exon")),
           Summary = shinySNPOutput(ns("best_snp")))
  })
  
  ## Downloads
  output$download_csv_plot <- renderUI({
    switch(req(input$button),
           Scan    =,
           Summary = tagList(fluidRow(
             column(6, shinySNPUI(ns("best_snp"))),
             column(6, shinyScan1SNPUI(ns("snp_scan"))))),
           Genes   = shinyGeneRegionUI(ns("gene_region")),
           Exons   = shinyGeneExonUI(ns("gene_exon")))
  })
  output$radio <- renderUI({
    radioButtons(ns("button"), "",
                 c("Scan", "Genes", "Exons", "Summary"),
                 input$button)
  })
  input
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinySNPAssoc
#' @export
shinySNPAssocInput <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("radio")),
    uiOutput(ns("snp_input"))
  )
}
#' @rdname shinySNPAssoc
#' @export
shinySNPAssocUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("download_csv_plot"))
}
#' @rdname shinySNPAssoc
#' @export
shinySNPAssocOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_output"))
}
