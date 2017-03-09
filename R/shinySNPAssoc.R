#' Shiny SNP Consequence
#'
#' @param input,output,session standard shiny arguments
#' @param snp_par,chr_pos,pheno_names,snp_scan_obj,snpinfo,top_snps_tbl,gene_exon_tbl,data_path,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#' 
#' @return tbl with top SNPs
#'
#' @export
#' @importFrom shiny callModule NS reactive req 
#'   radioButtons
#'   uiOutput
#'   renderUI
#'   fluidRow column tagList
shinySNPAssoc <- function(input, output, session,
                        snp_par, chr_pos, pheno_names,
                        snp_scan_obj, snpinfo, top_snps_tbl, 
                        gene_exon_tbl, data_path,
                        snp_action = shiny::reactive({"basic"})) {
  ns <- session$ns
  feature_file <- shiny::reactive({file.path(data_path(), 
                                      "mgi_db.sqlite")})
  
  ## Shiny Modules
  ## SNP Association Scan
  shiny::callModule(shinyScan1SNP, "snp_scan",
             snp_par, chr_pos, pheno_names,
             snp_scan_obj, snpinfo, snp_action)
  ## SNP Summary
  shiny::callModule(shinySNP, "best_snp", 
             chr_pos, top_snps_tbl, snp_action)
  ## Gene Region
  shiny::callModule(shinyGeneRegion, "gene_region",
             snp_par, 
             top_snps_tbl, feature_file,
             snp_action)
  ## Genes and Exons
  shiny::callModule(shinyGeneExon, "gene_exon",
             snp_par, chr_pos, 
             top_snps_tbl, gene_exon_tbl,
             snp_action)
  
  output$snp_check <- shiny::renderUI({
    switch(shiny::req(input$button),
           Genes   = shinyGeneRegionInput(ns("gene_region")))
  })
  output$snp_input <- shiny::renderUI({
    switch(shiny::req(input$button),
           Exons   = shinyGeneExonInput(ns("gene_exon")),
           Summary = shinySNPInput(ns("best_snp")))
  })
  output$snp_output <- shiny::renderUI({
    switch(shiny::req(input$button),
           Scan    = shinyScan1SNPOutput(ns("snp_scan")),
           Genes   = shinyGeneRegionOutput(ns("gene_region")),
           Exons   = shinyGeneExonOutput(ns("gene_exon")),
           Summary = shinySNPOutput(ns("best_snp")))
  })
  
  ## Downloads
  output$download_csv_plot <- shiny::renderUI({
    switch(shiny::req(input$button),
           Scan    =,
           Summary = shiny::tagList(shiny::fluidRow(
             shiny::column(6, shinySNPUI(ns("best_snp"))),
             shiny::column(6, shinyScan1SNPUI(ns("snp_scan"))))),
           Genes   = shinyGeneRegionUI(ns("gene_region")),
           Exons   = shinyGeneExonUI(ns("gene_exon")))
  })
  output$radio <- shiny::renderUI({
    shiny::radioButtons(ns("button"), "",
                 c("Scan", "Genes", "Exons", "Summary"),
                 input$button)
  })
  input
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinySNPAssoc
#' @export
shinySNPAssocInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(6, shiny::uiOutput(ns("radio"))),
      shiny::column(6, shiny::uiOutput(ns("snp_check")))),
    shiny::uiOutput(ns("snp_input"))
  )
}
#' @rdname shinySNPAssoc
#' @export
shinySNPAssocUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("download_csv_plot"))
}
#' @rdname shinySNPAssoc
#' @export
shinySNPAssocOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("snp_output"))
}
