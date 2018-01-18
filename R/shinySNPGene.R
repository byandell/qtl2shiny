#' Shiny SNP Association
#' 
#' Shiny module for SNP association mapping, with interfaces \code{shinySNPGeneInput}, \code{shinySNPGeneUI} and  \code{shinySNPGeneOutput}.
#'
#' @param input,output,session standard shiny arguments
#' @param snp_par,chr_pos,pheno_names,snp_scan_obj,snpinfo,top_snps_tbl,gene_exon_tbl,project_info,snp_action reactive arguments
#' @param id shiny identifier
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
shinySNPGene <- function(input, output, session,
                        snp_par, chr_pos, pheno_names,
                        snp_scan_obj, snpinfo, top_snps_tbl, 
                        gene_exon_tbl,
                        project_info,
                        snp_action = shiny::reactive({"basic"})) {
  ns <- session$ns

  ## Shiny Modules
  ## SNP Association Scan
  shiny::callModule(shinySNPPlot, "snp_scan",
             snp_par, chr_pos, pheno_names,
             snp_scan_obj, snpinfo, snp_action)
  ## SNP Summary
  shiny::callModule(shinySNPSum, "best_snp", 
             chr_pos, top_snps_tbl, snp_action)
  ## Gene Region
  shiny::callModule(shinyGeneRegion, "gene_region",
             snp_par, 
             top_snps_tbl,
             project_info,
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
           Summary = shinySNPSumInput(ns("best_snp")))
  })
  output$snp_output <- shiny::renderUI({
    switch(shiny::req(input$button),
           Scan    = shinySNPPlotOutput(ns("snp_scan")),
           Genes   = shinyGeneRegionOutput(ns("gene_region")),
           Exons   = shinyGeneExonOutput(ns("gene_exon")),
           Summary = shinySNPSumOutput(ns("best_snp")))
  })
  
  ## Downloads
  output$download_csv_plot <- shiny::renderUI({
    switch(shiny::req(input$button),
           Scan    =,
           Summary = shiny::tagList(shiny::fluidRow(
             shiny::column(6, shinySNPSumUI(ns("best_snp"))),
             shiny::column(6, shinySNPPlotUI(ns("snp_scan"))))),
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

shinySNPGeneInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(6, shiny::uiOutput(ns("radio"))),
      shiny::column(6, shiny::uiOutput(ns("snp_check")))),
    shiny::uiOutput(ns("snp_input"))
  )
}
shinySNPGeneUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("download_csv_plot"))
}
shinySNPGeneOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("snp_output"))
}
