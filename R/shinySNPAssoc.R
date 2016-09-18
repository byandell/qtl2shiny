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
                        snp_par, chr_pos, pheno_anal, 
                        snp_scan_obj, 
                        top_snps_tbl, gene_exon_tbl,
                        snp_action = reactive({"basic"})) {
  ns <- session$ns
  phename <- reactive({dimnames(snp_scan_obj()$lod)[[2]]})
  feature_file <- reactive({file.path(datapath, 
                                      "mgi_db.sqlite")})
  
  ## Shiny Modules
  ## SNP Association Scan
  callModule(shinyScan1SNP, "snp_scan",
             snp_par, chr_pos, pheno_anal, snp_scan_obj, 
             snp_action)
  ## SNP Summary
  callModule(shinySNP, "best_snp", 
             chr_pos, top_snps_tbl)
  ## Shiny Modules
  ## Gene Region
  callModule(shinyGeneRegion, "gene_region",
             snp_par, top_snps_tbl, feature_file)
  ## Genes and Exons
  callModule(shinyGeneExon, "gene_exon",
             chr_pos, top_snps_tbl, gene_exon_tbl)
  
  output$snp_choice <- renderUI({
    switch(req(input$button),
           Genes   = shinyGeneRegionUI(ns("gene_region")),
           Exons   = shinyGeneExonUI(ns("gene_exon")),
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
    fluidRow(
      column(6, shinySNPUI(ns("best_snp"))),
      column(6, uiOutput(ns("download_snp_assoc"))))
  })
  output$download_snp_assoc <- renderUI({
    switch(req(input$button),
           Scan    =, 
           Summary = shinyScan1SNPUI(ns("snp_scan")),
           Genes   =,
           Exons   = downloadButton(ns("downloadPlot"), "Plots"))
  })
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_scan_", chr_pos(), ".pdf")) },
    content = function(file) {
      scans <- req(snp_scan_obj())
      snp_w <- req(snp_par$scan_window)
      phenos <- req(phename())
      pdf(file)
      ## Plots over all phenotypes
      print(top_pat_plot(phenos, scans, snp_w,
                         group = "pheno", snp_action = snp_action()))
      print(top_pat_plot(phenos, scans, snp_w,
                         group = "pattern", snp_action = snp_action()))
      ## Plots by phenotype.
      for(pheno in phenos) {
        print(top_pat_plot(pheno, scans, snp_w, FALSE,
                           snp_action = snp_action()))
      }
      dev.off()
    }
  )
  input
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinySNPCsq
#' @export
shinySNPCsqInput <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("button"), "",
                 c("Scan", "Genes", "Exons", "Summary")),
    uiOutput(ns("snp_choice"))
  )
}
#' @rdname shinySNPCsq
#' @export
shinySNPCsqUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("download_csv_plot"))
}
#' @rdname shinySNPCsq
#' @export
shinySNPCsqOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_output"))
}
