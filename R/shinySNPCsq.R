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
                        snp_par, chr_pos, top_snps_tbl,
                        gene_exon_tbl,
                        snp_action = reactive({"basic"})) {
  ns <- session$ns

  ## Shiny Modules
  feature_file <- reactive({file.path(datapath, "mgi_db.sqlite")})
  ## Gene Region
  callModule(shinyGeneRegion, "gene_region",
             snp_par, top_snps_tbl, feature_file)
  ## Genes and Exons
  callModule(shinyGeneExon, "gene_exon",
             chr_pos, top_snps_tbl, gene_exon_tbl)
  
  output$snp_csq <- renderUI({
    switch(req(snp_par$snp_assoc),
           Genes = {
             shinyGeneRegionOutput(ns("gene_region"))
             },
           Exons = {
             shinyGeneExonOutput(ns("gene_exon"))
             })
  })
  output$snp_choice <- renderUI({
    switch(req(snp_par$snp_assoc),
           Genes = {
             shinyGeneRegionUI(ns("gene_region"))
           },
           Exons = {
             shinyGeneExonUI(ns("gene_exon"))
           })
  })
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_scan_", chr_pos(), ".pdf")) },
    content = function(file) {
      scans <- req(scan_obj())
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
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinySNPCsq
#' @export
shinySNPCsqInput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_choice"))
}
#' @rdname shinySNPCsq
#' @export
shinySNPCsqUI <- function(id) {
  ns <- NS(id)
  downloadButton(ns("downloadPlot"), "Plots")
}
#' @rdname shinySNPCsq
#' @export
shinySNPCsqOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_csq"))
}
