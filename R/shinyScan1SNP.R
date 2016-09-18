#' Shiny scan1 SNP analysis and plot module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,pheno_anal,probs_obj,K_chr,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyScan1SNP <- function(input, output, session,
                          snp_par, chr_pos, pheno_anal, snp_scan_obj,
                          snp_action = reactive({"basic"})) {
  ns <- session$ns

  pheno_id <- reactive({
    pheno_anal <- req(snp_par$pheno_assoc)
    pheno_anal()[pheno_anal]
  })

  ## Reactives for SNP analysis.
  ## Want to have slider for snp_w
  phename <- reactive({dimnames(snp_scan_obj()$lod)[[2]]})
  output$snpPlot <- renderPlot({
    if(is.null(pheno_id()) | is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP plots ...', value = 0, {
      setProgress(1)
      top_snp_asso(pheno_id(), snp_scan_obj(), snp_par$scan_window, snp_action())
    })
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_scan_", chr_pos(), ".pdf")) },
    content = function(file) {
      scans <- req(snp_scan_obj())
      snp_w <- req(snp_par$scan_window)
      phenos <- req(phename())
      pdf(file)
      ## Plots by phenotype.
      for(pheno in phenos) {
        top_snp_asso(pheno, scans, snp_w)
      }
      dev.off()
    }
  )
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyScan1SNP
#' @export
shinyScan1SNPUI <- function(id) {
  ns <- NS(id)
  downloadButton(ns("downloadPlot"), "Plots")
}
#' @rdname shinyScan1SNP
#' @export
shinyScan1SNPOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns("snpPlot"))
}
