#' Shiny scan1 SNP analysis and plot module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param snp_par,chr_pos,pheno_names,snp_scan_obj,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom shiny NS reactive req 
#'   plotOutput
#'   renderPlot
#'   withProgress setProgress
#'   downloadButton downloadHandler
shinyScan1SNP <- function(input, output, session,
                          snp_par, chr_pos, pheno_names,
                          snp_scan_obj,
                          snp_action = shiny::reactive({"basic"})) {
  ns <- session$ns

  output$snpPlot <- shiny::renderPlot({
    if(is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    shiny::withProgress(message = 'SNP plots ...', value = 0, {
      shiny::setProgress(1)
    top_snp_asso(snp_scan_obj(), 
                 snp_par$scan_window, snp_action())
    })
  })
  
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("snp_scan_", chr_pos(), "_", snp_action(), ".pdf")) },
    content = function(file) {
      scans <- shiny::req(snp_scan_obj())
      snp_w <- shiny::req(snp_par$scan_window)
      phenos <- shiny::req(pheno_names())
      pdf(file, width = 9)
      print(top_snp_asso(scans, snp_w))
      dev.off()
    }
  )
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyScan1SNP
#' @export
shinyScan1SNPUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::downloadButton(ns("downloadPlot"), "Plots")
}
#' @rdname shinyScan1SNP
#' @export
shinyScan1SNPOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::plotOutput(ns("snpPlot"))
}
