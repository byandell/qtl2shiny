#' Shiny coefficient analysis and plot module
#'
#' Shiny module for scan1 coefficient plots.
#'
#' @param input,output,session standard shiny arguments
#' @param chr_id,phe_df,cov_mx,scan_window,chr_pos,pheno_id,scan_obj,probs_chr,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return 2-element vector of scan window
#'
#' @export
shinyScan1Plot <- function(input, output, session,
                       chr_id, phe_df, cov_mx,
                       scan_window, chr_pos, pheno_id,
                       scan_obj, probs_chr, K_chr) {
  ## Scan1 plot
  output$scanPlot <- renderPlot({
    withProgress(message = 'Scan plots ...', value = 0,
    {
      setProgress(1)
      show_peaks(chr_id(), scan_obj(), mytitle="", xlim=scan_window())
    })
  })

  ## Coefficient Effects Reactives.
  eff_obj <- reactive({
    withProgress(message = 'Effect scans ...', value = 0,
    {
      setProgress(1)
      listof_scan1coefCC(phe_df(), cov_mx(), probs_chr(), K_chr())
    })
  })
  plot_eff <- function(pheno) {
    lodcol <- match(pheno, names(eff_obj()))
    plot_coefCC(eff_obj()[[lodcol]], xlim=scan_window())
    max_pos <- max(scan_obj(), lodcolumn=lodcol)$pos[1]
    abline(v=max_pos, lty=2, lwd=2, col=lodcol)
    at <- par("usr")[1:2]
    at <- seq(at[1],at[2],length.out=length(CCcolors))
    mtext(names(CCcolors),col=CCcolors,at=at)
    title(pheno)
  }
  output$effPlot <- renderPlot({
    req(pheno_id(),eff_obj())
    withProgress(message = 'Effect plots ...', value = 0,
    {
      setProgress(1)
      plot_eff(pheno_id())
    })
  })
  output$effSummary <- renderDataTable({
    req(eff_obj(),scan_obj())
    withProgress(message = 'Effect summary ...', value = 0,
                 {
                   setProgress(1)
                   summary(eff_obj(),scan_obj())
                 })
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("sum_effects_", chr_pos(), ".csv")) },
    content = function(file) {
      req(eff_obj(),scan_obj())
      write.csv(summary(eff_obj(),scan_obj()), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("genome_scan_", chr_pos(), ".pdf")) },
    content = function(file) {
      req(eff_obj())
      pdf(file)
      show_peaks(chr_id(), scan_obj(), mytitle="", xlim=scan_window())
      for(pheno in names(eff_obj()))
        plot_eff(pheno)
      dev.off()
    }
  )
}

#' UI for shinyScan1 Shiny Module
#'
#' UI for scan1 analyses and summaries to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyScan1Plot
#' @export
shinyScan1PlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    downloadButton(ns("downloadPlot"), "Download Plots"),
    plotOutput(ns("effPlot")),
    plotOutput(ns("scanPlot")))
}
#' Output for shinyScan1 Shiny Module
#'
#' Output for scan1 analyses and summaries to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyScan1Plot
#' @export
shinyScan1PlotOutput <- function(id) {
  ns <- NS(id)
  tagList(
    downloadButton(ns("downloadData"), "Download CSV"),
    dataTableOutput(ns("effSummary")))
}
