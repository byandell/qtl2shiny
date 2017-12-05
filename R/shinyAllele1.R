#' Shiny coefficient analysis and plot module
#'
#' Shiny module for scan1 coefficient plots.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,probs_obj,K_chr,analyses_df,patterns,scan_pat,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom qtl2pattern allele1
#' @importFrom shiny NS reactive req 
#'   radioButtons selectInput sliderInput updateSliderInput
#'   dataTableOutput plotOutput uiOutput
#'   renderDataTable renderPlot renderUI
#'   fluidRow column strong tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
#' @importFrom ggplot2 autoplot ggtitle
#' 
shinyAllele1 <- function(input, output, session,
                  win_par, 
                  phe_df, cov_mx, probs_obj, K_chr, analyses_df, 
                  patterns, scan_pat, snp_action) {
  ns <- session$ns
  
  # Scan Window slider
  output$pos_Mbp <- shiny::renderUI({
    chr_id <- shiny::req(win_par$chr_id)
    map <- shiny::req(probs_obj())$map[[chr_id]]
    rng <- round(2 * range(map)) / 2
    if(is.null(selected <- input$pos_Mbp))
      selected <- req(win_par$peak_Mbp)
    shiny::sliderInput(ns("pos_Mbp"), NULL, rng[1], rng[2],
                selected, step=.1)
  })
  ## Reset pos_Mbp if chromosome changes.
  observeEvent(win_par$chr_id, {
    map <- shiny::req(probs_obj()$map)
    chr <- shiny::req(win_par$chr_id)
    rng <- round(2 * range(map[[chr]])) / 2
    shiny::updateSliderInput(session, "pos_Mbp", NULL, 
                             req(win_par$peak_Mbp), 
                             rng[1], rng[2], step=.1)
  })
  
  ## Coefficient Effects.
  allele_obj <- shiny::reactive({
    shiny::req(snp_action())
    shiny::req(phe_df(), probs_obj(), K_chr(), cov_mx())
    blups <- attr(scan_pat(), "blups")
    shiny::withProgress(message = 'Effect scans ...', value = 0, {
      shiny::setProgress(1)
      qtl2pattern::allele1(phe_df(), cov_mx(), probs_obj()$probs, 
                           probs_obj()$map, K_chr(),
                           patterns = patterns(),
                           scan_pat = scan_pat(),
                           blups = blups)
      })
  })
  output$allele1Plot <- shiny::renderPlot({
    shiny::req(allele_obj(), input$pos_Mbp)
    shiny::withProgress(message = 'Allele plots ...', value = 0, {
      shiny::setProgress(1)
      p <- ggplot2::autoplot(allele_obj(), pos = input$pos_Mbp)
      if(is.null(p)) {
        plot_null()
      } else {
        p + ggplot2::ggtitle(names(phe_df()))
      }
    })
  })
  output$allele1Sum <- shiny::renderTable({
    shiny::req(allele_obj(), input$pos_Mbp)
    shiny::withProgress(message = 'Effect summary ...', value = 0, {
      shiny::setProgress(1)
      summary(allele_obj(), pos = input$pos_Mbp)
    })
  })
  
  ## Downloads.
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("allele1_", win_par$chr_id, "_", win_par$peak_Mbp, ".csv")) },
    content = function(file) {
      shiny::req(allele_obj())
      out <- tidyr::spread(
        dplyr::select(
          dplyr::mutate(allele_obj(),
                        allele = paste(source, allele, sep = ".")),
          -probe, -source),
        allele, effect)
      write.csv(out, file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("allele1_", win_par$chr_id, "_", win_par$peak_Mbp, ".pdf")) },
    content = function(file) {
      shiny::req(allele_obj(), input$pos_Mbp)
      pdf(file, width=9,height=9)
      print(ggplot2::autoplot(
        allele_obj(), pos = input$pos_Mbp) +
        ggplot2::ggtitle(names(phe_df())))
      dev.off()
    }
  )
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyAllele1
#' @export
shinyAllele1UI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("pos_Mbp")),
    shiny::fluidRow(
      shiny::column(6, shiny::downloadButton(ns("downloadData"), "CSV")),
      shiny::column(6, shiny::downloadButton(ns("downloadPlot"), "Plots"))))
}
#' @rdname shinyAllele1
#' @export
shinyAllele1Output <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::plotOutput(ns("allele1Plot")),
    shiny::tableOutput(ns("allele1Sum")))
}
