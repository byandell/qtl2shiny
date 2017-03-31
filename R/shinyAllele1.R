#' Shiny coefficient analysis and plot module
#'
#' Shiny module for scan1 coefficient plots.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,probs_obj,K_chr,analyses_df,patterns,snp_action reactive arguments
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
shinyAllele1 <- function(input, output, session,
                  win_par, 
                  phe_df, cov_mx, probs_obj, K_chr, analyses_df, patterns, 
                  snp_action) {
  ns <- session$ns
  
  # Scan Window slider
  output$scan_window <- shiny::renderUI({
    chr_id <- shiny::req(win_par$chr_id)
    map <- shiny::req(probs_obj())$map[[chr_id]]
    rng <- round(2 * range(map)) / 2
    if(is.null(selected <- input$scan_window))
      selected <- mean(rng)
    shiny::sliderInput(ns("scan_window"), NULL, rng[1], rng[2],
                selected, step=.1)
  })
  ## Reset scan_window if chromosome changes.
  observeEvent(win_par$chr_id, {
    map <- shiny::req(probs_obj()$map)
    chr <- shiny::req(win_par$chr_id)
    rng <- round(2 * range(map[[chr]])) / 2
    shiny::updateSliderInput(session, "scan_window", NULL, mean(rng), 
                      rng[1], rng[2], step=.1)
  })
  
  ## Select phenotype for plots.
  output$pheno_name <- shiny::renderUI({
    shiny::req(phe_df())
    shiny::selectInput(ns("pheno_name"), NULL,
                choices = names(phe_df()))
  })
  
  ## Coefficient Effects.
  allele_obj <- shiny::reactive({
    shiny::req(phe_df(), probs_obj(), K_chr(), cov_mx(), patterns())
    pheno <- shiny::req(input$pheno_name)
    shiny::withProgress(message = 'Effect scans ...', value = 0, {
      shiny::setProgress(1)
      qtl2pattern::allele1(phe_df()[, pheno, drop = FALSE], 
                           cov_mx(), probs_obj()$probs, 
                           probs_obj()$map, K_chr(), patterns())
      })
  })
  output$allele1Plot <- shiny::renderPlot({
    shiny::req(input$pheno_name, allele_obj(), input$scan_window)
    shiny::withProgress(message = 'Allele plots ...', value = 0, {
      shiny::setProgress(1)
      plot(allele_obj(), pos = input$scan_window)
    })
  })
  output$allele1Sum <- shiny::renderTable({
    shiny::req(input$pheno_name, allele_obj(), input$scan_window)
    shiny::withProgress(message = 'Effect summary ...', value = 0, {
      shiny::setProgress(1)
      summary(allele_obj(), pos = input$scan_window)
    })
  })
  
  ## Downloads.
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("sum_effects_", win_par$chr_id, ".csv")) },
    content = function(file) {
      shiny::req(eff_obj(), scan_obj(), probs_obj())
      write.csv(summary(eff_obj(), scan_obj(), probs_obj()$map), file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("scan_", win_par$chr_id, ".pdf")) },
    content = function(file) {
      shiny::req(win_par$chr_id)
      effs <- shiny::req(eff_obj())
      scans <- shiny::req(scan_obj())
      win <- shiny::req(input$scan_window)
      map <- shiny::req(probs_obj())$map
      pdf(file, width=9,height=9)
      print(plot(scans, map,
                 lodcolumn = seq_along(names(effs)),
                 chr = win_par$chr_id,
                 xlim = win))
      for(pheno in names(effs)) {
        plot_eff(pheno, effs, map, scans, win,
                 addlod = TRUE)
      }
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
    shiny::uiOutput(ns("pheno_name")),
    shiny::uiOutput(ns("scan_window")),
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
