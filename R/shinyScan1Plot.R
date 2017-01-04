#' Shiny coefficient analysis and plot module
#'
#' Shiny module for scan1 coefficient plots.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,probs_obj,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom doqtl2 listof_scan1coefCC plot_eff show_peaks
#' @importFrom qtl2scan scan1
#' @importFrom shiny NS reactive req 
#'   radioButtons selectInput sliderInput updateSliderInput
#'   dataTableOutput plotOutput uiOutput
#'   renderDataTable renderPlot renderUI
#'   fluidRow column strong tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
shinyScan1Plot <- function(input, output, session,
                  win_par, pmap_obj, 
                  phe_df, cov_mx, probs_obj, K_chr) {
  ns <- session$ns
  
  ## Scan1
  scan_obj <- shiny::reactive({
    shiny::req(phe_df(), probs_obj(), K_chr(), cov_mx())
    shiny::withProgress(message = "Genome Scan ...", value = 0, {
      shiny::setProgress(1)
      qtl2scan::scan1(probs_obj(), phe_df(), K_chr(), cov_mx())
    })
  })
  
  # Scan Window slider
  output$scan_window <- shiny::renderUI({
    chr_id <- shiny::req(win_par$chr_id)
    map <- shiny::req(pmap_obj())[[chr_id]]
    rng <- round(2 * range(map)) / 2
    if(is.null(selected <- input$scan_window))
      selected <- rng
    shiny::sliderInput(ns("scan_window"), NULL, rng[1], rng[2],
                selected, step=.5)
  })
  ## Reset scan_window if chromosome changes.
  observeEvent(win_par$chr_id, {
    map <- shiny::req(pmap_obj())
    chr <- shiny::req(win_par$chr_id)
    rng <- round(2 * range(map[[chr]])) / 2
    shiny::updateSliderInput(session, "scan_window", NULL, rng, 
                      rng[1], rng[2], step=.5)
  })
  
  ## Select phenotype for plots.
  output$pheno_name <- shiny::renderUI({
    shiny::req(phe_df())
    shiny::selectInput(ns("pheno_name"), NULL,
                choices = names(phe_df()))
  })

  ## Scan1 plot
  output$scanPlot <- shiny::renderPlot({
    shiny::req(win_par$chr_id, input$scan_window, scan_obj())
    shiny::withProgress(message = 'Genome LOD Plot ...', value = 0, {
      shiny::setProgress(1)
      doqtl2::show_peaks(win_par$chr_id, scan_obj(), mytitle="", 
                 xlim=input$scan_window)
    })
  })

  ## Coefficient Effects.
  eff_obj <- shiny::reactive({
    shiny::req(phe_df(), probs_obj(), K_chr(), cov_mx())
    shiny::withProgress(message = 'Effect scans ...', value = 0, {
      shiny::setProgress(1)
      doqtl2::listof_scan1coefCC(phe_df(), cov_mx(), probs_obj(), K_chr())
    })
  })
  output$effPlot <- shiny::renderPlot({
    shiny::req(input$pheno_name, scan_obj(), eff_obj())
    shiny::withProgress(message = 'Effect plots ...', value = 0, {
      shiny::setProgress(1)
      doqtl2::plot_eff(input$pheno_name, scan_obj(), eff_obj(), 
               input$scan_window)
    })
  })
  output$effSummary <- shiny::renderDataTable({
    shiny::req(eff_obj(), scan_obj())
    shiny::withProgress(message = 'Effect summary ...', value = 0, {
      shiny::setProgress(1)
      summary(eff_obj(), scan_obj())
    })
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  
  ## Effect and LOD Plot
  output$lod_effPlot <- shiny::renderPlot({
    shiny::req(win_par$chr_id, input$pheno_name, input$scan_window, 
        scan_obj(), eff_obj())
    eff_names <- names(eff_obj())
    shiny::withProgress(message = 'Effect & LOD plots ...', value = 0, {
      shiny::setProgress(1)
      par(mfrow=c(2,1))
      phenoi <- match(input$pheno_name, eff_names)
      print(doqtl2::show_peaks(win_par$chr_id, 
                 subset(scan_obj(), lodcolumn = phenoi),
                 mytitle="", 
                 xlim = input$scan_window))
      print(doqtl2::plot_eff(input$pheno_name, scan_obj(), eff_obj(), 
                     input$scan_window))
    })
  })
  
  output$pheno_choice <- shiny::renderUI({
    switch(shiny::req(input$button),
           "LOD & Effects" =,
           Effects = shiny::uiOutput(ns("pheno_name")))
  })
  output$win_choice <- shiny::renderUI({
    switch(shiny::req(input$button),
           LOD     =,
           "LOD & Effects" =,
           Effects = shiny::uiOutput(ns("scan_window")))
  })
  output$LOD <- shiny::renderUI({
    switch(shiny::req(input$button),
           LOD             =,
           "LOD & Effects" = shiny::plotOutput(ns("scanPlot")))
  })
  output$Effects <- shiny::renderUI({
    switch(shiny::req(input$button),
           Effects         =,
           "LOD & Effects" = shiny::plotOutput(ns("effPlot")))
  })
  output$Summary <- shiny::renderUI({
    switch(shiny::req(input$button),
           Summary = shiny::dataTableOutput(ns("effSummary")))
  })

  ## Downloads.
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("sum_effects_", win_par$chr_id, ".csv")) },
    content = function(file) {
      shiny::req(eff_obj(), scan_obj())
      write.csv(summary(eff_obj(), scan_obj()), file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("scan_", win_par$chr_id, ".pdf")) },
    content = function(file) {
      effs <- shiny::req(eff_obj())
      scans <- shiny::req(scan_obj())
      win <- shiny::req(input$scan_window)
      pdf(file, width=9,height=9)
      print(doqtl2::show_peaks(win_par$chr_id, scans, mytitle="", 
                       xlim=win))
      par(mfrow=c(2,1))
      for(phenoi in seq_along(effs)) {
        pheno <- names(effs)[phenoi]
        print(doqtl2::show_peaks(win_par$chr_id, 
                         subset(scans, lodcolumn = phenoi),
                         mytitle="", 
                         xlim = win))
        print(doqtl2::plot_eff(pheno, scans, effs, win))
      }
      dev.off()
    }
  )
  output$radio <- shiny::renderUI({
    shiny::radioButtons(ns("button"), "",
                 c("LOD","Effects","LOD & Effects","Summary"),
                 input$button)
  })
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyScan1Plot
#' @export
shinyScan1PlotUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::strong("Genome Scans"),
    shiny::uiOutput(ns("radio")),
    shiny::uiOutput(ns("pheno_choice")),
    shiny::uiOutput(ns("win_choice")),
    shiny::fluidRow(
      shiny::column(6, shiny::downloadButton(ns("downloadData"), "CSV")),
      shiny::column(6, shiny::downloadButton(ns("downloadPlot"), "Plots"))))
}
#' @rdname shinyScan1Plot
#' @export
shinyScan1PlotOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("LOD")),
    shiny::uiOutput(ns("Effects")),
    shiny::uiOutput(ns("Summary")))
}
