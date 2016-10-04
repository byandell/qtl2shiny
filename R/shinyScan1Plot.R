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
shinyScan1Plot <- function(input, output, session,
                  win_par, pmap_obj, 
                  phe_df, cov_mx, probs_obj, K_chr) {
  ns <- session$ns
  
  ## Scan1
  scan_obj <- reactive({
    req(phe_df(), probs_obj(), K_chr(), cov_mx())
    withProgress(message = "Genome Scan ...", value = 0, {
      setProgress(1)
      scan1(probs_obj(), phe_df(), K_chr(), cov_mx())
    })
  })
  
  # Scan Window slider
  output$scan_window <- renderUI({
    chr_id <- req(win_par$chr_id)
    map <- req(pmap_obj())[[chr_id]]
    rng <- round(2 * range(map)) / 2
    if(is.null(selected <- input$scan_window))
      selected <- rng
    sliderInput(ns("scan_window"), NULL, rng[1], rng[2],
                selected, step=.5)
  })
  ## Reset scan_window if chromosome changes.
  observeEvent(win_par$chr_id, {
    map <- req(pmap_obj())
    chr <- req(win_par$chr_id)
    rng <- round(2 * range(map[[chr]])) / 2
    updateSliderInput(session, "scan_window", NULL, rng, 
                      rng[1], rng[2], step=.5)
  })
  
  ## Select phenotype for plots.
  output$pheno_name <- renderUI({
    req(phe_df())
    selectInput(ns("pheno_name"), NULL,
                choices = names(phe_df()))
  })

  ## Scan1 plot
  output$scanPlot <- renderPlot({
    req(win_par$chr_id, input$scan_window, scan_obj())
    withProgress(message = 'Genome LOD Plot ...', value = 0, {
      setProgress(1)
      show_peaks(win_par$chr_id, scan_obj(), mytitle="", 
                 xlim=input$scan_window)
    })
  })

  ## Coefficient Effects.
  eff_obj <- reactive({
    req(phe_df(), probs_obj(), K_chr(), cov_mx())
    withProgress(message = 'Effect scans ...', value = 0, {
      setProgress(1)
      listof_scan1coefCC(phe_df(), cov_mx(), probs_obj(), K_chr())
    })
  })
  output$effPlot <- renderPlot({
    req(input$pheno_name,eff_obj())
    withProgress(message = 'Effect plots ...', value = 0, {
      setProgress(1)
      plot_eff(input$pheno_name, scan_obj(), eff_obj(), input$scan_window)
    })
  })
  output$effSummary <- renderDataTable({
    req(eff_obj(), scan_obj())
    withProgress(message = 'Effect summary ...', value = 0, {
      setProgress(1)
      summary(eff_obj(), scan_obj())
    })
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))

  output$scan_choice <- renderUI({
    switch(req(input$button),
           Effects = uiOutput(ns("pheno_name")))
  })
  output$win_choice <- renderUI({
    switch(req(input$button),
           LOD     =,
           Effects = uiOutput(ns("scan_window")))
  })
  output$scan_output <- renderUI({
    switch(req(input$button),
           LOD     = plotOutput(ns("scanPlot")),
           Effects = plotOutput(ns("effPlot")),
           Summary = dataTableOutput(ns("effSummary")))
  })

  ## Downloads.
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("sum_effects_", win_par$chr_id, ".csv")) },
    content = function(file) {
      req(eff_obj(),scan_obj())
      write.csv(summary(eff_obj(),scan_obj()), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("scan_", win_par$chr_id, ".pdf")) },
    content = function(file) {
      req(eff_obj())
      pdf(file)
      show_peaks(win_par$chr_id, scan_obj(), mytitle="", 
                 xlim=req(input$scan_window))
      for(pheno in names(eff_obj()))
        plot_eff(pheno)
      dev.off()
    }
  )
  output$radio <- renderUI({
    radioButtons(ns("button"), "",
                 c("LOD","Effects","Summary"),
                 input$button)
  })
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyScan1Plot
#' @export
shinyScan1PlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4(strong("Genome Scans")),
    uiOutput(ns("radio")),
    uiOutput(ns("scan_choice")),
    uiOutput(ns("win_choice")),
    fluidRow(
      column(6, downloadButton(ns("downloadData"), "CSV")),
      column(6, downloadButton(ns("downloadPlot"), "Plots"))))
}
#' @rdname shinyScan1Plot
#' @export
shinyScan1PlotOutput <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("scan_output")))
}
