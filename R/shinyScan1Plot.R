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
    req(input$pheno_name, scan_obj(), eff_obj())
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
  
  ## Effect and LOD Plot
  output$lod_effPlot <- renderPlot({
    req(win_par$chr_id, input$pheno_name, input$scan_window, 
        scan_obj(), eff_obj())
    eff_names <- names(eff_obj())
    withProgress(message = 'Effect & LOD plots ...', value = 0, {
      setProgress(1)
      par(mfrow=c(2,1))
      phenoi <- match(input$pheno_name, eff_names)
      show_peaks(win_par$chr_id, 
                 subset(scan_obj(), lodcolumn = phenoi),
                 mytitle="", 
                 xlim = input$scan_window)
      plot_eff(input$pheno_name, scan_obj(), eff_obj(), input$scan_window)
    })
  })
  
  output$pheno_choice <- renderUI({
    switch(req(input$button),
           "LOD & Effects" =,
           Effects = uiOutput(ns("pheno_name")))
  })
  output$win_choice <- renderUI({
    switch(req(input$button),
           LOD     =,
           "LOD & Effects" =,
           Effects = uiOutput(ns("scan_window")))
  })
  output$LOD <- renderUI({
    switch(req(input$button),
           LOD             =,
           "LOD & Effects" = plotOutput(ns("scanPlot")))
  })
  output$Effects <- renderUI({
    switch(req(input$button),
           Effects         =,
           "LOD & Effects" = plotOutput(ns("effPlot")))
  })
  output$Summary <- renderUI({
    switch(req(input$button),
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
      effs <- req(eff_obj())
      scans <- req(scan_obj())
      win <- req(input$scan_window)
      pdf(file, width=9,height=9)
      show_peaks(win_par$chr_id, scans, mytitle="", 
                 xlim=win)
      par(mfrow=c(2,1))
      for(phenoi in seq_along(effs)) {
        pheno <- names(effs)[phenoi]
        show_peaks(win_par$chr_id, 
                   subset(scans, lodcolumn = phenoi),
                   mytitle="", 
                   xlim = win)
        plot_eff(pheno, scans, effs, win)
      }
      dev.off()
    }
  )
  output$radio <- renderUI({
    radioButtons(ns("button"), "",
                 c("LOD","Effects","LOD & Effects","Summary"),
                 input$button)
  })
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyScan1Plot
#' @export
shinyScan1PlotUI <- function(id) {
  ns <- NS(id)
  tagList(
    strong("Genome Scans"),
    uiOutput(ns("radio")),
    uiOutput(ns("pheno_choice")),
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
    uiOutput(ns("LOD")),
    uiOutput(ns("Effects")),
    uiOutput(ns("Summary")))
}
