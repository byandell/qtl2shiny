#' Shiny coefficient analysis and plot module
#'
#' Shiny module for scan1 coefficient plots.
#'
#' @param input,output,session standard shiny arguments
#' @param chr_id,phe_df,cov_mx,pheno_anal,probs_obj,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return 2-element vector of scan window
#'
#' @export
shinyScan1Plot <- function(input, output, session,
                       chr_id, phe_df, cov_mx, 
                       pheno_anal, probs_obj, K_chr) {
  ns <- session$ns
  
  ## Scan1
  scan_obj <- reactive({
    req(pheno_anal())
    withProgress(message = "Genome Scan ...", value = 0, {
      setProgress(1)
      scan1(probs_obj(), phe_df(), K_chr(), cov_mx())
    })
  })
  
  # Scan Window slider
  output$scan_window <- renderUI({
    req(chr_id(), pheno_anal(), probs_obj())
    rng <- round(range(probs_obj()$map[[chr_id()]]), 2)
    sliderInput(ns("scan_window"), NULL, rng[1], rng[2],
                rng, step=.1)
  })
  
  ## Select phenotype for plots.
  output$pheno_anal <- renderUI({
    req(pheno_anal())
    selectInput(ns("pheno_anal"), NULL,
                choices = names(pheno_anal()))
  })
  pheno_id <- reactive({
    pheno_anal()[req(input$pheno_anal)]
  })
  
  ## Scan1 plot
  output$scanPlot <- renderPlot({
    withProgress(message = 'Genome LOD Plot ...', value = 0, {
      setProgress(1)
      show_peaks(chr_id(), scan_obj(), mytitle="", 
                 xlim=input$scan_window)
    })
  })

  ## Coefficient Effects.
  eff_obj <- reactive({
    withProgress(message = 'Effect scans ...', value = 0, {
      setProgress(1)
      listof_scan1coefCC(phe_df(), cov_mx(), probs_obj(), K_chr())
    })
  })
  output$effPlot <- renderPlot({
    req(pheno_id(),eff_obj())
    withProgress(message = 'Effect plots ...', value = 0, {
      setProgress(1)
      plot_eff(pheno_id(), scan_obj(), eff_obj(), input$scan_window)
    })
  })
  output$effSummary <- renderDataTable({
    req(eff_obj(),scan_obj())
    withProgress(message = 'Effect summary ...', value = 0, {
      setProgress(1)
      summary(eff_obj(),scan_obj())
    })
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  
  ## Downloads.
  chr_pos <- reactive({
    make_chr_pos(chr_id(), range = input$scan_window)
  })
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
      show_peaks(chr_id(), scan_obj(), mytitle="", xlim=input$scan_window)
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
    fluidRow(
      column(3, h4(strong("Genome Plots"))),
      column(6, uiOutput(ns("scan_window"))),
      column(3, 
             downloadButton(ns("downloadData"), "CSV"),
             downloadButton(ns("downloadPlot"), "Plots"))),
    tabsetPanel(
      tabPanel("LOD", 
               plotOutput(ns("scanPlot"))),
      tabPanel("Effects",
               plotOutput(ns("effPlot")),
               uiOutput(ns("pheno_anal"))),
      tabPanel("summary",
               dataTableOutput(ns("effSummary"))))
  )
}
