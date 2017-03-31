#' Shiny coefficient analysis and plot module
#'
#' Shiny module for scan1 coefficient plots.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,probs_obj,K_chr,analyses_df,datapath reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom qtl2pattern mediate1
#' @importFrom qtl2scan scan1
#' @importFrom qtl2ggplot plot_scan1
#' @importFrom shiny NS reactive req 
#'   radioButtons selectInput sliderInput updateSliderInput
#'   dataTableOutput plotOutput uiOutput
#'   renderDataTable renderPlot renderUI
#'   fluidRow column strong tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
#' @importFrom plotly renderPlotly plotlyOutput
#' 
shinyMediate1Plot <- function(input, output, session,
                  win_par,  
                  phe_df, cov_mx, probs_obj, K_chr, analyses_df,
                  datapath) {
  ns <- session$ns
  
  ## Mediate1
  mediate_obj <- shiny::reactive({
    chr_id <- shiny::req(win_par$chr_id)
    shiny::req(phe1_df(), probs_obj(), K_chr(), cov_mx(),
               win_par$peak_Mbp, win_par$window_Mbp)
    shiny::withProgress(message = "Mediation Scan ...", value = 0, {
      shiny::setProgress(1)
      qtl2pattern::mediate1(chr_id, win_par$peak_Mbp, 2^win_par$window_Mbp,
                            phe1_df(), cov_mx(), probs_obj()$probs, K_chr(), 
                            probs_obj()$map, datapath())
    })
  })
  
  phe1_df <- reactive({
    phe_df()[, shiny::req(input$pheno_name), drop = FALSE]
  })
  ## Select phenotype for plots.
  output$pheno_name <- shiny::renderUI({
    shiny::req(phe_df())
    shiny::selectInput(ns("pheno_name"), NULL,
                choices = names(phe_df()))
  })

  ## Mediate1 plot
  output$medPlot <- shiny::renderPlot({
    shiny::req(mediate_obj())
    shiny::withProgress(message = 'Mediation Plot ...', value = 0, {
      shiny::setProgress(1)
      plot(mediate_obj())
    })
  })
  ## Mediate1 plotly
  output$medPlotly <- plotly::renderPlotly({
    shiny::req(mediate_obj())
    shiny::withProgress(message = 'Mediation Plot ...', value = 0, {
      shiny::setProgress(1)
      plot(mediate_obj())
    })
  })

  output$medSummary <- shiny::renderDataTable({
    mediate_obj()
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))

  output$out_choice <- shiny::renderUI({
    switch(shiny::req(input$button),
           Plot        = shiny::plotOutput(ns("medPlot")),
           Interactive = plotly::plotlyOutput(ns("medPlotly")),
           Summary     = shiny::dataTableOutput(ns("medSummary")))
  })
  output$radio <- shiny::renderUI({
    shiny::radioButtons(ns("button"), "",
                        c("Plot","Interactive","Summary"),
                        input$button)
  })
  
  ## Downloads.
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("mediate_", win_par$chr_id, "_", win_par$peak_Mbp, ".csv")) },
    content = function(file) {
      shiny::req(mediate_obj())
      write.csv(mediate_obj(), file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("mediate_", win_par$chr_id, "_", win_par$peak_Mbp, ".pdf")) },
    content = function(file) {
      chr_id <- shiny::req(win_par$chr_id)
      shiny::req(phe_df(), probs_obj(), K_chr(), cov_mx(),
                 win_par$peak_Mbp, win_par$window_Mbp)
      pdf(file, width=9,height=9)
      for(pheno in names(phe_df())) {
        med <- qtl2pattern::mediate1(chr_id, win_par$peak_Mbp, 2^win_par$window_Mbp,
                              phe_df()[, pheno, drop = FALSE],
                              cov_mx(), probs_obj()$probs, K_chr(), 
                              probs_obj()$map, datapath())
        print(plot(med))
      }
      dev.off()
    }
  )
}
#' @param id identifier for shiny use
#' @rdname shinyMediate1Plot
#' @export
shinyMediate1PlotUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::strong("Mediation"),
    shiny::uiOutput(ns("radio")),
    shiny::uiOutput(ns("pheno_name")),
    shiny::fluidRow(
      shiny::column(6, shiny::downloadButton(ns("downloadData"), "CSV")),
      shiny::column(6, shiny::downloadButton(ns("downloadPlot"), "Plots"))))
}
#' @rdname shinyMediate1Plot
#' @export
shinyMediate1PlotOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("out_choice"))
}
