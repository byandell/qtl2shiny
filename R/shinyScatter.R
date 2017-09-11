#' Shiny scatter plot module
#'
#' Shiny module for scatter plots.
#'
#' @param input,output,session standard shiny arguments
#' @param med_par,patterns,geno_max,peak_mar,med_ls,mediate_obj,phe_df,cov_mx,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom qtl2pattern mediate1
#' @importFrom qtl2scan scan1
#' @importFrom qtl2ggplot plot_scan1
#' @importFrom shiny NS reactive req isTruthy
#'   radioButtons selectInput sliderInput updateSliderInput
#'   dataTableOutput plotOutput uiOutput
#'   renderDataTable renderPlot renderUI
#'   fluidRow column strong tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
#' @importFrom plotly renderPlotly plotlyOutput
#' 
shinyScatterPlot <- function(input, output, session,
                  med_par, patterns, 
                  geno_max, peak_mar, med_ls, mediate_obj,
                  phe_df, cov_mx, K_chr) {
  ns <- session$ns

  ## Select triad for plots.
  output$triad <- shiny::renderUI({
    triad <- shiny::req(mediate_obj())$best$triad
    choices <- levels(triad)
    choices <- choices[choices %in% unique(triad)]
    shiny::selectInput(ns("triad"), NULL,
                       choices = choices)
  })
  medID <- reactive({
    ifelse("symbol" %in% names(shiny::req(mediate_obj())$best), "symbol", "longname")
  })
  ## Select mediator for plots.
  output$med_name <- shiny::renderUI({
    shiny::req(mediate_obj(), input$triad)
    choices <- dplyr::filter(mediate_obj()$best, triad == input$triad)[[medID()]]
    shiny::selectInput(ns("med_name"), NULL,
                       choices = choices)
  })
  # Select pattern 
  sdps <- shiny::reactive({
    shiny::req(patterns(), med_par$pheno_name)
    unique(dplyr::filter(patterns(), pheno == med_par$pheno_name)$sdp)
  })
  output$pattern <- shiny::renderUI({
    choices <- CCSanger::sdp_to_pattern(sdps())
     shiny::selectInput(ns("pattern"), NULL,
                       choices = choices, input$pattern)
  })
  
  scat_dat <- reactive({
    shiny::req(geno_max(), phe_df(), med_ls(), medID(),
               input$med_name, input$pattern)
    med_scat(med_ls(), geno_max(), phe_df(), K_chr()[[1]], cov_mx(), 
             input$pattern, input$med_name, medID())
  })
  
  ## Select plot format.
  output$med_plot <- shiny::renderUI({
    shiny::selectInput(ns("med_plot"), NULL,
                       choices = c("by_mediator", 
                                   "by_target", 
                                   "driver_offset", 
                                   "driver"))
  })

  ## Scatter plot
  output$scatPlot <- shiny::renderPlot({
    if(!shiny::isTruthy(patterns())) {
      return(plot_null("first run\nAllele Patterns"))
    }
    if(!shiny::isTruthy(scat_dat())) {
      plot_null("too much\nmissing data\nin mediators\nreduce window width")
    } else {
      shiny::req(input$med_plot, input$med_name, phe_df())
      shiny::withProgress(message = 'Scatter Plot ...', value = 0, {
        shiny::setProgress(1)
        plot(scat_dat(), type = input$med_plot,
             dname = peak_mar(),
             mname = input$med_name,
             tname = names(phe_df()))
      })
    }
  })

  ## Downloads.
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("scatter.csv")) },
    content = function(file) {
      shiny::req(mediate_obj())
      write.csv(mediate_obj()$best, file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("scatter.pdf")) },
    content = function(file) {
      shiny::req(phe_df(), scat_dat(), input$med_plot)
      pdf(file, width=9,height=9)
      for(types in c("by_mediator", 
                     "by_target", 
                     "driver_offset", 
                     "driver")) {
        print(plot(scat_dat(), type = types,
                   dname = peak_mar(),
                   mname = input$med_name,
                   tname = names(phe_df())))
      }
      dev.off()
    })
  med_ls
}
#' @param id identifier for shiny use
#' @rdname shinyScatterPlot
#' @export
shinyScatterPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("triad")),
    shiny::uiOutput(ns("med_name")),
    shiny::uiOutput(ns("pattern")),
    shiny::uiOutput(ns("med_plot")),
    shiny::fluidRow(
      shiny::column(6, shiny::downloadButton(ns("downloadData"), "CSV")),
      shiny::column(6, shiny::downloadButton(ns("downloadPlot"), "Plots"))))
}
#' @rdname shinyScatterPlot
#' @export
shinyScatterPlotOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::plotOutput(ns("scatPlot"))
}
