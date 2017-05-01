#' Shiny scatter plot module
#'
#' Shiny module for scatter plots.
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
                  med_par,  
                  geno_max, med_ls, mediate_obj,
                  pattern,
                  phe_df, cov_mx, K_chr) {
  ns <- session$ns
  
  ## Select triad for plots.
  output$triad <- shiny::renderUI({
    shiny::req(mediate_obj())
    choices <- levels(mediate_obj()$triad)
    choices <- choices[choices %in% unique(mediate_obj()$triad)]
    shiny::selectInput(ns("triad"), NULL,
                       choices = choices)
  })
  medID <- reactive({
    ifelse("symbol" %in% names(shiny::req(mediate_obj())), "symbol", "longname")
  })
  ## Select mediator for plots.
  output$med_name <- shiny::renderUI({
    shiny::req(mediate_obj(), input$triad)
    choices <- dplyr::filter(med_ls(), triad == input$triad)[[medID()]]
    shiny::selectInput(ns("med_name"), NULL,
                       choices = choices)
  })
  # Select pattern 
  sdps <- shiny::reactive({
    shiny::req(pattern(), med_par$pheno_name)
    unique(dplyr::filter(pattern(), pheno == med_par$pheno_name)$sdp)
  })
  output$pattern <- shiny::renderUI({
    shiny::selectInput(ns("pattern"), NULL,
                       choices = CCSanger::sdp_to_pattern(sdps()))
  })
  
  scat_dat <- reactive({
    shiny::req(geno_max(), phe_df(), med_ls(), input$med_name, input$pattern)
    sdp <- sdps()[CCSanger::sdp_to_pattern(sdps()) == input$pattern]
    CausalMST:::med_scatter(geno_max(), phe_df(), med_ls()[[1]][, input$med_name, drop = FALSE],
                            K_chr[[1]], cov_mx(), med_ls()$cov_med,
                            qtl2scan::fit1,
                            sdp = sdp, allele = TRUE)
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
  output$medPlot <- shiny::renderPlot({
    if(!shiny::isTruthy(scat_dat())) {
      plot_null("too much\nmissing data\nin mediators\nreduce window width")
    } else {
      shiny::req(input$med_plot)
      shiny::withProgress(message = 'Scatter Plot ...', value = 0, {
        shiny::setProgress(1)
        plot(scat_dat(), type = input$med_plot)
      })
    }
  })

  ## Downloads.
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("scatter_", win_par$chr_id, "_", win_par$peak_Mbp, ".csv")) },
    content = function(file) {
      shiny::req(mediate_obj())
      write.csv(mediate_obj(), file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("scatter_", win_par$chr_id, "_", win_par$peak_Mbp, ".pdf")) },
    content = function(file) {
      chr_id <- shiny::req(win_par$chr_id)
      shiny::req(phe_df(), probs_obj(), K_chr(), cov_mx(),
                 input$pos_Mbp, win_par$window_Mbp)
      pdf(file, width=9,height=9)
      for(pheno in names(phe_df())) {
        print(plot(scat_dat(), type = input$med_plot))
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
    shiny::strong("Scatter Plot"),
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
  shiny::uiOutput(ns("medPlot"))
}
