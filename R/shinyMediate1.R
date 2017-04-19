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
  
  ## Expression data
  expr_ls <- reactive({
    shiny::req(win_par, datapath())
    scan_window <- win_par$peak_Mbp + c(-1,1) * 2 ^ win_par$window_Mbp
    indID <- rownames(shiny::req(phe1_df()))
    
    # Get expression mRMNA measurements.
    out <- DOread::read_mrna(indID, win_par$chr_id,
                      scan_window[1], scan_window[2],
                      datapath())
    # Covariate matrix covar is global, but reget it here to be sure.
    out$cov_med <- readRDS(file.path(datapath(), "covar.rds"))[, c("sex", paste0("DOwave", 2:4))]
    out
  })
  
  
  ## Mediate1
  mediate_obj <- shiny::reactive({
    chr_id <- shiny::req(win_par$chr_id)
    shiny::req(phe1_df(), probs_obj(), K_chr(), cov_mx(),
               input$pos_Mbp, win_par$window_Mbp, expr_ls())
    shiny::withProgress(message = "Mediation Scan ...", value = 0, {
      shiny::setProgress(1)
      # qtl2pattern::mediate1
      med_test(chr_id, input$pos_Mbp, expr_ls(),
                            phe1_df(), cov_mx(), probs_obj()$probs, K_chr(), 
                            probs_obj()$map)
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
  ## Select phenotype for plots.
  output$med_plot <- shiny::renderUI({
    shiny::selectInput(ns("med_plot"), NULL,
                       choices = c("pos_lod","pos_pv","pv_lod"))
  })
  
  ## Mediate1 plot
  output$medPlot <- shiny::renderPlot({
    shiny::req(mediate_obj(), input$med_plot)
    shiny::withProgress(message = 'Mediation Plot ...', value = 0, {
      shiny::setProgress(1)
      plot(mediate_obj(), input$med_plot)
    })
  })
  ## Mediate1 plotly
  output$medPlotly <- plotly::renderPlotly({
    shiny::req(mediate_obj())
    shiny::withProgress(message = 'Mediation Plot ...', value = 0, {
      shiny::setProgress(1)
      plot(mediate_obj(), input$med_plot)
    })
  })

  output$medSummary <- shiny::renderDataTable({
    mediate_obj()
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))

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

  output$out_choice <- shiny::renderUI({
    switch(shiny::req(input$button),
           Static      = shiny::plotOutput(ns("medPlot")),
           Interactive = plotly::plotlyOutput(ns("medPlotly")),
           Summary     = shiny::dataTableOutput(ns("medSummary")))
  })
  output$radio <- shiny::renderUI({
    shiny::radioButtons(ns("button"), "",
                        c("Static","Interactive","Summary"),
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
                 input$pos_Mbp, win_par$window_Mbp)
      pdf(file, width=9,height=9)
      for(pheno in names(phe_df())) {
        med <- med_test(chr_id, input$pos_Mbp, expr_ls(),
                              phe_df()[, pheno, drop = FALSE],
                              cov_mx(), probs_obj()$probs, K_chr(), 
                              probs_obj()$map)
        print(plot(med), "pos_lod")
        print(plot(med), "pos_pv")
        print(plot(med), "pv_lod")
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
    shiny::uiOutput(ns("med_plot")),
    shiny::uiOutput(ns("pos_Mbp")),
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
