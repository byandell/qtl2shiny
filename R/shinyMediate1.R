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
#' @importFrom shiny NS reactive req isTruthy
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
  
  chr_id <- reactive({
    shiny::req(win_par$chr_id)
  })
  scan_window <- reactive({
    shiny::req(win_par)
    win_par$peak_Mbp + c(-1,1) * 2 ^ win_par$window_Mbp
  })
  ## Expression data
  expr_ls <- reactive({
    shiny::req(win_par)
    # Covariate matrix covar is global.
    expr_region(chr_id(), scan_window(), datapath(), covar)
  })
  
  ## Comediator data
  comed_ls <- reactive({
    shiny::req(input$pheno_name, win_par)
    # Objects covar, analyses_tbl, pheno_data, peaks are global.
    comediator_region(input$pheno_name, chr_id(), scan_window(), 
                      covar, analyses_tbl, pheno_data, peaks)
  })
  med_ls <- reactive({
    out <- switch(shiny::req(input$med_type),
           expression = expr_ls(),
           phenotype = comediator_type(comed_ls(), shiny::isTruthy(input$other)))
    out
  })
  
  # Get genotype matrix and map at 
  geno_max <- reactive({
    shiny::req(input$pos_Mbp, probs_obj())
    peak_mar <- qtl2geno::find_marker(probs_obj()$map, chr_id(), input$pos_Mbp)
    subset(probs_obj()$probs, chr = chr_id(), mar = peak_mar)[[1]][,,1]
  })
  
  
  ## Scatter Plots
  shiny::callModule(shinyScatterPlot, "scatter",
                    win_par,
                    phe1_df, cov_mx, probs_obj, K_chr, analyses_df,
                    data_path)
  
  ## Mediate1
  mediate_obj <- shiny::reactive({
    shiny::req(phe1_df(), probs_obj(), K_chr(), cov_mx(), geno_max(), 
               input$pos_Mbp, input$med_type, med_ls())
    shiny::withProgress(message = "Mediation Scan ...", value = 0, {
      shiny::setProgress(1)
      # qtl2pattern::mediate1
      med_test(med_ls(), geno_max(), phe1_df(), K_chr(), cov_mx(),
               input$pos_Mbp, data_type = input$med_type)
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
  ## Select plot format.
  output$med_plot <- shiny::renderUI({
    shiny::selectInput(ns("med_plot"), NULL,
                       choices = c("Position by LOD", 
                                   "Position by P-value", 
                                   "P-value by LOD"))
  })
  ## Select type of mediation.
  output$med_type <- shiny::renderUI({
    shiny::selectInput(ns("med_type"), NULL,
                       choices = c("expression","phenotype"))
  })
  
  med_plot_type <- reactive({
    switch(shiny::req(input$med_plot),
           "Position by LOD" = "pos_lod",
           "Position by P-value" = "pos_pvalue",
           "P-value by LOD" = "pvalue_lod")
  })
  ## Mediate1 plot
  output$medPlot <- shiny::renderPlot({
    if(!shiny::isTruthy(mediate_obj())) {
      plot_null("too much\nmissing data\nin mediators\nreduce window width")
    } else {
      shiny::withProgress(message = 'Mediation Plot ...', value = 0, {
        shiny::setProgress(1)
        plot(mediate_obj(), med_plot_type(),
             local_only = input$local, 
             significant = input$signif)
      })
    }
  })
  ## Mediate1 plotly
  output$medPlotly <- plotly::renderPlotly({
    shiny::req(mediate_obj())
    shiny::withProgress(message = 'Mediation Plot ...', value = 0, {
      shiny::setProgress(1)
      plot(mediate_obj(), med_plot_type(),
           local_only = input$local, 
           significant = TRUE)
    })
  })

  output$medSummary <- shiny::renderDataTable({
    mediate_obj()
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))

  # Scan Window slider
  output$pos_Mbp <- shiny::renderUI({
    map <- shiny::req(probs_obj())$map[[chr_id()]]
    rng <- round(2 * range(map)) / 2
    if(is.null(selected <- input$pos_Mbp))
      selected <- req(win_par$peak_Mbp)
    shiny::sliderInput(ns("pos_Mbp"), NULL, rng[1], rng[2],
                       selected, step=.1)
  })
  ## Reset pos_Mbp if chromosome changes.
  observeEvent(chr_id(), {
    map <- shiny::req(probs_obj()$map)
    rng <- round(2 * range(map[[chr_id()]])) / 2
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
  output$local_other <- shiny::renderUI({
    switch(shiny::req(input$med_type),
           expression = shiny::checkboxInput(ns("local"), "Local?"),
           phenotype  = shiny::checkboxInput(ns("other"), "Other types?"))
  })
  
  
  ## Downloads.
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("mediate_", chr_id(), "_", win_par$peak_Mbp, ".csv")) },
    content = function(file) {
      shiny::req(mediate_obj())
      write.csv(mediate_obj(), file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("mediate_", chr_id(), "_", win_par$peak_Mbp, ".pdf")) },
    content = function(file) {
      shiny::req(phe_df(), geno_max(), K_chr(), cov_mx(),
                 input$pos_Mbp, input$med_type)
      pdf(file, width=9,height=9)
      for(pheno in names(phe_df())) {
        med <- med_test(med_ls(), geno_max(),
                        phe_df()[, pheno, drop = FALSE],
                        K_chr(), cov_mx(), input$pos_Mbp,
                        data_type = input$med_type)
        print(plot(med), "pos_lod",
              local_only = input$local, 
              significant = input$signif)
        print(plot(med), "pos_pvalue",
              local_only = input$local, 
              significant = input$signif)
        print(plot(med), "pvalue_lod",
              local_only = input$local, 
              significant = input$signif)
      }
      dev.off()
    })
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
    shiny::uiOutput(ns("med_type")),
    shiny::fluidRow(
      shiny::column(6, shiny::checkboxInput(ns("signif"), "Significant?")),
      shiny::column(6, shiny::uiOutput(ns("local_other")))),
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
