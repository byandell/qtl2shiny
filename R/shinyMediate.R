#' Shiny coefficient analysis and plot module
#'
#' Shiny module for scan1 coefficient plots, with interfaces \code{shinyMediateUI} and  \code{shinyMediateOutput}.
#'
#' @param input,output,session standard shiny arguments
#' @param job_par,win_par,patterns,phe_mx,cov_df,probs_obj,K_chr,analyses_df,pmap_obj,covar,pheno_data,analyses_tbl,peaks,project_info,allele_info reactive arguments
#' @param id shiny identifier
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom CausalMST mediate1_test
#' @importFrom qtl2pattern covar_df_mx pheno_trans pheno_region expr_region
#' @importFrom qtl2 find_marker scan1
#' @importFrom ggplot2 autoplot
#' @importFrom shiny NS reactive req isTruthy
#'   radioButtons selectInput sliderInput updateSliderInput
#'   dataTableOutput plotOutput uiOutput
#'   renderDataTable renderPlot renderUI
#'   fluidRow column strong tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler callModule
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom dplyr filter
#' 
shinyMediate <- function(input, output, session,
                              job_par, win_par, patterns,
                              phe_mx, cov_df, probs_obj, K_chr, analyses_df,
                              pmap_obj, covar, pheno_data, analyses_tbl, peaks,
                              project_info, allele_info) {
  ns <- session$ns
  
  chr_id <- reactive({
    shiny::req(win_par$chr_id)
  })
  scan_window <- reactive({
    shiny::req(win_par)
    win_par$peak_Mbp + c(-1,1) * win_par$window_Mbp
  })
  ## Expression data
  expr_ls <- reactive({
    shiny::req(win_par)
    expr_region(chr_id(), scan_window(), covar(), 
                shiny::req(input$qtls), shiny::req(pmap_obj()),
                project_info())
  })
  
  ## Comediator data
  comed_ls <- reactive({
    shiny::req(input$pheno_name, win_par)
    comediator_region(input$pheno_name, chr_id(), scan_window(), 
                      covar(), analyses_tbl(), pheno_data(), peaks(), 
                      shiny::req(input$qtls), shiny::req(pmap_obj()))
  })
  med_ls <- reactive({
    out <- switch(shiny::req(input$med_type),
           expression = expr_ls(),
           phenotype = comediator_type(comed_ls(), shiny::req(peaks()),
                                       shiny::isTruthy(input$other)))
    out
  })
  
  # Get genotype matrix and map at 
  peak_mar <- reactive({
    qtl2::find_marker(probs_obj()$map, chr_id(), input$pos_Mbp)
  })
  geno_max <- reactive({
    shiny::req(input$pos_Mbp, probs_obj())
    subset(probs_obj()$probs, chr = chr_id(), mar = peak_mar())[[1]][,,1]
  })
  
  ## Scatter Plots
  shiny::callModule(shinyScatter, "scatter",
                    input, patterns, 
                    geno_max, peak_mar, med_ls, mediate_signif,
                    phe1_mx, cov_df, K_chr,
                    allele_info)

  ## Mediate1
  probs_chr <- reactive({
    probs_obj()$probs[[chr_id()]]
  })
  mediate_obj <- shiny::reactive({
    shiny::req(phe1_mx(), probs_obj(), K_chr(), cov_df(), geno_max(), 
               input$pos_Mbp, input$med_type, med_ls())
    shiny::withProgress(message = "Mediation Scan ...", value = 0, {
      shiny::setProgress(1)
      med_test(med_ls(), geno_max(), phe1_mx(), K_chr(), cov_df(),
               input$pos_Mbp, data_type = input$med_type,
               probs_chr())
    })
  })
  mediate_signif <- shiny::reactive({
    out <- shiny::req(mediate_obj())
    out$best <- dplyr::filter(out$best, pvalue <= 0.1)
    class(out) <- class(mediate_obj())
    out
  })

  phe1_mx <- reactive({
    shiny::req(phe_mx())
    phename <- shiny::req(input$pheno_name)
    if(phename %in% colnames(phe_mx())) {
      phe_mx()[, phename, drop = FALSE]
    } else {
      NULL
    }
  })
  ## Select phenotype for plots.
  output$pheno_name <- shiny::renderUI({
    shiny::req(phe_mx())
    shiny::selectInput(ns("pheno_name"), NULL,
                choices = colnames(phe_mx()),
                selected = input$pheno_name)
  })
  ## Select plot format.
  output$med_plot <- shiny::renderUI({
    shiny::selectInput(ns("med_plot"), NULL,
                       choices = c("Position by LOD", 
                                   "Position by P-value", 
                                   "P-value by LOD",
                                   "Allele Effects",
                                   "Mediator Effects"),
                       selected = input$med_plot)
  })
  ## Select type of mediation.
  output$med_type <- shiny::renderUI({
    shiny::selectInput(ns("med_type"), NULL,
                       choices = c("phenotype","expression"),
                       selected = input$med_type)
  })
  
  med_plot_type <- reactive({
    switch(shiny::req(input$med_plot),
           "Position by LOD" = "pos_lod",
           "Position by P-value" = "pos_pvalue",
           "P-value by LOD" = "pvalue_lod",
           "Allele Effects" = "alleles",
           "Mediator Effects" = "mediator")
  })
  ## Mediate1 plot
  output$medPlot <- shiny::renderPlot({
    if(!shiny::isTruthy(med_ls()) || !shiny::isTruthy(mediate_obj())) {
      plot_null("too much\nmissing data\nin mediators\nreduce window width")
    } else {
      shiny::req(med_plot_type())
      shiny::withProgress(message = 'Mediation Plot ...', value = 0, {
        shiny::setProgress(1)
        ggplot2::autoplot(
          mediate_obj(), med_plot_type(),
          local_only = input$local, 
          significant = input$signif)
      })
    }
  })
  ## Mediate1 plotly
  output$medPlotly <- plotly::renderPlotly({
    shiny::req(mediate_obj())
    shiny::withProgress(message = 'Mediation Plotly ...', value = 0, {
      shiny::setProgress(1)
      ggplot2::autoplot(
        mediate_signif(), med_plot_type(),
        local_only = input$local, 
        significant = TRUE)
    })
  })

  output$medSummary <- shiny::renderDataTable({
    shiny::req(mediate_obj())$best
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
  output$qtls <- shiny::renderUI({
    if(is.null(selected <- input$qtls))
      selected <- 2
    shiny::radioButtons(ns("qtls"), "",
                        c("1 QTL" = 1, "2 QTLs" = 2),
                        selected, inline = TRUE)
  })
  output$radio <- shiny::renderUI({
    shiny::radioButtons(ns("button"), "",
                        c("Static","Interactive","Summary"),
                        "Static")
  })
  output$checkplot <- shiny::renderUI({
    shiny::checkboxInput(ns("checkplot"), "Scatter Plot", input$checkplot)
  })
  output$local_other <- shiny::renderUI({
    switch(shiny::req(input$med_type),
           expression = shiny::checkboxInput(ns("local"), "Local?", input$local),
           phenotype  = shiny::checkboxInput(ns("other"), "Other types?", input$other))
  })
  output$signif <- shiny::renderUI({
    if(shiny::isTruthy(input$signif)) {
      value <- input$signif
    } else {
      value <- TRUE
    }
    shiny::checkboxInput(ns("signif"), "Significant?", value)
  })

  ## Downloads.
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("mediate_", chr_id(), "_", win_par$peak_Mbp, ".csv")) },
    content = function(file) {
      shiny::req(mediate_obj())
      write.csv(mediate_obj()$best, file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("mediate_", chr_id(), "_", win_par$peak_Mbp, ".pdf")) },
    content = function(file) {
      shiny::req(phe_mx(), geno_max(), K_chr(), cov_df(),
                 input$pos_Mbp, input$med_type)
      pdf(file, width=9,height=9)
      for(pheno in colnames(phe_mx())) {
        med <- med_test(med_ls(), geno_max(),
                        phe_mx()[, pheno, drop = FALSE],
                        K_chr(), cov_df(), input$pos_Mbp,
                        data_type = input$med_type,
                        probs_chr())
        print(ggplot2::autoplot(
          med, "pos_lod",
          local_only = input$local, 
          significant = input$signif))
        print(ggplot2::autoplot(
          med, "pos_pvalue",
          local_only = input$local, 
          significant = TRUE))
        print(ggplot2::autoplot(
          med, "pvalue_lod",
          local_only = input$local, 
          significant = TRUE))
        print(ggplot2::autoplot(
          med, "mediator",
          local_only = input$local, 
          significant = TRUE))
      }
      dev.off()
    })
  output$mediation <- renderUI({
    shiny::tagList(
      shiny::uiOutput(ns("qtls")),
      shiny::uiOutput(ns("radio")),
      shiny::uiOutput(ns("pheno_name")),
      shiny::uiOutput(ns("med_type")),
      shiny::fluidRow(
        shiny::column(6, shiny::uiOutput(ns("signif"))),
        shiny::column(6, shiny::uiOutput(ns("local_other")))),
      shiny::uiOutput(ns("med_plot")),
      shiny::uiOutput(ns("pos_Mbp")),
      shiny::fluidRow(
        shiny::column(6, shiny::downloadButton(ns("downloadData"), "CSV")),
        shiny::column(6, shiny::downloadButton(ns("downloadPlot"), "Plots"))))
  })
  output$medUI <- shiny::renderUI({
    switch(1 + shiny::isTruthy(input$checkplot),
           {
             shiny::uiOutput(ns("mediation"))
           },
           {
             shinyScatterUI(ns("scatter"))
           })
  })
  output$medOutput <- shiny::renderUI({
    if(shiny::isTruthy(input$checkplot))
      shinyScatterOutput(ns("scatter"))
    else
      shiny::uiOutput(ns("out_choice"))
  })
}
shinyMediateUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::strong("Mediation"),
    shiny::uiOutput(ns("checkplot")),
    shiny::uiOutput(ns("medUI")))
}
shinyMediateOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("medOutput"))
}
