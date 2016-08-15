#' Shiny scan1 analysis and plot module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param plot_type,chr_id,phe_df,cov_mx,pheno_anal,probs_obj,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return object of class \code{\link[qtl2scan]{scan1}}
#'
#' @export
shinyScan1 <- function(input, output, session,
                       plot_type, chr_id, phe_df, cov_mx, pheno_anal,
                       probs_obj, K_chr) {
  ns <- session$ns

  ## Scan1
  scan_obj <- reactive({
    req(pheno_anal())
    withProgress(message = paste("Scan1 for", plot_type(), "..."),
                 value = 0,
    {
      setProgress(1)
      scan1(probs_obj(), phe_df(), K_chr(), cov_mx())
    })
  })

  # Scan Window slider
  output$choose_scan_window <- renderUI({
    req(plot_type(),pheno_anal())
    rng <- round(range(probs_obj()$map[[chr_id()]]), 2)
    sliderInput(ns("scan_window"), NULL, rng[1], rng[2],
                rng, step=.1)
  })

  snp_coef <- reactive({
    req(plot_type(),pheno_anal())
    phename <- dimnames(scan_obj()$lod)[[2]]
    coef_topsnp(scan_obj(),probs_obj(),phe_df(),K_chr(),cov_mx())
  })
  output$scanTable <- renderDataTable({
    snp_coef()
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))

  ## Select phenotype for plots.
  output$choose_pheno <- renderUI({
    req(pheno_anal())
    selectInput(ns("pheno_anal"), NULL,
                choices = names(pheno_anal()))
  })
  pheno_id <- reactive({
    pheno_anal <- req(input$pheno_anal)
    pheno_anal()[pheno_anal]
  })
  chr_pos <- reactive({
    chr <- req(chr_id())
    scan_w <- req(input$scan_window)
    scan_w <- round(scan_w, 2)
    paste(chr, scan_w[1], scan_w[2], sep = "_")
  })

  ## Module for scan1 plots.
  scan_window <- reactive({input$scan_window})
  callModule(shinyScan1Plot, "scan1plot", chr_id, phe_df, cov_mx,
             scan_window, chr_pos, pheno_id, scan_obj, probs_obj, K_chr)

  ## Module for SNP plots.
  callModule(shinyScan1SNP, "scan1snp", chr_id,
             scan_window, chr_pos, pheno_id, scan_obj, probs_obj)

  output$scan1_plot <- renderUI({
    req(pheno_anal())
    withProgress(message = paste("Plots for", plot_type(), "..."),
                 value = 0,
    {
      setProgress(1)
      if(plot_type() == "scans") {
        shinyScan1PlotUI(ns("scan1plot"))
      } else {
        shinyScan1SNPUI(ns("scan1snp"))
      }
    })
  })
  output$effectTable <- renderUI({
    req(pheno_anal(),plot_type())
    withProgress(message = paste("Effects for", plot_type(), "..."),
                 value = 0,
                 {
                   setProgress(1)
                   if(plot_type() == "scans") {
                     shinyScan1PlotOutput(ns("scan1plot"))
                   } else { "snsp"
                     dataTableOutput(ns("scanTable"))
                   }
                 })
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_effects_", chr_pos(), ".csv")) },
    content = function(file) {
      write.csv(snp_coef(), file)
    }
  )
  scan_obj
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
#' @rdname shinyScan1
#' @export
shinyScan1UI <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      tabPanel("summary",
               downloadButton(ns("downloadData"), "Download SNP Effects"),
               uiOutput(ns("effectTable"))),
      tabPanel("plots",
               uiOutput(ns("choose_scan_window")),
               uiOutput(ns("choose_pheno")),
               uiOutput(ns("scan1_plot")))
    )
  )
}
