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
                       chr_id, phe_df, cov_mx,
                       pheno_anal, probs_obj, K_chr) {
  
  ## This is basically obsolete now.
  
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

  ## Select phenotype for plots.
  output$choose_pheno <- renderUI({
    req(pheno_anal())
    selectInput(ns("pheno_anal"), NULL,
                choices = names(pheno_anal()))
  })

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
