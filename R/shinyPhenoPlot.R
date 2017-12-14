#' Shiny Phenotype Plot module
#'
#' Shiny module to plot phenotypes.
#'
#' @param input,output,session standard shiny arguments
#' @param phe_mx,cov_df reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return 2-element vector of scan window
#'
#' @export
#' @importFrom shiny NS 
#'   plotOutput tableOutput
#'   renderPlot renderTable
#'   tagList
#'   withProgress setProgress
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes_string geom_density geom_rug
#' @importFrom GGally ggscatmat
shinyPhenoPlot <- function(input, output, session,
                           phe_mx, cov_df) {
  ns <- session$ns
  
  ## Scatter plot or density
  output$phe_sum <- shiny::renderTable({
    shiny::req(phe_mx())
    shiny::withProgress(message = 'Pheno Summary ...', value = 0, {
      shiny::setProgress(1)
      summary_na(phe_mx())
    })
  })
  output$phePlot <- shiny::renderPlot({
    # If no df, gives blank area.
    # Better to either do nothing or do plot_null()
    # use uiOutput and renderUI?
    shiny::req(phe_mx(), cov_df())
    shiny::withProgress(message = 'Pheno Plot ...', value = 0, {
      shiny::setProgress(1)
      plot_sex(phe_mx(), cov_df())
    })
  })
  output$phePlotTable <- shiny::renderUI({
    shiny::req(phe_mx(), cov_df())
    shiny::tagList(
      shiny::plotOutput(ns("phePlot")),
      shiny::tableOutput(ns("phe_sum")))
  })
}

#' UI for shinyPhenoPlot Shiny Module
#'
#' UI for phenotype plots and summaries to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyPhenoPlot
#' @export
shinyPhenoPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("phePlotTable"))
}
