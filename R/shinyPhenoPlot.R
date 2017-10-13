#' Shiny Phenotype Plot module
#'
#' Shiny module to plot phenotypes.
#'
#' @param input,output,session standard shiny arguments
#' @param phe_df,cov_mx reactive arguments
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
shinyPhenoPlot <- function(input, output, session,
                           phe_df, cov_mx) {
  ## Scatter plot or density
  output$phe_sum <- shiny::renderTable({
    shiny::req(phe_df())
    shiny::withProgress(message = 'Pheno Summary ...', value = 0, {
      shiny::setProgress(1)
      summary_na(phe_df())
    })
  })
  output$phePlot <- shiny::renderPlot({
    # If no df, gives blank area.
    # Better to either do nothing or do plot_null()
    # use uiOutput and renderUI?
    shiny::req(phe_df(), cov_mx())
    shiny::withProgress(message = 'Pheno Plot ...', value = 0, {
      shiny::setProgress(1)
      plot_sex(phe_df(), cov_mx())
    })
  })
  output$phePlotTable <- shiny::renderUI({
    shiny::req(phe_df(), cov_mx())
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
