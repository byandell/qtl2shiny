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
shinyPhenoPlot <- function(input, output, session,
                           phe_df, cov_mx) {
  ## Scatter plot or density
  output$phe_sum <- renderTable({
    withProgress(message = 'Pheno Summary ...', value = 0, {
      setProgress(1)
      summary_na(phe_df())
    })
  })
  output$phePlot <- renderPlot({
    withProgress(message = 'Pheno Plot ...', value = 0, {
      setProgress(1)
      plot_sex(phe_df(), cov_mx())
    })
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
  ns <- NS(id)
  tagList(
    plotOutput(ns("phePlot")),
    tableOutput(ns("phe_sum"))
  )
}
