#' Shiny app for qtl2
#'
#' Run shiny app for qtl2 data.
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#' 
#' @return No return value; called for side effects.
#'
#' @details 
#' See qtl2shinyData vignette for data setup.
#' 
#' @examples
#' \dontrun{qtl2shinyApp()}
#'
#' @export
#' @importFrom shiny runApp
qtl2shinyApp <- function() {
  shiny::runApp(system.file(file.path("qtl2shinyApp", "app.R"), package='qtl2shiny'))
}
