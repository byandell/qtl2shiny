#' Shiny app for doqtl2
#'
#' run shiny app for doqtl2
#'
#' @param appfile name of app file for shiny
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{doqtl2_app()}
#'
#' @export
doqtl2_app <- function(appfile = "app.R") {
  shiny::runApp(system.file(file.path("shiny", appfile), package='doqtl2'))
}
