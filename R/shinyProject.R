#' Shiny project selection
#'
#' Shiny module for selection of project, with interface \code{shinyProjectUI}.
#'
#' @param input,output,session,projects_info standard shiny arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom shiny callModule NS isTruthy 
#'   renderUI uiOutput selectInput
#' @importFrom dplyr distinct filter
#' @importFrom rlang .data
#' 
shinyProject <- function(input, output, session, projects_info) {
  ns <- session$ns
  
  output$project <- shiny::renderUI({
    shiny::req(projects_info())
    choices <- unique(projects_info()$project)
    if(is.null(selected <- input$project)) {
      selected <- choices[1]
    }
    shiny::selectInput(ns("project"), "Project", choices, selected)
  })
  
  shiny::reactive({
    shiny::req(projects_info())
    project_id <- NULL
    if(shiny::isTruthy(input$project)) {
      project_id <- input$project
    }
    if(is.null(project_id)) {
      project_id <- projects_info()$project[1]
    }
    dplyr::distinct(
      dplyr::filter(
        projects_info(),
        .data$project == project_id),
      .data$project, .keep_all = TRUE)
  })
}
shinyProjectUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("project"))
}
