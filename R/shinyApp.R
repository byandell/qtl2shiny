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

#' @export
#' @rdname qtl2shinyApp
#' 
shinyDash <- function() {
  shinydashboard::dashboardPage(
    skin="red",
    shinydashboard::dashboardHeader(
      title = "qtl2shiny"),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        qtl2shiny::shinyMainInput("qtl2shiny"),
        shinydashboard::menuItem(
          "Phenotypes and Region",
          tabName = "phenos",
          icon = icon("dashboard", verify_fa = FALSE)),
        shinydashboard::menuItem(
          "Haplotype Scans",
          tabName = "hap_scan",
          icon = icon("dashboard", verify_fa = FALSE)),
        shinydashboard::menuItem(
          "SNP/Gene Action",
          tabName = "dip_scan",
          icon = icon("dashboard", verify_fa = FALSE)),
        tags$div(
          id = "popup",
          helpPopup(
            "qtl2shiny help",
            shiny::includeMarkdown(system.file(file.path("qtl2shinyApp", "about.md"), package='qtl2shiny')),
            placement = "right", trigger = "click"))
      )
    ),
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        ## Phenotypes and Region
        shinydashboard::tabItem(
          tabName = "phenos",
          qtl2shiny::shinyMainUI("qtl2shiny")),
        ## Scans
        shinydashboard::tabItem(
          tabName="hap_scan",
          qtl2shiny::shinyMainOutput("qtl2shiny")),
        ## Diploid Analysis
        shinydashboard::tabItem(
          tabName="dip_scan", 
          qtl2shiny::shinyMainOutput2("qtl2shiny"))
      )
    )
  )
}
