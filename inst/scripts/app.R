## app.R ##

# Master control file for projects
projects <- read.csv("projects.csv", stringsAsFactors = FALSE)
# Assumes all project data in RDS format
read_project_data <- qtl2shiny::read_project_rds

#####################################################
## Function from Joe Cheng
## https://gist.github.com/jcheng5/5913297
helpPopup <- function(title, content,
                      placement = c('right', 'top', 'left', 'bottom'),
                      trigger = c('click', 'hover', 'focus', 'manual')) {
  shiny::tagList(
    shiny::singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover()})"),
        tags$style(type = "text/css", ".popover{max-width:500px; position: fixed;}")
      )
    ),
    tags$a(
      href = "#", class = "btn btn-link",
      `data-toggle` = "popover", `data-html` = "true",
      title = title, `data-content` = content, `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok = TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok = TRUE)[1],
      "Help ..."
    )
  )
}

#####################################################
ui <- shinydashboard::dashboardPage(
  skin="red",
  shinydashboard::dashboardHeader(
    title = "qtl2shiny"),
  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      qtl2shiny::shinyMainInput("qtl2shiny"),
      shinydashboard::menuItem(
        "Phenotypes and Region",
        tabName = "phenos",
        icon = icon("dashboard")),
      shinydashboard::menuItem(
        "Haplotype Scans",
        tabName = "hap_scan",
        icon = icon("dashboard")),
      shinydashboard::menuItem(
        "SNP/Gene Action",
        tabName = "dip_scan",
        icon = icon("dashboard")),
      tags$div(
        id = "popup",
        helpPopup(
          NULL,
          shiny::includeMarkdown("about.md"),
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

#####################################################
server <- function(input, output, session) {
  
  projects_info <- shiny::reactive({projects})
  shiny::callModule(qtl2shiny::shinyMain, "qtl2shiny", projects_info)

  # Allow reconnect with Shiny Server.
  session$allowReconnect(TRUE)
}

shinyApp(ui, server)
