## app.R ##

# Used here in conjunction with RStudio Connect and UW ResearchDrive
# See https://docs.rstudio.com/connect/user/git-backed/

# Needed libraries
devtools::install_github("byandell/intermediate")
devtools::install_github("byandell/qtl2mediate")
devtools::install_github("byandell/qtl2shiny")

# Master control file for projects
projects <- read.csv("qtl2shinyData/projects.csv", stringsAsFactors = FALSE)

ui <- qtl2shiny::shinyDash()
server <- function(input, output, session) {
  
  projects_info <- shiny::reactive({projects})
  shiny::callModule(qtl2shiny::shinyMain, "qtl2shiny", projects_info)
  
  # Allow reconnect with Shiny Server.
  session$allowReconnect(TRUE)
}

browser()

shinyApp(ui, server)
