## app.R ##

# Master control file for projects
projects <- read.csv("qtl2shinyData/projects.csv", stringsAsFactors = FALSE)

ui <- qtl2shiny::shinyDash()
server <- function(input, output, session) {
  
  projects_info <- shiny::reactive({projects})
  shiny::callModule(qtl2shiny::shinyMain, "qtl2shiny", projects_info)
  
  # Allow reconnect with Shiny Server.
  session$allowReconnect(TRUE)
}

shinyApp(ui, server)
