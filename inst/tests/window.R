
suppressPackageStartupMessages({
  library(doqtl2)
  library(shiny)
  library(dplyr)
})

dirpath <- "~/Documents/Research/attie_alan/DO/data"
datapath <- file.path(dirpath, "DerivedData")

pmap <- readRDS(file.path(datapath, "pmap.rds"))
#source("inst/tests/shiny/shinyWindow.R")

library(shiny)

ui <- fluidPage(
  shinyWindowUI("window"),
  uiOutput("chr_id"),
  uiOutput("peak_Mbp"),
  uiOutput("window_Mbp")
)

server <- function(input, output, session) {
  ## Reactives for testPhenos.
  pmap_obj <- reactive({pmap})

  win_par <- callModule(shinyWindow, "window", pmap_obj)

  ## Show selections
  output$chr_id <- renderText({
    paste("Chr:", win_par$chr_id)
  })
  output$peak_Mbp <- renderText({
    paste("Peak:", win_par$peak_Mbp)
  })
  output$window_Mbp <- renderText({
    paste("Window:", win_par$window_Mbp)
  })
}

shinyApp(ui = ui, server = server)
