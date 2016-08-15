
suppressPackageStartupMessages({
  library(doqtl2)
  library(shiny)
  library(dplyr)
})

dirpath <- "~/Documents/Research/attie_alan/DO/data"
datapath <- file.path(dirpath, "DerivedData")

pmap <- readRDS(file.path(datapath, "pmap.rds"))
peaks <- readRDS(file.path(datapath, "peaks.rds"))
peak_info <- peaks$output
analyses_tbl <- readRDS(file.path(datapath, "analyses.rds")) %>%
  filter(output %in% peak_info)
rm(peak_info)

pheno_type <- c("all", sort(unique(analyses_tbl$pheno_type)))
#source("inst/tests/shiny/shinyPhenos.R")

library(shiny)

ui <- fluidPage(
  shinyPhenosUI("phenos"),
  shinyWindowUI("window"),
  uiOutput("sel_phenos"),
  tableOutput("peaks_tbl")
)

server <- function(input, output, session) {
  ## Reactives for testPhenos.
  pheno_typer <- reactive({pheno_type})
  analyses_tblr <- reactive({analyses_tbl})
  peaks_tbl <- reactive({peaks})
  pmap_obj <- reactive({pmap})

  win_par <- callModule(shinyWindow, "window", pmap_obj)

  pheno_re <- callModule(shinyPhenos, "phenos",
             pheno_typer, peaks_tbl, analyses_tblr,
             win_par)

  ## Show selection
  output$sel_phenos <- renderText({
    paste("Selected:", paste(pheno_re(), collapse=", "))
  })

  ## Set up peaks data frame.
  peaks_df <- reactive({
    peaks %>%
      filter(output %in% pheno_re())
  })

  # Output the peaks table
  output$peaks_tbl <- renderTable({peaks_df()})
}

shinyApp(ui = ui, server = server)
