## Fairly clunky with "all" option to chr and selecting ind chr
## Want multiple colors by pheno_type
## Want to select one or more pheno_type
## Want ggplot version of plot_scan1
## Want to return best peak per chr for those selected, and window.

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

library(shiny)
library(qtl2plot)

ui <- fluidPage(
  shinyWindowUI("window"),
  shinyPeaksInput("shinypeaks"),
  shinyPeaksOutput("shinypeaks"),
  tableOutput("peak_tbl"),
  uiOutput("chr_id"),
  uiOutput("peak_Mbp"),
  uiOutput("window_Mbp")
)

server <- function(input, output, session) {
  ## Reactives for testPhenos.
  pheno_typer <- reactive({pheno_type})
  peaks_tbl <- reactive({peaks})
  pmap_obj <- reactive({pmap})

  peaks_re <- callModule(shinyPeaks, "shinypeaks",
             pheno_typer, peaks_tbl, pmap_obj)
  output$peak_tbl <- renderTable({peaks_re()})

  ## Use peaks as input to shinyWindow.
  win_par <- callModule(shinyWindow, "window", pmap_obj,
                        peaks_re)

}

shinyApp(ui = ui, server = server)
