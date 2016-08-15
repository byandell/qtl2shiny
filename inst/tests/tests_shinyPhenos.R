
suppressPackageStartupMessages({
#  library(doqtl2)
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

pheno_type <- c(sort(unique(analyses_tbl$pheno_type)), "all")

source("R/shinyPhenos.R")

library(shiny)

ui <- fluidPage(
  shinyPhenosUI("phenos"),
  uiOutput("choose_phenoanal"),
  uiOutput("choose_chr"),
  uiOutput("choose_peak"),
  sliderInput("window_Mbp", "Window over Peak if Positive",
              0, 6, 0, step=0.5)
)

server <- function(input, output, session) {
  ## Reactives for testPhenos.
  pheno_typer <- reactive({pheno_type})
  analyses_tblr <- reactive({analyses_tbl})
  peaks_tbl <- reactive({peaks})

  # Select chromosome
  output$choose_chr <- renderUI({
    selectInput("chr_id", strong("Key Chromosome"),
                choices = c(as.character(1:19),"X"))
  })
  chr_id <- reactive({input$chr_id})

  # Position slider
  output$choose_peak <- renderUI({
    req(input$chr_id)
    rng <- round(range(pmap[[input$chr_id]]), 2)
    sliderInput("peak_Mbp", "Peak Position (Mbp)", rng[1], rng[2], mean(rng),
                step=.5)
  })
  peak_Mbp <- reactive({input$peak_Mbp})
  window_Mbp <- reactive({input$window_Mbp})

  pheno_re <- callModule(shinyPhenos, "phenos",
                         pheno_typer, peaks_tbl, analyses_tblr,
                         chr_id, peak_Mbp, window_Mbp)

  ## Stuff below works, but want more of this in the module.
  ## Could likely put peak stuff in easily.
  ## Probably length limit (so don't pass long pheno_re() back and forth).

  # Choose Phenotypes for Analysis.
  output$choose_phenoanal <- renderUI({
    phenames <- pheno_re()
    choices <- c("all","none", phenames)
    checkboxGroupInput("pheno_anal", "Choose phenotypes",
                       choices = choices, inline=TRUE)
  })

  ## Update list of pheno_anal based on analyses_df()$output.
  observe({
    ## pheno_anal from analyses_df().
    cat(file=stderr(),"choose_phenoanal\n")
    phenames <- pheno_re()
    cat(file=stderr(), "phe:", paste(phenames, collapse=","), "\n")
    ## Could put something in here to shorten length of phenames with more ...
    ## but low priority for now.

    ## Selected from most recent input.
    selected <- input$pheno_anal
    cat(file=stderr(), "sel:", paste(selected, collapse=","), "\n")

    if("all" %in% selected)
      selected <- c(selected[!(selected %in% c("all","none"))], phenames)
    if("none" %in% selected)
      selected <- ""

    ## Update phenames to include selected (but not "")
    phenames <- unique(c(selected, sort(phenames)))
    phenames <- phenames[phenames != ""]

    ## Change label if filtering by peak.
    newlabel <- "Choose Phenotypes"
    if(length(input$use_pos)) {
      if(input$use_pos & req(window_Mbp()) > 0) {
        newlabel <- "Phenotypes near Peak"
      }
    }

    ## Now update list of phenotypes.
    updateCheckboxGroupInput(session, "pheno_anal",
                             label=newlabel,
                             choices=c("all","none", phenames),
                             selected=selected, inline=TRUE)
  })
  output$sel_phenos <-
    renderText({paste("Phenos:", paste(pheno_re(), collapse=","))})
}

shinyApp(ui = ui, server = server)
