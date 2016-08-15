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
  selectInput("plotType", "Plot Type",
              c("density", "scans", "snps")),
  conditionalPanel(condition="input.plotType == 'density'",
                   plotOutput("distPlot")),
  conditionalPanel(condition="input.plotType == 'scans'",
                   shinyScan1UI("genome_scan")),
  conditionalPanel(condition="input.plotType == 'snps'",
                   shinyScan1UI("snp_scan"))
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
  chr_id <- reactive({as.character(req(input$chr_id))})

  ## Reactive for phenotypes.
  phe_df <- reactive({
    get_pheno(pheno_data,
              analyses_df() %>%
                distinct(pheno, .keep_all=TRUE))
  })

  ## Reactive for covariates
  cov_mx <- reactive({
    get_covar(covar, analyses_df())
  })

  plot_type <- reactive({input$plotType})

  # Choose Phenotypes for Analysis.
  output$choose_phenoanal <- renderUI({
    phenames <- pheno_re()
    choices <- c("all","none", phenames)
    checkboxGroupInput("pheno_anal", "Choose phenotypes",
                       choices = choices, inline=TRUE)
  })
  pheno_anal_list <- reactive({req(input$pheno_anal)})

  ## Genome scan.
  K_chr <- reactive({K[chr_id()]})
  pmap_obj <- reactive({pmap})
  probs_obj <- reactive({
    read_probs(chr_id(), datapath)
  })
  callModule(shinyScan1, "genome_scan", plot_type,
             chr_id, phe_df, cov_mx, pheno_anal_list, analyses_df,
             probs_obj, K_chr)
}

shinyApp(ui = ui, server = server)
