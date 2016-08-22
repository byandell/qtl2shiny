## app.R ##
suppressPackageStartupMessages({
  library(qtl2shiny)
  library(GGally)
})

source("setup.R")

#####################################################
## Function from Joe Cheng
## https://gist.github.com/jcheng5/5913297
helpPopup <- function(title, content,
                      placement = c('right', 'top', 'left', 'bottom'),
                      trigger = c('click', 'hover', 'focus', 'manual')) {
  tagList(
    singleton(
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
ui <- dashboardPage(skin="red",
  dashboardHeader(title = "DOQTL2 Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      includeMarkdown("attieDO.md"),
      textOutput("num_pheno"),
      uiOutput("chr_pos"),
      menuItem("Phenotypes and Region", tabName = "phenos",
               icon = icon("dashboard")),
      menuItem("Identify Hotspot", tabName = "peak",
               icon = icon("dashboard")),
      menuItem("Phenotype Info", tabName = "pheno_info",
               icon = icon("dashboard")),
      menuItem("Genome Scans", tabName = "scans",
               icon = icon("dashboard")),
      menuItem("SNP Scans", tabName = "snps",
               icon = icon("dashboard")),
      menuItem("SNP Consequence", tabName = "snp_detail", 
               icon = icon("dashboard")),
      menuItem("Diploid Scan", tabName = "dip_scan", 
               icon = icon("dashboard")),
      tags$div(id = "popup",
               helpPopup(NULL,
                         includeMarkdown("about.md"),
                         placement = "right", trigger = "click"))
    )
  ),
  dashboardBody(
    tabItems(
      ## Genome Region
      tabItem(tabName = "phenos", fluidRow(
        column(6, shinyPhenosUI("phenos")),
        column(4, shinyWindowUI("window")))),
      tabItem(tabName = "peak", fluidRow(
        column(4, shinyPeaksInput("shinypeaks")),
        column(8, shinyPeaksOutput("shinypeaks")))),

      ## Phenotype Information
      tabItem(tabName="pheno_info", 
              tabsetPanel(
                tabPanel("Peaks",
                         dataTableOutput("peaks_tbl")),
                tabPanel("Raw Data",
                         shinyPhenoPlotUI("PhenoPlotRaw")),
                tabPanel("Trans Data",
                         shinyPhenoPlotUI("PhenoPlotTrans")),
                tabPanel("Covariates",
                         dataTableOutput("analyses_tbl")))),

      ## Scans
      tabItem(tabName="scans", fluidRow(
        shinyScan1PlotUI("hap_scan"))),

      tabItem(tabName="snps", fluidRow(
        shinyScan1SNPUI("snp_scan"))),

      ## Diploid Scan
      tabItem(tabName="dip_scan",
              shinyDominanceUI("dip_scan")),

      ## SNP Detail
      tabItem(tabName="snp_detail",
              shinySNPCsqUI("snp_csq"))
    )
  )
)

server <- function(input, output, session) {
  ## Reactives for shinyPhenos.
  pheno_typer <- reactive({pheno_type})
  peaks_tbl <- reactive({peaks})
  pmap_obj <- reactive({pmap})
  analyses_tblr <- reactive({analyses_tbl})

  hot_peak <- callModule(shinyPeaks, "shinypeaks",
                         pheno_typer, peaks_tbl, pmap_obj)

  ## Use peaks as input to shinyWindow.
  win_par <- callModule(shinyWindow, "window",
                        pmap_obj, hot_peak)

  chr_pos <- reactive({
    make_chr_pos(win_par$chr_id, 
                 win_par$peak_Mbp, win_par$window_Mbp)
  })
  output$chr_pos <- renderText({
    paste0("Region: ", chr_pos(), "Mbp")
  })
  output$num_pheno <- renderText({
    num_pheno(character(), analyses_tbl)
  })
  observeEvent(pheno_anal(), {
    output$num_pheno <- renderText({
      pheno <- req(pheno_anal())
      num_pheno(pheno, analyses_tbl)
    })
  })

  ## Use window as input to shinyPhenos.
  pheno_anal <- callModule(shinyPhenos, "phenos",
                         pheno_typer, peaks_tbl, analyses_tblr,
                         hot_peak, win_par)

  ##############################################
  ## Continue with Plots and Analysis.

  ## Set up peaks data frame.
  peaks_df <- reactive({
    phenos <- req(pheno_anal())
    peaks %>%
      filter(output %in% names(phenos))
  })

  ## Set up analyses data frame.
  analyses_df <- reactive({
    phenos <- req(pheno_anal())
    analyses_tbl %>%
      filter(output %in% names(phenos))
  })

  # Output the analyses table
  output$analyses_tbl <- renderDataTable({
    collapse_covar(analyses_df())
  }, options = list(scrollX = TRUE, pageLength = 10))

  # Output the peaks table
  output$peaks_tbl <- renderDataTable({
    peaks_df() %>%
      select(pheno,chr,pos,lod)
  }, options = list(scrollX = TRUE, pageLength = 10))

  ## Reactive for phenotypes.
  phe_df <- reactive({
    get_pheno(pheno_data,
              analyses_df() %>%
                distinct(pheno, .keep_all=TRUE))
  })
  raw_phe_df <- reactive({
    get_pheno(pheno_data,
              analyses_df() %>%
                distinct(pheno, .keep_all=TRUE),
              FALSE)
  })

  ## Reactive for covariates
  cov_mx <- reactive({
    get_covar(covar, analyses_df())
  })

  ## Density or scatter plot of phenotypes.
  callModule(shinyPhenoPlot, "PhenoPlotRaw", raw_phe_df, cov_mx)
  callModule(shinyPhenoPlot, "PhenoPlotTrans", phe_df, cov_mx)
  
  ## Set up reactives for scan1 module.
  chr_id <- reactive({as.character(req(win_par$chr_id))})
  K_chr <- reactive({K[chr_id()]})

  ## Genome scans by haplotype alleles and SNPs.
  probs_obj <- reactive({
    req(chr_id())
    withProgress(message = 'Read probs ...', value = 0, {
      setProgress(1)
      read_probs(chr_id(), datapath)
    })
  })
  callModule(shinyScan1Plot, "hap_scan", 
                             chr_id, phe_df, cov_mx,
                             pheno_anal, probs_obj, K_chr)
  
  snp_scan_obj <- callModule(shinyScan1SNP, "snp_scan",
                             win_par, phe_df, cov_mx, 
                             pheno_anal, probs_obj, K_chr)
  
  ## Dominance
  callModule(shinyDominance, "dip_scan",
                             win_par, phe_df, cov_mx,
                             pheno_anal, probs_obj, K_chr)
  
  ## SNP Consequence
  callModule(shinySNPCsq, "snp_csq", 
             snp_scan_obj, chr_pos)
}

shinyApp(ui, server)
