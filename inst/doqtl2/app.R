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
      shinySetupOutput("setup"),
      menuItem("Phenotypes and Region", tabName = "phenos",
               icon = icon("dashboard")),
      menuItem("Genome Scans", tabName = "scans",
               icon = icon("dashboard")),
      menuItem("SNP Scans", tabName = "snps",
               icon = icon("dashboard")),
      menuItem("SNP Consequence", tabName = "snp_detail", 
               icon = icon("dashboard")),
      menuItem("SNP/Gene Action", tabName = "dip_scan", 
               icon = icon("dashboard")),
      tags$div(id = "popup",
               helpPopup(NULL,
                         includeMarkdown("about.md"),
                         placement = "right", trigger = "click"))
    )
  ),
  dashboardBody(
    tabItems(
      ## Phenotypes and Region
      tabItem(tabName = "phenos", 
              shinySetupUI("setup")),

      ## Scans
      tabItem(tabName="scans", fluidRow(
        shinyScan1PlotUI("hap_scan"))),

      tabItem(tabName="snps", tagList(
        sidebarPanel(shinyScan1SNPUI("snp_scan")),
        mainPanel(shinyScan1SNPOutput("snp_scan")))),

      ## Diploid Scan
      tabItem(tabName="dip_scan", tagList(
        sidebarPanel(shinyDominanceUI("dip_scan")),
        mainPanel(shinyDominanceOutput("dip_scan")))),

      ## SNP Detail
      tabItem(tabName="snp_detail", tagList(
        sidebarPanel(shinySNPCsqUI("snp_csq")),
        mainPanel(shinySNPCsqOutput("snp_csq"))))
    )
  )
)

server <- function(input, output, session) {
  
  ## Data Setup
  pheno_typer <- reactive({pheno_type})
  peaks_tbl <- reactive({peaks})
  pmap_obj <- reactive({pmap})
  analyses_tblr <- reactive({analyses_tbl})
  
  out <- callModule(shinySetup, "setup", 
                    pheno_typer, peaks_tbl, 
                    pmap_obj, analyses_tblr, cov_mx)
  
  pheno_anal <- reactive({out()$pheno_anal})
  win_par <- reactive({out()$win_par})

  ## Continue with Plots and Analysis.

  ## Phenotypes and Covariates.
  analyses_df <- reactive({
    phenos <- req(pheno_anal())
    analyses_tbl %>%
      filter(output %in% names(phenos))
  })
  phe_df <- reactive({
    get_pheno(pheno_data,
              analyses_df() %>%
                distinct(pheno, .keep_all=TRUE))
  })
  cov_mx <- reactive({
    get_covar(covar, analyses_df())
  })
  
  ## Set up reactives for scan1 module.
  chr_id <- reactive({as.character(req(win_par()$chr_id))})
  K_chr <- reactive({K[chr_id()]})

  ## Genotype Probabilities.
  probs_obj <- reactive({
    req(chr_id())
    withProgress(message = 'Read probs ...', value = 0, {
      setProgress(1)
      read_probs(chr_id(), datapath)
    })
  })
  
  ## Genome Scan.
  callModule(shinyScan1Plot, "hap_scan", 
                             chr_id, phe_df, cov_mx,
                             pheno_anal, probs_obj, K_chr)
  
  ## SNP Scan.
  snp_scan_obj <- callModule(shinyScan1SNP, "snp_scan",
                             win_par, phe_df, cov_mx, 
                             pheno_anal, probs_obj, K_chr)
  
  ## SNP Consequence.
  callModule(shinySNPCsq, "snp_csq", 
             win_par, snp_scan_obj)
  
  ## SNP/Gene Action.
  callModule(shinyDominance, "dip_scan",
             win_par, phe_df, cov_mx,
             pheno_anal, probs_obj, K_chr)
}

shinyApp(ui, server)
