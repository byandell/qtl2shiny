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
      menuItem("Haplotype Scans", tabName = "hap_scan",
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
      tabItem(tabName="hap_scan",
              shinyHaploUI("hap_scan")),
      ## Diploid Analysis
      tabItem(tabName="dip_scan", 
              shinyDiploUI("dip_scan"))
    )
  )
)

#####################################################
server <- function(input, output, session) {
  
  ## Data Setup
  pheno_typer <- reactive({pheno_type})
  peaks_tbl <- reactive({peaks})
  pmap_obj <- reactive({pmap})
  analyses_tblr <- reactive({analyses_tbl})
  
  set_par <- callModule(shinySetup, "setup", 
                    pheno_typer, peaks_tbl, 
                    pmap_obj, analyses_tblr, cov_mx)
  
  ## Continue with Plots and Analysis.

  ## Phenotypes and Covariates.
  analyses_df <- reactive({
    phename <- req(set_par()$phe_par$pheno_names)
    ## The analyses_tbl should only have one row per pheno.
    analyses_tbl %>%
      filter(pheno %in% phename)
  })
  phe_df <- reactive({
    ## Make sure we get only one column per distinct pheno.
    get_pheno(pheno_data,
              analyses_df() %>%
                distinct(pheno, .keep_all=TRUE))
  })
  cov_mx <- reactive({
    get_covar(covar, analyses_df())
  })
  
  ## Set up reactives for scan1 module.
  K_chr <- reactive({K[set_par()$win_par$chr_id]})

  ## Haplotype Analysis.
  callModule(shinyHaplo, "hap_scan", 
             set_par()$win_par, pmap_obj, phe_df, cov_mx, K_chr)

  ## Diplotype Analysis.
  callModule(shinyDiplo, "dip_scan",
             set_par()$win_par, phe_df, cov_mx, K_chr)
}

shinyApp(ui, server)
