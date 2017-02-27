## app.R ##
suppressPackageStartupMessages({
  library(qtl2shiny)
  library(GGally)
  library(grid)
  library(gridBase)
  library(qtl2geno)
  library(qtl2scan)
})

source("setup4.R")

#####################################################
## Function from Joe Cheng
## https://gist.github.com/jcheng5/5913297
helpPopup <- function(title, content,
                      placement = c('right', 'top', 'left', 'bottom'),
                      trigger = c('click', 'hover', 'focus', 'manual')) {
  shiny::tagList(
    shiny::singleton(
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
ui <- shinydashboard::dashboardPage(skin="red",
  shinydashboard::dashboardHeader(title = "DOQTL2 Wave 4"),
  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      shinySetupOutput("setup"),
      shinydashboard::menuItem("Phenotypes and Region", tabName = "phenos",
               icon = icon("dashboard")),
      shinydashboard::menuItem("Haplotype Scans", tabName = "hap_scan",
               icon = icon("dashboard")),
      shinydashboard::menuItem("SNP/Gene Action", tabName = "dip_scan", 
               icon = icon("dashboard")),
      tags$div(id = "popup",
               helpPopup(NULL,
                         includeMarkdown("about.md"),
                         placement = "right", trigger = "click"))
    )
  ),
  shinydashboard::dashboardBody(
    shinydashboard::tabItems(
      ## Phenotypes and Region
      shinydashboard::tabItem(tabName = "phenos", 
              shinySetupUI("setup")),
      ## Scans
      shinydashboard::tabItem(tabName="hap_scan",
              shinyHaploUI("hap_scan")),
      ## Diploid Analysis
      shinydashboard::tabItem(tabName="dip_scan", 
              shinyDiploUI("dip_scan"))
    )
  )
)

#####################################################
server <- function(input, output, session) {
  
  ## Data Setup
  pheno_typer <- shiny::reactive({pheno_type})
  peaks_tbl <- shiny::reactive({peaks})
  pmap_obj <- shiny::reactive({pmap})
  analyses_tblr <- shiny::reactive({analyses_tbl})
  data_path <- shiny::reactive({datapath})
  
  set_par <- shiny::callModule(shinySetup, "setup", 
                    pheno_typer, peaks_tbl, 
                    pmap_obj, analyses_tblr, cov_mx)
  
  ## Continue with Plots and Analysis.

  ## Phenotypes and Covariates.
  analyses_df <- shiny::reactive({
    phename <- shiny::req(set_par()$phe_par$pheno_names)
    ## The analyses_tbl should only have one row per pheno.
    dplyr::filter(analyses_tbl, pheno %in% phename)
  })
  phe_df <- shiny::reactive({
    ## Make sure we get only one column per distinct pheno.
    DOread::get_pheno(pheno_data,
                      dplyr::distinct(analyses_df(),
                                      pheno, .keep_all=TRUE))
  })
  cov_mx <- shiny::reactive({
    DOread::get_covar(covar, analyses_df())
  })
  
  ## Set up shiny::reactives for scan1 module.
  K_chr <- shiny::reactive({K[set_par()$win_par$chr_id]})

  ## Haplotype Analysis.
  shiny::callModule(shinyHaplo, "hap_scan", 
             set_par()$win_par, pmap_obj, 
             phe_df, cov_mx, K_chr,
             data_path)

  ## Diplotype Analysis.
  shiny::callModule(shinyDiplo, "dip_scan",
             set_par()$win_par, 
             phe_df, cov_mx, K_chr,
             data_path)
}

shinyApp(ui, server)
