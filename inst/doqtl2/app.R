## app.R ##

suppressPackageStartupMessages({
  library(qtl2shiny)
  library(qtl2feather)
})

source("setup.R")

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
  shinydashboard::dashboardHeader(title = "qtl2shiny"),
  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      shinySetupOutput("setup"),
      shinydashboard::menuItem("Project", tabName = "project",
                               icon = icon("dashboard")),
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
              shinyDiploUI("dip_scan")),
      ## Project
      shinydashboard::tabItem(tabName="project", 
                              shinyProjectUI("project"))
    )
  )
)

#####################################################
server <- function(input, output, session) {
  
  ## Data Setup
  projects_info <- shiny::reactive({projects})
  project_info <- shiny::callModule(shinyProject, "project", projects_info)
  pheno_type <- shiny::reactive({
    shiny::req(project_info())
    read_project_data(project_info(), "pheno_type")
    })
  peaks <- shiny::reactive({
    shiny::req(project_info())
    read_project_data(project_info(), "peaks")
  })
  analyses_tbl <- shiny::reactive({
    shiny::req(project_info())
    ## The analyses_tbl should only have one row per pheno.
    read_project_data(project_info(), "analyses")
  })
  pheno_data <- shiny::reactive({
    shiny::req(project_info())
    read_project_data(project_info(), "pheno_data")
  })
  covar <- shiny::reactive({
    shiny::req(project_info())
    read_project_data(project_info(), "covar")
  })
  
  pmap_obj <- shiny::reactive({
    shiny::req(project_info())
    read_project_data(project_info(), "pmap")
  })
  kinship <- shiny::reactive({
    shiny::req(project_info())
    read_project_data(project_info(), "kinship")
  })

  set_par <- shiny::callModule(shinySetup, "setup", 
                    pheno_type, peaks, 
                    pmap_obj, analyses_tbl, 
                    cov_mx, pheno_data, project_info)
  
  ## Continue with Plots and Analysis.

  ## Phenotypes and Covariates.
  analyses_df <- shiny::reactive({
    phename <- shiny::req(set_par()$phe_par$pheno_names)
    dplyr::filter(analyses_tbl(), pheno %in% phename)
  })
  phe_df <- shiny::reactive({
    shiny::req(analyses_df())
    pheno_read(pheno_data(), analyses_df())
  })
  cov_mx <- shiny::reactive({
    DOread::get_covar(covar(), analyses_df())
  })
  
  ## Set up shiny::reactives for scan1 module.
  K_chr <- shiny::reactive({
    kinship()[set_par()$win_par$chr_id]
  })

  ## Haplotype Analysis.
  shiny::callModule(shinyHaplo, "hap_scan", 
             set_par()$win_par, pmap_obj, 
             phe_df, cov_mx, K_chr, analyses_df, 
             covar, pheno_data, analyses_tbl, peaks,
             project_info)

  ## Diplotype Analysis.
  shiny::callModule(shinyDiplo, "dip_scan",
             set_par()$win_par, 
             phe_df, cov_mx, K_chr, analyses_df,
             project_info)
  
  # Allow reconnect with Shiny Server.
  session$allowReconnect(TRUE)
}

shinyApp(ui, server)
