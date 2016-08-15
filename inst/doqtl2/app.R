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
      menuItem("Phenotype Info", tabName = "tables", icon = icon("dashboard"),
               collapsible =
                 menuSubItem("Phenotype Peaks", tabName = "peaks_tbl"),
                 menuSubItem("Analyses Details", tabName = "analyses_tbl"),
                 menuSubItem("Phenotype Plots", tabName = "pheno_plot")
      ),
      menuItem("Genome Scans", tabName = "scans",
               icon = icon("dashboard")),
      menuItem("SNP Scans", tabName = "snps",
               icon = icon("dashboard")),
      menuItem("SNP Detail", tabName = "snp_detail", icon = icon("dashboard"),
               collapsible =
               menuSubItem("Top SNPs", tabName = "top_snps_tbl"),
               menuSubItem("Genes in SNP Region", tabName = "gene_region"),
               menuSubItem("Top Features in SNP Region", tabName = "top_feature"),
               menuSubItem("Genes and Exons in SNP Region", tabName = "gene_exon")
      ),
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
        shinyPhenosUI("phenos"),
        shinyWindowUI("window"))),
      tabItem(tabName = "peak", fluidRow(
        shinyPeaksInput("shinypeaks"),
        shinyPeaksOutput("shinypeaks"))),

      ## Tables and Plots
      tabItem(tabName="analyses_tbl", fluidRow(
        dataTableOutput("analyses_tbl"))),
      tabItem(tabName="peaks_tbl", fluidRow(
        dataTableOutput("peaks_tbl"))),
      tabItem(tabName="pheno_plot", fluidRow(
        shinyPhenoPlotUI("PhenoPlot"))),
      tabItem(tabName="scans", fluidRow(
        shinyScan1UI("genome_scan"))),
      tabItem(tabName="snps", fluidRow(
        shinyScan1UI("snp_scan"))),
      tabItem(tabName="top_snps_tbl", fluidRow(
        shinySNPUI("best_snp"))),
      tabItem(tabName="gene_region", fluidRow(
        shinyGeneRegionUI("gene_region"))),
      tabItem(tabName="top_feature", fluidRow(
        shinyTopFeatureUI("top_feature"))),
      tabItem(tabName="gene_exon", fluidRow(
        shinyGeneExonUI("gene_exon")))
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
    peaks_df()
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
  callModule(shinyPhenoPlot, "PhenoPlot", raw_phe_df, phe_df, cov_mx)

  ## Set up reactives for scan1 module.
  chr_id <- reactive({as.character(req(win_par$chr_id))})
  K_chr <- reactive({K[chr_id()]})

  ## Genome scan.
  probs_obj <- reactive({
    req(chr_id())
    withProgress(message = 'Read probs ...', value = 0, {
      setProgress(1)
      read_probs(chr_id(), datapath)
    })
  })
  plot_scans <- reactive({"scans"})
  callModule(shinyScan1, "genome_scan", plot_scans,
                         chr_id, phe_df, cov_mx, pheno_anal,
                         probs_obj, K_chr)

  ## SNP analyses.
  snpprobs_obj <- reactive({
    window_Mbp <- req(win_par$window_Mbp)
    peak_Mbp <- req(win_par$peak_Mbp)
    withProgress(message = 'SNP Scan ...', value = 0, {
      setProgress(1)
      get_snpprobs(chr_id(), peak_Mbp, window_Mbp,
                   names(phe_df()), probs_obj(),
                   pattern = "AB1NZCPW", datapath)
    })
  })
  plot_snps <- reactive({"snps"})
  snp_scan_obj <- callModule(shinyScan1, "snp_scan", plot_snps,
                             chr_id, phe_df, cov_mx, pheno_anal,
                             snpprobs_obj, K_chr)
  top_snps_tbl <- reactive({
    withProgress(message = 'Get Top SNPs ...', value = 0, {
      setProgress(1)
      get_top_snps_tbl(snp_scan_obj())
    })
  })
  callModule(shinySNP, "best_snp", chr_pos, top_snps_tbl)
  feature_file <- reactive({file.path(datapath, "mgi_db.sqlite")})
  callModule(shinyGeneRegion, "gene_region",
             top_snps_tbl, feature_file)
  gene_exon_tbl <- reactive({
    withProgress(message = 'Gene Exon Calc ...', value = 0, {
      setProgress(1)
      sql_filename <- req(feature_file())
      top_snps_set <- req(top_snps_tbl())
      get_gene_exon_snp(top_snps_set, sql_filename)
    })
  })
  callModule(shinyTopFeature, "top_feature",
             chr_pos, snp_scan_obj, top_snps_tbl, gene_exon_tbl)
  callModule(shinyGeneExon, "gene_exon",
             chr_pos, top_snps_tbl, gene_exon_tbl)
}

shinyApp(ui, server)
