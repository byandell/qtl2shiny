#' Shiny setup for DOQTL selection
#'
#' Shiny module for phenotype selection.
#'
#' @param input,output,session standard shiny arguments
#' @param pheno_typer,peaks_tbl,pmap_obj,analyses_tbl,cov_mx reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom dplyr distinct filter select
#' @importFrom DOread get_pheno
#' @importFrom shiny callModule NS reactive req 
#'   radioButtons selectInput
#'   dataTableOutput textOutput uiOutput
#'   renderDataTable renderText renderUI
#'   observeEvent
#'   strong tagList
shinySetup <- function(input, output, session,
                        pheno_typer, peaks_tbl, pmap_obj, analyses_tbl, cov_mx) {
  ns <- session$ns
  
  # Select phenotype dataset
  pheno_group <- shiny::reactive({
    sort(unique(shiny::req(analyses_tbl())$pheno_group))
  })
  pheno_type <- shiny::reactive({
    phe_gp <- shiny::req(input$pheno_group)
    analyses_group <- 
      dplyr::filter(
        shiny::req(analyses_tbl()),
        pheno_group %in% phe_gp)
    sort(unique(analyses_group$pheno_type))
  })
  output$pheno_group <- shiny::renderUI({
    shiny::req(choices <- pheno_group())
    if(is.null(selected <- input$pheno_group)) {
      selected <- choices[1]
    }
    shiny::selectInput(ns("pheno_group"), "Phenotype Group",
                       choices = as.list(choices),
                       selected = selected,
                       multiple = TRUE)
  })
  output$dataset <- shiny::renderUI({
    choices <- c("all", shiny::req(pheno_type()))
    if(is.null(selected <- input$dataset))
      selected <- NULL
    shiny::selectInput(ns("dataset"), "Phenotype Set",
                choices = as.list(choices),
                selected = selected,
                multiple = TRUE)
  })
  
  ## Locate Peak.
  win_par <- shiny::callModule(shinyPeaks, "shinypeaks",
                         input, pheno_type, peaks_tbl, pmap_obj)
  
  chr_pos <- shiny::reactive({
    make_chr_pos(win_par$chr_id, 
                 win_par$peak_Mbp, win_par$window_Mbp)
  })
  output$chr_pos <- shiny::renderText({
    paste0("Region: ", chr_pos(), "Mbp")
  })
  output$num_pheno <- shiny::renderText({
    num_pheno(character(), analyses_tbl())
  })
  shiny::observeEvent(phe_par$pheno_names, {
    output$num_pheno <- shiny::renderText({
      num_pheno(phe_par$pheno_names, analyses_tbl())
    })
  })
  
  ## Use window as input to shinyPhenos.
  phe_par <- shiny::callModule(shinyPhenos, "phenos",
             input, peaks_tbl, analyses_tbl, win_par)
  
  ## Density or scatter plot of phenotypes.
  analyses_df <- shiny::reactive({
    dplyr::filter(analyses_tbl(), pheno %in% shiny::req(phe_par$pheno_names))
  })
  phe_df <- shiny::reactive({
    read_pheno(pheno_data, analyses_df())
  })
  raw_phe_df <- shiny::reactive({
    read_pheno(pheno_data, analyses_df(), FALSE)
  })
  shiny::callModule(shinyPhenoPlot, "PhenoPlotRaw", raw_phe_df, cov_mx)
  shiny::callModule(shinyPhenoPlot, "PhenoPlotTrans", phe_df, cov_mx)
  
  # Output the analyses table
  output$analyses_tbl <- shiny::renderDataTable({
    collapse_covar(analyses_df())
  }, options = list(scrollX = TRUE, pageLength = 10))
  
  # Output the peaks table
  
  ## Set up peaks data frame.
  peaks_df <- shiny::reactive({
    phename <- shiny::req(phe_par$pheno_names) 
    dplyr::filter(peaks, pheno %in% phename)
  })
  output$peaks_tbl <- shiny::renderDataTable({
    dplyr::select(peaks_df(), pheno, chr, pos, lod)
  }, options = list(scrollX = TRUE, pageLength = 10))
  
  ## Setup input logic.
  output$title <- shiny::renderUI({
    shiny::strong(shiny::req(input$setup))
  })
  output$sidebar_setup <- shiny::renderUI({
    switch(shiny::req(input$setup),
           Phenotypes = shiny::tagList(
             shinyPhenosUI(ns("phenos")),
             shiny::radioButtons(ns("show_data"), NULL,
                          c("LOD Peaks","Covariates",
                            "Trans Data","Raw Data"),
                          input$show_data)),
           Region = shinyPeaksInput(ns("shinypeaks")))
  })
  output$main_setup <- shiny::renderUI({
    switch(shiny::req(input$setup),
           Phenotypes = {
             switch(shiny::req(input$show_data),
                    "LOD Peaks"  = shiny::dataTableOutput(ns("peaks_tbl")),
                    "Raw Data"   = shinyPhenoPlotUI(ns("PhenoPlotRaw")),
                    "Trans Data" = shinyPhenoPlotUI(ns("PhenoPlotTrans")),
                    "Covariates" = shiny::dataTableOutput(ns("analyses_tbl")))
             },
           Region = shinyPeaksOutput(ns("shinypeaks")))
  })
  
  output$radio <- shiny::renderUI({
    shiny::radioButtons(ns("setup"), NULL,
                 c("Region", "Phenotypes"),
                 input$setup,
                 inline=TRUE)
  })
  
  ## Return.
  shiny::reactive({
    list(phe_par = phe_par,
         win_par = win_par)
  })
}
#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinySetup
#' @export
shinySetupOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::textOutput(ns("num_pheno")),
    shiny::uiOutput(ns("chr_pos"))
  )
}
#' @rdname shinySetup
#' @export
shinySetupUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    sidebarPanel(
      shiny::uiOutput(ns("title")),
      shiny::uiOutput(ns("radio")),
      shiny::uiOutput(ns("sidebar_setup")),
      shiny::uiOutput(ns("pheno_group")),
      shiny::uiOutput(ns("dataset"))
    ),
    mainPanel(shiny::uiOutput(ns("main_setup")))
  )
}
