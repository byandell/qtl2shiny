#' Shiny setup for DOQTL selection
#'
#' Shiny module for phenotype selection.
#'
#' @param input,output,session standard shiny arguments
#' @param pheno_typer,peaks_tbl,pmap_obj,analyses_tbl,cov_df,pheno_data,projects_info reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom gdata humanReadable
#' @importFrom dplyr arrange desc distinct filter mutate one_of select 
#' @importFrom tidyr unite
#' @importFrom qtl2pattern pheno_trans
#' @importFrom shiny callModule NS reactive req 
#'   radioButtons selectInput
#'   dataTableOutput textOutput uiOutput
#'   renderDataTable renderText renderUI
#'   observeEvent
#'   strong tagList
shinySetup <- function(input, output, session,
                       pheno_typer, peaks_tbl, pmap_obj, analyses_tbl, 
                       cov_df, pheno_data, projects_info) {
  ns <- session$ns

  project_info <- shiny::callModule(shinyProject, "project", projects_info)
  
  # Select phenotype dataset
  pheno_group <- shiny::reactive({
    cat("pheno_group\n", file = stderr())
    shiny::req(project_info())
    sort(unique(shiny::req(analyses_tbl())$pheno_group))
  })
  pheno_type <- shiny::reactive({
    cat("pheno_type\n", file = stderr())
    shiny::req(project_info())
    phe_gp <- shiny::req(input$pheno_group)
    analyses_group <- 
      dplyr::filter(
        shiny::req(analyses_tbl()),
        pheno_group %in% phe_gp)
    sort(unique(analyses_group$pheno_type))
  })
  output$pheno_group <- shiny::renderUI({
    cat("output$pheno_group\n", file = stderr())
    shiny::req(choices <- pheno_group())
    if(is.null(selected <- input$pheno_group)) {
      selected <- choices[1]
    }
    shiny::selectInput(ns("pheno_group"), "",
                       choices = as.list(choices),
                       selected = selected,
                       multiple = TRUE)
  })
  output$dataset <- shiny::renderUI({
    cat("output$dataset\n", file = stderr())
    shiny::req(project_info())
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
                         input, pheno_type, peaks_tbl, pmap_obj, 
                         project_info)
  
  chr_pos <- shiny::reactive({
    shiny::req(project_info())
    make_chr_pos(win_par$chr_id, 
                 win_par$peak_Mbp, win_par$window_Mbp)
  })
  output$chr_pos <- shiny::renderText({
    paste0("Region: ", chr_pos(), "Mbp")
  })
  output$num_pheno <- shiny::renderText({
    shiny::req(project_info())
    num_pheno(character(), analyses_tbl())
  })
  shiny::observeEvent(project_info(), {
    cat("observeEvent project_info\n", file = stderr())
    output$num_pheno <- shiny::renderText({
      num_pheno(phe_par$pheno_names, analyses_tbl())
    })
  })
  shiny::observeEvent(phe_par$pheno_names, {
    cat("observeEvent phe_par$pheno_names\n", file = stderr())
    output$num_pheno <- shiny::renderText({
      num_pheno(phe_par$pheno_names, analyses_tbl())
    })
  })
  
  ## Use window as input to shinyPhenos.
  phe_par <- shiny::callModule(shinyPhenos, "phenos",
             input, win_par, peaks_tbl, analyses_tbl, pheno_data, cov_df,
             project_info)
  
  ## Setup input logic.
  output$title <- shiny::renderUI({
    shiny::tagList(
      switch(shiny::req(input$radio),
             Region = shinyProjectUI(ns("project")),
             Phenotypes = shiny::strong(paste("Project:", 
                                              shiny::req(project_info()$project),
                                              "\n"))),
      shiny::strong(shiny::req(input$radio)))
  })
  output$sidebar_setup <- shiny::renderUI({
    switch(shiny::req(input$radio),
           Phenotypes = shinyPhenosUI(ns("phenos")),
           Region     = shinyPeaksInput(ns("shinypeaks")))
  })
  output$main_setup <- shiny::renderUI({
    switch(shiny::req(input$radio),
           Phenotypes = shinyPhenosOutput(ns("phenos")),
           Region     = shinyPeaksOutput(ns("shinypeaks")))
  })
  
  output$radio <- shiny::renderUI({
    shiny::radioButtons(ns("radio"), NULL,
                 c("Region", "Phenotypes"),
                 input$radio,
                 inline=TRUE)
  })
  
  ## Return.
  shiny::reactive({
    list(project_info = project_info(),
         phe_par = phe_par,
         win_par = win_par)
  })
}
#' @param id identifier for \code{\link{shinyScan1}}
#' @rdname shinySetup
#' @export
shinySetupInput <- function(id) {
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
