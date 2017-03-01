#' Shiny Pattern module
#'
#' @param input,output,session standard shiny arguments
#' @param chr_pos,phe_df,cov_mx,probs36_obj,K_chr,analyses_df,patterns,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' 
#' @importFrom dplyr distinct
#' @importFrom stringr str_split
#' @importFrom grid plotViewport pushViewport
#' @importFrom gridBase baseViewports
#' @importFrom CCSanger sdp_to_pattern
#' @importFrom shiny NS reactive req 
#'   observeEvent
#'   radioButtons selectInput updateSelectInput
#'   dataTableOutput plotOutput uiOutput
#'   renderDataTable renderPlot renderUI
#'   fluidRow column tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
#'   
shinyPattern <- function(input, output, session,
                         chr_pos, 
                         phe_df, cov_mx, probs36_obj, K_chr, analyses_df,
                         patterns, snp_action = shiny::reactive({NULL})) {
  ns <- session$ns
  
  ###### Want title to have phenotype name.

  ## Select phenotype for plots.
  output$pheno_name <- shiny::renderUI({
    shiny::selectInput(ns("pheno_name"), NULL,
                choices = names(phe_df()),
                selected = input$pheno_name)
  })
  ## Select pattern for plots.
  pats <- shiny::reactive({
    shiny::req(input$pheno_name, patterns())
    dplyr::filter(patterns(), pheno == input$pheno_name)
  })
  pattern_choices <- shiny::reactive({
    CCSanger::sdp_to_pattern(pats()$sdp)
  })
  output$pattern <- shiny::renderUI({
    shiny::req(pattern_choices(), snp_action())
    choices <- pattern_choices()
    if(!length(choices)) {
      choices <- input$pattern
    }
    shiny::selectInput(ns("pattern"), NULL,
                choices = choices,
                selected = input$pattern)
  })
  shiny::observeEvent(patterns(), update_patterns())
  shiny::observeEvent(input$pheno_name, update_patterns())
  update_patterns <- function() {
    shiny::req(snp_action(), input$pheno_name, patterns())
    pats <- dplyr::filter(patterns(), pheno == input$pheno_name)
    if(nrow(pats)) {
      choices <- CCSanger::sdp_to_pattern(pats$sdp)
    } else {
      choices <- input$pattern
    }
    if(!is.null(selected <- input$pattern)) {
      if(!(selected %in% choices))
        selected <- choices[1]
    }
    shiny::updateSelectInput(session, "pattern", NULL,
                      choices, selected)
  }
  
  ## Names of haplos and diplos in terms of founders.
  haplos <- shiny::reactive({
    unique(unlist(stringr::str_split(diplos(), "")))
  })
  diplos <- shiny::reactive({
    dimnames(req(probs36_obj())$probs[[1]])[[2]]
  })

  scan_pat <- shiny::reactive({
    req(snp_action())
    pheno_in <- shiny::req(input$pheno_name)
    shiny::req(phe_df(), cov_mx(), probs36_obj(), K_chr(),
               patterns())
    withProgress(message = 'Scan Patterns ...', value = 0, {
      setProgress(1)
      scan1_pattern(pheno_in, phe_df(), cov_mx(), 
                    probs36_obj(), K_chr(), analyses_df(),
                                pats(), haplos(), diplos())
    })
  })

  output$scan_pat_lod <- shiny::renderPlot({
    shiny::req(scan_pat(), pattern_choices(), input$pheno_name)
    withProgress(message = 'Pattern LODs ...', value = 0, {
      setProgress(1)
      scan_pat_type(scan_pat(), "lod", pattern_choices(), input$pheno_name)
    })
  })
  output$scan_pat_coef <- shiny::renderPlot({
    shiny::req(scan_pat(), pattern_choices(), input$pheno_name)
    withProgress(message = 'Pattern Effects ...', value = 0, {
      setProgress(1)
      scan_pat_type(scan_pat(), "coef", pattern_choices(), input$pheno_name)
    })
  })
  output$scanSummary <- shiny::renderDataTable({
    shiny::req(scan_pat())
    withProgress(message = 'Pattern summary ...', value = 0, {
      setProgress(1)
      summary(scan_pat())
    })
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))

  output$eff_lodPlot <- shiny::renderPlot({
    shiny::req(pattern_choices(), input$pheno_name)
    withProgress(message = 'Pattern Effects & LOD ...', value = 0, {
      setProgress(1)
      scan_pat_type(scan_pat(), "coef_and_lod", pattern_choices(), input$pheno_name)
    })
  })
  
  output$pattern_choice <- shiny::renderUI({
    switch(shiny::req(input$button),
           LOD =,
           "LOD & Effects" =,
           Effects = shiny::uiOutput(ns("pattern")))
  })

  output$LOD <- shiny::renderUI({
    switch(shiny::req(input$button),
           LOD             = shiny::plotOutput(ns("scan_pat_lod")))
  })
  output$Effects <- shiny::renderUI({
    switch(shiny::req(input$button),
           Effects         = shiny::plotOutput(ns("scan_pat_coef")))
  })
  output$Both <- shiny::renderUI({ ## does not work yet
    ## needs work on qtl2ggplot::plot_coef_and_lod for listof_scan1coef object
    switch(shiny::req(input$button),
           "LOD & Effects" = shiny::plotOutput(ns("eff_lodPlot")))
  })
  output$Summary <- shiny::renderUI({
    switch(shiny::req(input$button),
           Summary = shiny::dataTableOutput(ns("scanSummary")))
  })

  ## Downloads
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      pheno_in <- shiny::req(input$pheno_name)
      file.path(paste0(pheno_in, "_", snp_action(), "_effects_", chr_pos(), ".csv"))
    },
    content = function(file) {
      scan_in <- shiny::req(scan_pat())
      write.csv(summary(scan_in), file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      shiny::req(input$pheno_name)
      file.path(paste0("pattern", "_", snp_action(), "_scan_", chr_pos(), ".pdf"))
    },
    content = function(file) {
      shiny::req(scan_pat(), pattern_choices())
      scan_in <- shiny::req(scan_pat())
      top_panel_prop <- 0.65
      pdf(file, width = 9, height = 9)
      for(pheno_in in names(phe_df())) {
        pats <- dplyr::filter(patterns(), pheno == pheno_in)
        pat_choices <- CCSanger::sdp_to_pattern(pats$sdp)
        
        scan_now <- scan1_pattern(pheno_in, phe_df(), cov_mx(), 
                                  probs36_obj(), K_chr(), analyses_df(),
                                  pats, haplos(), diplos())
        
        scan_pat_type(scan_now, "coef_and_lod", pat_choices, pheno_in)
      }
      dev.off()
    }
  )
  output$radio <- shiny::renderUI({
    shiny::radioButtons(ns("button"), "",
                 c("LOD","Effects","LOD & Effects","Summary"),
                 input$button)
  })
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyPattern
#' @export
shinyPatternUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("radio")),
    shiny::uiOutput(ns("pheno_name")),
#    shiny::uiOutput(ns("pattern_choice")),
    shiny::fluidRow(
      shiny::column(6, shiny::downloadButton(ns("downloadData"), "CSV")),
      shiny::column(6, shiny::downloadButton(ns("downloadPlot"), "Plots")))
  )
}
#' @rdname shinyPattern
#' @export
shinyPatternOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("LOD")),
    shiny::uiOutput(ns("Effects")),
    shiny::uiOutput(ns("Both")),
    shiny::uiOutput(ns("Summary"))
  )
}
