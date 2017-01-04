#' Shiny Pattern module
#'
#' @param input,output,session standard shiny arguments
#' @param chr_pos,phe_df,K_chr,cov_mx,probs36_obj,patterns,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom dplyr distinct
#' @importFrom grid plotViewport pushViewport
#' @importFrom gridBase baseViewports
#' @importFrom qtl2ggplot sdp_to_pattern
#' @importFrom doqtl2 scan_pattern
#' @importFrom shiny NS reactive req 
#'   observeEvent
#'   radioButtons selectInput updateSelectInput
#'   dataTableOutput plotOutput uiOutput
#'   renderDataTable renderPlot renderUI
#'   fluidRow column tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
shinyPattern <- function(input, output, session,
                         chr_pos, phe_df, K_chr, cov_mx,
                         probs36_obj, patterns,
                         snp_action = shiny::reactive({NULL})) {
  ns <- session$ns
  
  ###### Want title to have phenotype name.

  ## Select phenotype for plots.
  output$pheno_name <- shiny::renderUI({
    shiny::selectInput(ns("pheno_name"), NULL,
                choices = names(phe_df()),
                selected = input$pheno_name)
  })
  ## Select pattern for plots.
  output$pattern <- shiny::renderUI({
    shiny::req(input$pheno_name, snp_action())
    pats <- dplyr::filter(patterns(), pheno == input$pheno_name)
    if(nrow(pats)) {
      choices <- sdp_to_pattern(pats$sdp)
    } else {
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
      choices <- qtl2ggplot::sdp_to_pattern(pats$sdp)
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
    unique(unlist(str_split(diplos(), "")))
  })
  diplos <- shiny::reactive({
    dimnames(probs36_obj()$probs[[1]])[[2]]
  })

  scan_pat <- shiny::reactive({
    pheno_in <- shiny::req(input$pheno_name)
    pats <- patterns()
    withProgress(message = 'Scan Patterns ...', value = 0, {
      setProgress(1)
      doqtl2::scan_pattern(probs36_obj(),
                   phe_df()[,pheno_in, drop=FALSE],
                   K_chr(), cov_mx(),
                   pats,
                   haplos(), diplos())
    })
  })

  scan_pat_type <- function(scan_pat, type, pattern) {
    pattern_cont <- 
      dplyr::filter(scan_pat$patterns,
                    qtl2ggplot::sdp_to_pattern(sdp) == pattern)$contrast[1]
    plot(scan_pat, type, pattern_cont) 
  }

  output$scan_pat_lod <- shiny::renderPlot({
    shiny::req(scan_pat(), input$pattern, input$pheno_name)
    withProgress(message = 'Pattern LODs ...', value = 0, {
      setProgress(1)
      scan_pat_type(scan_pat(), "lod", input$pattern)
    })
  })
  output$scan_pat_coef <- shiny::renderPlot({
    shiny::req(scan_pat(), input$pattern, input$pheno_name)
    withProgress(message = 'Pattern Effects ...', value = 0, {
      setProgress(1)
      scan_pat_type(scan_pat(), "coef", input$pattern)
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
    shiny::req(scan_pat(), input$pattern)
    withProgress(message = 'Pattern Effects & LOD ...', value = 0, {
      setProgress(1)
      par(mfrow=c(2,1))
      scan_pat_type(scan_pat(), "coef", input$pattern)
      plot.new()
      vps <- gridBase::baseViewports()
      grid::pushViewport(vps$figure)
      vp1 <- grid::plotViewport(c(1,1,1,1)) 
      print(scan_pat_type(scan_pat(), "lod", input$pattern),
            vp = vp1)
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
           LOD             =,
           "LOD & Effects" = shiny::plotOutput(ns("scan_pat_lod")))
  })
  output$Effects <- shiny::renderUI({
    switch(shiny::req(input$button),
           Effects         =,
           "LOD & Effects" = shiny::plotOutput(ns("scan_pat_coef")))
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
      pheno_in <- shiny::req(input$pheno_name)
      file.path(paste0(pheno_in, "_", snp_action(), "_scan_", chr_pos(), ".pdf"))
    },
    content = function(file) {
      scan_in <- shiny::req(scan_pat())
      pats <- dplyr::filter(patterns(), pheno == input$pheno_name)
      if(nrow(pats)) {
        pdf(file, width = 9, height = 9)
        choices <- qtl2ggplot::sdp_to_pattern(pats$sdp)
        for(pattern in choices) {
          par(mfrow=c(2,1))
          scan_pat_type(scan_in, "coef", pattern)
          plot.new()
          if(pattern == choices[1]) {
            vps <- gridBase::baseViewports()
            grid::pushViewport(vps$figure)
            vp1 <- grid::plotViewport(c(0,0,0,0)) 
          }
          print(scan_pat_type(scan_in, "lod", pattern),
                vp = vp1)
        }
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
    shiny::uiOutput(ns("pattern_choice")),
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
    shiny::uiOutput(ns("Summary"))
  )
}
