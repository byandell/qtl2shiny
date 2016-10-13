#' Shiny Pattern module
#'
#' @param input,output,session standard shiny arguments
#' @param chr_pos,phe_df,K_chr,cov_mx,probs36_obj,patterns,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyPattern <- function(input, output, session,
                         chr_pos, phe_df, K_chr, cov_mx,
                         probs36_obj, patterns,
                         snp_action = reactive({NULL})) {
  ns <- session$ns
  
  ###### Want title to have phenotype name.

  ## Select phenotype for plots.
  output$pheno_name <- renderUI({
    selectInput(ns("pheno_name"), NULL,
                choices = names(phe_df()),
                selected = input$pheno_name)
  })
  ## Select pattern for plots.
  output$pattern <- renderUI({
    req(input$pheno_name,snp_action())
    pats <- patterns() %>%
      filter(pheno == input$pheno_name)
    if(nrow(pats)) {
      choices <- sdp_to_pattern(pats$sdp)
    } else {
      choices <- input$pattern
    }
    selectInput(ns("pattern"), NULL,
                choices = choices,
                selected = input$pattern)
  })

  ## Names of haplos and diplos in terms of founders.
  haplos <- reactive({
    unique(unlist(str_split(diplos(), "")))
  })
  diplos <- reactive({
    dimnames(probs36_obj()$probs[[1]])[[2]]
  })

  scan_pat <- reactive({
    pheno_in <- req(input$pheno_name)
    pats <- patterns()
    withProgress(message = 'Scan Patterns ...', value = 0, {
      setProgress(1)
      scan_pattern(probs36_obj(),
                   phe_df()[,pheno_in, drop=FALSE],
                   K_chr(), cov_mx(),
                   pats,
                   haplos(), diplos())
    })
  })

  scan_pat_type <- function(scan_pat, type, pattern) {
    pattern_cont <- 
      (scan_pat$patterns %>%
         filter(sdp_to_pattern(sdp) == pattern))$contrast
    plot(scan_pat, type, pattern_cont) 
  }

  output$scan_pat_lod <- renderPlot({
    req(scan_pat(),input$pattern)
    withProgress(message = 'Pattern LODs ...', value = 0, {
      setProgress(1)
      scan_pat_type(scan_pat(), "lod", input$pattern)
    })
  })
  output$scan_pat_coef <- renderPlot({
    req(scan_pat(),input$pattern)
    withProgress(message = 'Pattern Effects ...', value = 0, {
      setProgress(1)
      scan_pat_type(scan_pat(), "coef", input$pattern)
    })
  })
  output$scanSummary <- renderDataTable({
    req(scan_pat())
    withProgress(message = 'Pattern summary ...', value = 0, {
      setProgress(1)
      summary(scan_pat())
    })
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))

  output$scan_choice <- renderUI({
    switch(req(input$button),
           LOD =,
           Effects = uiOutput(ns("pattern")))
  })

  output$pat_output <- renderUI({
    switch(req(input$button),
           LOD     = plotOutput(ns("scan_pat_lod")),
           Effects = plotOutput(ns("scan_pat_coef")),
           Summary = dataTableOutput(ns("scanSummary")))
  })

  ## Downloads
  output$downloadData <- downloadHandler(
    filename = function() {
      pheno_in <- req(input$pheno_name)
      file.path(paste0(pheno_in, "_", snp_action(), "_effects_", chr_pos(), ".csv"))
    },
    content = function(file) {
      scan_in <- req(scan_pat())
      write.csv(summary(scan_in), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      pheno_in <- req(input$pheno_name)
      file.path(paste0(pheno_in, "_", snp_action(), "_scan_", chr_pos(), ".pdf"))
    },
    content = function(file) {
      scan_in <- req(scan_pat())
      pdf(file, width = 9)
      print(plot(scan_in, "lod"))
      plot(scan_in, "coef")
      dev.off()
    }
  )
  output$radio <- renderUI({
    radioButtons(ns("button"), "",
                 c("LOD","Effects","Summary"),
                 input$button)
  })
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyPattern
#' @export
shinyPatternUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("radio")),
    uiOutput(ns("pheno_name")),
    uiOutput(ns("scan_choice")),
    fluidRow(
      column(6, downloadButton(ns("downloadData"), "CSV")),
      column(6, downloadButton(ns("downloadPlot"), "Plots")))
  )
}
#' @rdname shinyPattern
#' @export
shinyPatternOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("pat_output"))
}
