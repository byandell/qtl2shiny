#' Shiny Pattern module
#'
#' @param input,output,session standard shiny arguments
#' @param probs1,patterns,phe_df,K_chr,cov_mx,chr_pos,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyPattern <- function(input, output, session,
                         probs1, patterns,
                         phe_df, K_chr, cov_mx, chr_pos,
                         snp_action = reactive({NULL})) {
  ## lod values too small--getting wrong one?
  ## error with replacement having 1 row briefly when switch pheno
  ns <- session$ns
  
  ## Select phenotype for plots.
  output$pheno <- renderUI({
    selectInput(ns("pheno"), NULL,
                choices = names(phe_df()),
                selected = input$pheno)
  })
  ## Select pattern for plots.
  output$pattern <- renderUI({
    req(input$pheno,snp_action())
    pats <- patterns() %>%
      filter(pheno == input$pheno)
    if(nrow(pats)) {
      choices <- pats$AB1NZCPW
    } else {
      choices <- input$pattern
    }
    selectInput(ns("pattern"), NULL,
                choices = choices,
                selected = input$pattern)
  })

  ## Names of haplos and diplos in terms of founders.
  haplos <- reactive({
    str_split("AB1NZCPW","")[[1]]
  })
  diplos <- reactive({
    dimnames(probs1()$probs[[1]])[[2]]
  })

  scan_pat <- reactive({
    pheno_in <- req(input$pheno)
    pats <- patterns()
    withProgress(message = 'Scan Patterns ...', value = 0, {
      setProgress(1)
      scan_pattern(probs1(),
                   phe_df()[,pheno_in, drop=FALSE],
                   K_chr(), cov_mx(),
                   pats,
                   haplos(), diplos())
    })
  })
  
  scan_pat_lod <- reactive({
    req(scan_pat())
    p <- plot(scan_pat(), "lod")
    if(length(patterns()) == 1)
      p <- p + theme(legend.position="none")
    p
    
  })

  output$scan_pat_lod <- renderPlot({
    withProgress(message = 'Pattern LODs ...', value = 0, {
      setProgress(1)
      scan_pat_lod()
    })
  })
  output$scan_pat_coef <- renderPlot({
    scan_in <- req(scan_pat())
    pattern_in <- req(input$pattern)
    withProgress(message = 'Pattern Effects ...', value = 0, {
      setProgress(1)
      pattern_cont <- (scan_in$patterns %>%
                         filter(AB1NZCPW == pattern_in))$contrast
      plot(scan_in, "coef", pattern_cont)
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
  
  output$genome_scan <- renderUI({
    switch(req(input$genome_scan),
           LOD = {
             plotOutput(ns("scan_pat_lod"))
           },
           Effects = {
             tagList(
               uiOutput(ns("pattern")),
               plotOutput(ns("scan_pat_coef"))
             )
           },
           Summary = {
             dataTableOutput(ns("scanSummary"))
           })
  })
  
  ## Downloads
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("sum_effects_", chr_pos(), ".csv")) },
    content = function(file) {
      req(eff_obj(),scan_obj())
      write.csv(summary(scan_pat()), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("genome_scan_", chr_pos(), ".pdf")) },
    content = function(file) {
      pdf(file)
      scan_pat_lod()
      dev.off()
    }
  )
}
#' @rdname shinyPattern
#' @export
shinyPatternUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(3, 
             h4(strong("Genome Scans")),
             radioButtons(ns("genome_scan"), "",
                          c("LOD","Effects","Summary")),
             uiOutput(ns("pheno")),
             fluidRow(
               column(6, downloadButton(ns("downloadData"), "CSV")),
               column(6, downloadButton(ns("downloadPlot"), "Plots")))),
      column(9,
             uiOutput(ns("genome_scan"))))
  )
}
