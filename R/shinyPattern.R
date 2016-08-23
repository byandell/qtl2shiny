#' Shiny Pattern module
#'
#' @param input,output,session standard shiny arguments
#' @param plot_type,chr_id,phe_df,cov_mx,pheno_anal,probs_obj,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyPattern <- function(input, output, session,
                         probs1, top_snps_tbl,
                         phe_df, K_chr, cov_mx) {
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
    req(input$pheno)
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
    founders <- str_split("AB1NZCPW","")[[1]]
    names(founders) <- LETTERS[seq_along(founders)]
    founders
  })
  diplos <- reactive({
    str_replace_all(dimnames(probs1()$probs[[1]])[[2]], haplos())
  })

  patterns <- reactive({
    summary(top_snps_tbl())
  })
  scan_pat <- reactive({
    pheno_in <- req(input$pheno)
    pattern_in <- req(input$pattern)
    pats <- patterns() %>%
      filter(AB1NZCPW == pattern_in,
             pheno == pheno_in)
    scan_pattern(probs1(),
                 phe_df()[,pheno_in, drop=FALSE],
                 K_chr(), cov_mx(),
                 pats,
                 haplos(), diplos())
  })

  output$scan_pat_lod <- renderPlot({
    plot(scan_pat(), "lod") + 
      theme(legend.position="none")
  })
  output$scan_pat_coef <- renderPlot({
    plot(scan_pat(), "coef")
  })
}
#' @rdname shinyPattern
#' @export
shinyPatternUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(6, uiOutput(ns("pheno"))),
      column(6, uiOutput(ns("pattern")))),
    plotOutput(ns("scan_pat_lod")),
    plotOutput(ns("scan_pat_coef"))
  )
}
