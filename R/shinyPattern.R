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
  ns <- session$ns

  ## Names of haplos and diplos in terms of founders.
  founders <- reactive({
    founders <- str_split("AB1NZCPW","")[[1]]
    names(founders) <- LETTERS[seq_along(founders)]
    founders
  })
  haplos <- reactive({
    str_replace_all(dimnames(probs_obj()$probs[[1]])[[2]], founders())
  })
  diplos <- reactive({
    str_replace_all(dimnames(probs1()$probs[[1]])[[2]], founders())
  })
  
  patterns <- reactive({
    summary(top_snps_tbl()) %>%
      mutate(pheno = names(phe_df())[1])
  })
  
  output$scan_pattern <- renderPlot({
    ## Want to use idea of plot_snpscan except:
    ## split out LOD and effect
    ## get rid of pheno column use
    ## link up to shinyScan1Plot
    plot_snpscan(probs1(), phe_df(), K_chr(), cov_mx(),
                 patterns(),
                 haplos(), diplos())
  })
}
#' @rdname shinyPattern
#' @export
shinyPatternUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(3,h4(strong("SNP Gene Action"))),
      column(6, 
             selectInput(ns("snp_action"), "",
                         c("both","additive","dominance",
                           "B6-recessive","B6-dominance"))),
      column(3, 
             selectInput(ns("snp_dip"), "",
                         c("Scan","Consequence")))),
    uiOutput(ns("snp_dip"))
  )
}
