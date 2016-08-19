#' Shiny Dominance module
#'
#' @param input,output,session standard shiny arguments
#' @param plot_type,chr_id,phe_df,cov_mx,pheno_anal,probs_obj,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyDominance <- function(input, output, session,
                       win_par, phe_df, cov_mx,
                       pheno_anal, probs_obj, K_chr) {
  
  ## This is basically obsolete now.
  
  ns <- session$ns
  chr_id <- reactive({win_par$chr_id})
  range_val <- reactive({
    req(win_par$peak_Mbp, win_par$window_Mbp)
    c(win_par$peak_Mbp + c(-1,1) * win_par$window_Mbp)
  })
  
  ## Probs object for 36 diplotypes.
  probs1 <- reactive({
    req(chr_id())
    withProgress(message = 'Diplotype Probs ...', value = 0, {
      setProgress(1)
      read_probs36(chr_id(), range_val()[1], range_val()[2],
                   file.path(datapath, "wave4", "DOQTL"))
    })
  })
  
  snp_scan_obj <- callModule(shinyScan1SNP, "dip_scan",
                             win_par, phe_df, cov_mx, 
                             pheno_anal, probs1, K_chr,
                             snp_action)
  
  ## Dominance type.
  output$snp_action <- renderUI({
    ## Need to modify snpprobs based on this.
    ## Cannot do at this level.
    selectInput(ns("snp_action"), "", 
                c("both","additive","dominance",
                  "B6-recessive","B6-dominance"))
  })
  snp_action <- reactive({input$snp_action})
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
}
#' @rdname shinyDominance
#' @export
shinyDominanceUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(6,h4(strong("SNP Gene Action"))),
      column(6, uiOutput(ns("snp_action")))),
    shinyScan1SNPUI(ns("dip_scan"))
  )
}
