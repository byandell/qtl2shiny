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
  
  chr_pos <- reactive({
    make_chr_pos (win_par$chr_id, win_par$peak_Mbp, win_par$window_Mbp)
  })
  
  patterns <- callModule(shinySNPCsq, "dip_csq",
                             snp_scan_obj, chr_pos, 
                             snp_action)
  
  callModule(shinyPattern, "dip_pat", 
             probs1, patterns,
             phe_df, K_chr, cov_mx, chr_pos,
             snp_action)
  
  
  snp_action <- reactive({input$snp_action})
  output$snp_dip <- renderUI({
    switch(input$snp_dip,
           Scan = shinyScan1SNPUI(ns("dip_scan")),
           Consequence  = shinySNPCsqUI(ns("dip_csq")),
           Pattern = shinyPatternUI(ns("dip_pat")))
  })
}
#' @rdname shinyDominance
#' @export
shinyDominanceUI <- function(id) {
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
                         c("Scan","Consequence","Pattern")))),
    uiOutput(ns("snp_dip"))
  )
}
