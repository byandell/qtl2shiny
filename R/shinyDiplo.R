#' Shiny Diplotype module
#'
#' Shiny diplotype SNP/Gene action analysis.
#' 
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyDiplo <- function(input, output, session,
                       win_par, phe_df, cov_mx, K_chr) {
  ns <- session$ns

  chr_id <- reactive({win_par()$chr_id})
  chr_pos <- reactive({
    make_chr_pos (win_par()$chr_id, win_par()$peak_Mbp, win_par()$window_Mbp)
  })

  range_val <- reactive({
    req(win_par()$peak_Mbp, win_par()$window_Mbp)
    c(win_par()$peak_Mbp + c(-1,1) * win_par()$window_Mbp)
  })

  ## Probs object for 36 diplotypes.
  probs36_obj <- reactive({
    req(chr_id())
    withProgress(message = 'Diplotype Probs ...', value = 0, {
      setProgress(1)
      read_probs36(chr_id(), range_val()[1], range_val()[2],
                   datapath)
    })
  })

  snp_action <- reactive({input$snp_action})
  
  ## SNP Association
  patterns <- callModule(shinySNPAllele, "snp_allele",
                         input, win_par, phe_df, cov_mx, 
                         probs36_obj, K_chr, snp_action)
  
  callModule(shinyPattern, "dip_pat",
             chr_pos, phe_df, K_chr, cov_mx,
             probs36_obj, patterns,
             snp_action)

  ## CC names
  output$cc_names <- renderText({
    cc <- CCcolors
    paste(LETTERS[seq_along(cc)], names(cc), sep = "=", collapse = ", ")
  })
  output$dip_input <- renderUI({
    switch(input$button,
           "SNP Plots"       = shinyScan1SNPUI(ns("dip_scan")),
           "SNP Association" =,
           "Allele Pattern"  = shinySNPAlleleUI(ns("snp_allele")))
  })
  output$dip_output <- renderUI({
    switch(input$button,
           "SNP Plots"       = shinyScan1SNPOutput(ns("dip_scan")),
           "SNP Association" = ,
           "Allele Pattern"  = shinySNPAlleleOutput(ns("snp_allele")))
  })
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyDiplo
#' @export
shinyDiploUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarPanel(
      h4(strong("SNP/Gene Action")),
      radioButtons(ns("button"), "",
                   c("SNP Plots","Consequence","Genome Scans")),
      selectInput(ns("snp_action"), "",
                  c("add+dom","additive","non-add",
                    "recessive","dominant")),
      uiOutput(ns("dip_input")),
      textOutput(ns("cc_names"))),
  mainPanel(
    uiOutput(ns("dip_output")))
  )
}
