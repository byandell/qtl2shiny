#' Shiny Diplotype module
#'
#' Shiny diplotype SNP/Gene action analysis.
#' 
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,K_chr,data_path reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyDiplo <- function(input, output, session,
                       win_par, 
                       phe_df, cov_mx, K_chr,
                       data_path) {
  ns <- session$ns

  chr_pos <- reactive({
    make_chr_pos(win_par$chr_id, 
                 win_par$peak_Mbp, 
                 win_par$window_Mbp)
  })

  ## Probs object for 36 diplotypes.
  probs36_obj <- callModule(shinyProbs36, "probs36", 
                            win_par, data_path)

  snp_action <- reactive({input$snp_action})
  
  ## SNP Association
  patterns <- callModule(shinySNPAllele, "snp_allele",
                         input, win_par, phe_df, cov_mx, 
                         probs36_obj, K_chr, data_path,
                         snp_action)
  
  callModule(shinyPattern, "dip_pat",
             chr_pos, phe_df, K_chr, cov_mx,
             probs36_obj, patterns,
             snp_action)

  ## CC names
  output$cc_names <- renderText({
    cc <- CCcolors
    paste(LETTERS[seq_along(cc)], names(cc), 
          sep = "=", collapse = ", ")
  })
  output$dip_input <- renderUI({
    switch(req(input$button),
           "Genome Scans"    = shinyPatternUI(ns("dip_pat")),
           "SNP Association" =,
           "Allele Pattern"  = shinySNPAlleleUI(ns("snp_allele")))
  })
  output$dip_output <- renderUI({
    switch(req(input$button),
           "Genome Scans"    = shinyPatternOutput(ns("dip_pat")),
           "SNP Association" = ,
           "Allele Pattern"  = shinySNPAlleleOutput(ns("snp_allele")))
  })
  output$radio <- renderUI({
    radioButtons(ns("button"), "",
                 c("SNP Association","Allele Pattern","Genome Scans"),
                 input$button)
  })
  output$select <- renderUI({
    selectInput(ns("snp_action"), "",
                c("add+dom","additive","non-add",
                  "recessive","dominant"),
                input$select)
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
      uiOutput(ns("radio")),
      uiOutput(ns("select")),
      uiOutput(ns("dip_input")),
      textOutput(ns("cc_names"))),
  mainPanel(
    uiOutput(ns("dip_output")))
  )
}
