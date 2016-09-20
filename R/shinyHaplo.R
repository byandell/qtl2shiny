#' Shiny haplotype analysis
#'
#' Shiny module for analysis based on haplotype alleles.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyHaplo <- function(input, output, session,
                       win_par, phe_df, cov_mx, K_chr) {
  ns <- session$ns

  chr_pos <- reactive({
    make_chr_pos(win_par$chr_id, 
                 win_par$peak_Mbp, win_par$window_Mbp)
  })

  ## Genotype Probabilities.
  probs_obj <- callModule(shinyProbs, "probs", 
                          win_par)

  ## Genome Scan.
  callModule(shinyScan1Plot, "hap_scan", 
             chr_pos, phe_df, cov_mx, probs_obj, K_chr)
  
  ## SNP Association
  patterns <- callModule(shinySNPAllele, "snp_allele",
              input, win_par, phe_df, cov_mx, probs_obj, K_chr)

  ## CC names
  output$cc_names <- renderText({
    cc <- CCcolors
    paste(LETTERS[seq_along(cc)], names(cc), sep = "=", collapse = ", ")
  })

  output$hap_input <- renderUI({
    switch(req(input$button),
           "Genome Scans"    = shinyScan1PlotUI(ns("hap_scan")),
           "SNP Association" =,
           "Allele Pattern"  = shinySNPAlleleUI(ns("snp_allele")))
  })
  output$hap_output <- renderUI({
    switch(req(input$button),
           "Genome Scans"    = shinyScan1PlotOutput(ns("hap_scan")),
           "SNP Association" = ,
           "Allele Pattern"  = shinySNPAlleleOutput(ns("snp_allele")))
  })
  output$radio <- renderUI({
    radioButtons(ns("button"), "",
                 c("Genome Scans","SNP Association","Allele Pattern"),
                 input$button)
  })
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyHaplo
#' @export
shinyHaploUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarPanel(
      h4(strong("SNP/Gene Haplo Analysis")),
      uiOutput(ns("radio")),
      uiOutput(ns("hap_input")),
      textOutput(ns("cc_names"))),
    mainPanel(
      uiOutput(ns("hap_output")))
  )
}
