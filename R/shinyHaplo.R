#' Shiny haplotype analysis
#'
#' Shiny module for analysis based on haplotype alleles.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,pheno_anal,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyHaplo <- function(input, output, session,
                       win_par, phe_df, cov_mx, pheno_anal, K_chr) {
  ns <- session$ns

  chr_id <- reactive({as.character(req(win_par()$chr_id))})
  
  ## Genotype Probabilities.
  probs_obj <- reactive({
    req(chr_id())
    withProgress(message = 'Read probs ...', value = 0, {
      setProgress(1)
      read_probs(chr_id(), datapath)
    })
  })
  
  ## Genome Scan.
  callModule(shinyScan1Plot, "hap_scan", 
             chr_id, phe_df, cov_mx,
             pheno_anal, probs_obj, K_chr)
  
  ## SNP Scan.
  snp_scan_obj <- callModule(shinyScan1SNP, "snp_scan",
                             win_par, phe_df, cov_mx, 
                             pheno_anal, probs_obj, K_chr)
  
  ## SNP Consequence.
  callModule(shinySNPCsq, "snp_csq", 
             win_par, snp_scan_obj)
  
  ## CC names
  output$cc_names <- renderText({
    cc <- CCcolors
    paste(LETTERS[seq_along(cc)], names(cc), sep = "=", collapse = ", ")
  })

  output$hap_choice <- renderUI({
    switch(input$snp_hap,
           "Genome Scans"    = shinyScan1PlotUI(ns("hap_scan")),
           "SNP Association" = shinyScan1SNPUI(ns("snp_scan")),
           "Allele Pattern"  = shinySNPCsqUI(ns("snp_csq")))
  })
  output$snp_hap <- renderUI({
    switch(input$snp_hap,
           "Genome Scans"    = shinyScan1PlotOutput(ns("hap_scan")),
           "SNP Association" = shinyScan1SNPOutput(ns("snp_scan")),
           "Allele Pattern"  = shinySNPCsqOutput(ns("snp_csq")))
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
      radioButtons(ns("snp_hap"), "",
                   c("Genome Scans","SNP Association","Allele Pattern")),
      uiOutput(ns("hap_choice")),
      textOutput(ns("cc_names"))),
    mainPanel(
      uiOutput(ns("snp_hap")))
  )
}
