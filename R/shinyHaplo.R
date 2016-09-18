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
  chr_pos <- reactive({
    make_chr_pos(win_par()$chr_id, 
                 win_par()$peak_Mbp, win_par()$window_Mbp)
  })

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
  
  ## SNP Association
  top_snps_tbl <- callModule(shinySNPAllele, "snp_assoc",
                             win_par, phe_df, cov_mx, 
                             pheno_anal, probs_obj, K_chr)
  
  
  
  ## CC names
  output$cc_names <- renderText({
    cc <- CCcolors
    paste(LETTERS[seq_along(cc)], names(cc), sep = "=", collapse = ", ")
  })

  output$hap_choice <- renderUI({
    switch(input$snp_hap,
           "Genome Scans"    = shinyScan1PlotUI(ns("hap_scan")),
           "SNP Association" = tagList(
             shinyScan1SNPUI(ns("snp_scan")),
             shinySNPCsqUI(ns("snp_csq"))),
           "Allele Pattern"  = shinyTopSNPUI(ns("top_snp")))
  })
  output$snp_hap <- renderUI({
    switch(input$snp_hap,
           "Genome Scans"    = shinyScan1PlotOutput(ns("hap_scan")),
           "SNP Association" = tagList(
             shinyScan1SNPOutput(ns("snp_scan")),
             shinySNPCsqOutput(ns("snp_csq"))),
           "Allele Pattern"  = shinyTopSNPOutput(ns("top_snp")))
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
