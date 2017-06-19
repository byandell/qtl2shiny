#' Shiny haplotype analysis
#'
#' Shiny module for analysis based on haplotype alleles.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,pmap_obj,phe_df,cov_mx,K_chr,analyses_df,data_path reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom shiny callModule NS req 
#'   radioButtons
#'   textOutput uiOutput
#'   renderText renderUI
#'   mainPanel sidebarPanel strong tagList
shinyHaplo <- function(input, output, session,
                       win_par, pmap_obj, 
                       phe_df, cov_mx, K_chr, analyses_df,
                       data_path) {
  ns <- session$ns

  ## Genotype Probabilities.
  probs_obj <- shiny::callModule(shinyProbs, "probs", 
                          win_par, data_path)

  ## Genome Scan.
  shiny::callModule(shinyScan1Plot, "hap_scan", 
             input, win_par, 
             phe_df, cov_mx, probs_obj, K_chr, analyses_df)
  
  ## SNP Association
  patterns <- shiny::callModule(shinySNPAllele, "snp_allele",
              input, win_par, 
              phe_df, cov_mx, probs_obj, K_chr, analyses_df,
              data_path)

  ## Mediation
  shiny::callModule(shinyMediate1Plot, "mediate",
                    input, win_par, patterns,
                    phe_df, cov_mx, probs_obj, K_chr, analyses_df,
                    pmap_obj, data_path)

  ## CC names
  output$cc_names <- shiny::renderText({
    cc <- CCSanger::CCcolors
    paste(LETTERS[seq_along(cc)], names(cc), sep = "=", collapse = ", ")
  })

  output$hap_input <- shiny::renderUI({
    switch(shiny::req(input$button),
           "Genome Scans"    = shinyScan1PlotUI(ns("hap_scan")),
           "SNP Association" =,
           "Allele Pattern"  = shinySNPAlleleUI(ns("snp_allele")),
           "Mediation"       = shinyMediate1PlotUI(ns("mediate")))
  })
  output$hap_output <- shiny::renderUI({
    switch(shiny::req(input$button),
           "Genome Scans"    = shinyScan1PlotOutput(ns("hap_scan")),
           "SNP Association" = ,
           "Allele Pattern"  = shinySNPAlleleOutput(ns("snp_allele")),
           "Mediation"       = shinyMediate1PlotOutput(ns("mediate")))
  })
  output$radio <- shiny::renderUI({
    shiny::radioButtons(ns("button"), "",
                 c("Genome Scans","SNP Association","Allele Pattern","Mediation"),
                 input$button)
  })
  output$sex_type <- shiny::renderUI({
    choices <- c("A","I","F","M","all")
    if(ncol(shiny::req(phe_df())) > 1 & shiny::req(input$button) != "Mediation") {
      choices <- choices[1:4]
    }
    shiny::radioButtons(ns("sex_type"), "Sex:",
                        choices,
                        input$sex_type, inline = TRUE)
  })
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyHaplo
#' @export
shinyHaploUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::sidebarPanel(
      shiny::strong("SNP/Gene Additive"),
      shiny::uiOutput(ns("radio")),
      shiny::uiOutput(ns("sex_type")),
      shiny::uiOutput(ns("hap_input")),
      shiny::textOutput(ns("cc_names"))),
    shiny::mainPanel(
      shiny::uiOutput(ns("hap_output")))
  )
}
