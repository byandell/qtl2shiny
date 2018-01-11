#' Shiny haplotype analysis
#'
#' Shiny module for analysis based on haplotype alleles.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,pmap_obj,phe_mx,cov_df,K_chr,analyses_df,covar,pheno_data,analyses_tbl,peaks,project_info,allele_names reactive arguments
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
                       phe_mx, cov_df, K_chr, analyses_df, 
                       covar, pheno_data, analyses_tbl, peaks,
                       project_info, allele_names) {
  ns <- session$ns

  ## Genotype Probabilities.
  probs_obj <- shiny::callModule(shinyProbs, "probs", 
                          win_par, project_info)

  ## Genome Scan.
  shiny::callModule(shinyScan1Plot, "hap_scan", 
             input, win_par, 
             phe_mx, cov_df, probs_obj, K_chr, analyses_df,
             project_info)
  
  ## SNP Association
  patterns <- shiny::callModule(shinySNPAllele, "snp_allele",
              input, win_par, 
              phe_mx, cov_df, probs_obj, K_chr, analyses_df,
              project_info)

  ## Mediation
  shiny::callModule(shinyMediate1Plot, "mediate",
                    input, win_par, patterns,
                    phe_mx, cov_df, probs_obj, K_chr, analyses_df,
                    pmap_obj, 
                    covar, pheno_data, analyses_tbl, peaks,
                    project_info)

  output$allele_names <- shiny::renderText({
    shiny::req(allele_names())
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
    if(ncol(shiny::req(phe_mx())) > 1 & shiny::req(input$button) != "Mediation") {
      choices <- choices[1:4]
    }
    shiny::radioButtons(ns("sex_type"), "Sex:",
                        choices,
                        input$sex_type, inline = TRUE)
  })
  output$project <- shiny::renderUI({
    shiny::strong(shiny::req(paste("Project:",
                                   project_info()$project,
                                   "\n")))
  })
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyHaplo
#' @export
shinyHaploUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::sidebarPanel(
      shiny::uiOutput(ns("project")),
      shiny::strong("SNP/Gene Additive"),
      shiny::uiOutput(ns("radio")),
      shiny::uiOutput(ns("sex_type")),
      shiny::uiOutput(ns("hap_input")),
      shiny::textOutput(ns("allele_names"))),
    shiny::mainPanel(
      shiny::uiOutput(ns("hap_output")))
  )
}
