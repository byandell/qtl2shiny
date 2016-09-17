#' Shiny SNP Association analysis and plot module
#'
#' Shiny module for SNP association analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,pheno_anal,probs_obj,K_chr,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' 
# Still work in progress. Idea is that
# server: callModule's Scan1SNP, TopSNP, SNPCsq, TopFeature
# ui screens:
#    1 Scan1SNP, SNPCsq
#    2 TopFeature, TopSNP
# Need to figure out download stuff as well.
## Need to worry about what extra options are for each item.
## Need to worry about download.
shinySNPAssoc <- function(input, output, session,
                          win_par, phe_df, cov_mx,
                          pheno_anal, probs_obj, K_chr,
                          ui_task,
                          snp_action = reactive({"basic"})) {
  ns <- session$ns

  chr_pos <- reactive({
    make_chr_pos(req(win_par()$chr_id), range = req(input$scan_window))
  })
  
  ## Reactives
  ## SNP Probabilities.
  snpprobs_obj <- reactive({
    withProgress(message = 'SNP Probs ...', value = 0, {
      setProgress(1)
      get_snpprobs(chr_id(), 
                   win_par()$peak_Mbp, 
                   win_par()$window_Mbp,
                   names(phe_df()), 
                   probs_obj(),
                   datapath)
    })
  })
  snpprobs_act <- reactive({
    snpprob_collapse(req(snpprobs_obj()), snp_action())
  })
  ## SNP Scan.
  snp_scan_obj <- reactive({
    req(pheno_anal())
    withProgress(message = "SNP Scan ...", value = 0, {
      setProgress(1)
      scan1(snpprobs_act(), phe_df(), K_chr(), cov_mx())
    })
  })
  ## Top SNPs table (used in several places)
  top_snps_tbl <- reactive({
    req(snp_action())
    withProgress(message = 'Get Top SNPs ...', value = 0, {
      setProgress(1)
      get_top_snps_tbl(req(snp_scan_obj()))
    })
  })
  ## Genes and Exons.
  gene_exon_tbl <- reactive({
    req(snp_action())
    withProgress(message = 'Gene Exon Calc ...', value = 0, {
      setProgress(1)
      get_gene_exon_snp(req(top_snps_tbl()),
                        file.path(datapath, "mgi_db.sqlite"))
    })
  })
  
  ## Shiny Modules
  ## SNP Association
  ## SNP Association Scan
  callModule(shinyScan1SNP, "snp_scan",
             input, chr_pos, snp_scan_obj, snp_action)
  ## SNP Summary
  callModule(shinySNP, "best_snp", 
             chr_pos, top_snps_tbl)
  ## SNP Consequence: Genes & Exons & Summary Table
  callModule(shinySNPCsq, "snp_csq", 
             input, chr_pos, top_snps_tbl, snp_action)

  ## Allele Patterns
  ## Top SNPs = Allele Patterns
  callModule(shinyTopSNP, "top_snp",
             input, chr_pos, snp_scan_obj, snp_action)
  ## Top Feature
  callModule(shinyTopFeature, "top_feature",
             chr_pos, snp_scan_obj, top_snps_tbl, gene_exon_tbl)

  # Scan Window slider
  output$scan_window <- renderUI({
    req(win_par())
    rng <- round(win_par()$peak_Mbp + 
                   c(-1,1) * win_par()$window_Mbp, 
                 1)
    if(is.null(selected <- input$scan_window))
      selected <- rng
    sliderInput(ns("scan_window"), NULL, rng[1], rng[2],
                selected, step=.1)
  })

  ## Select phenotype for plots.
  output$pheno_assoc <- renderUI({
    req(pheno_anal())
    selectInput(ns("pheno_assoc"), NULL,
                choices = names(pheno_anal()),
                selected = input$pheno_assoc)
  })
  
  ## Button Options.
  output$phe_choice <- renderUI({
    ## Show Phenotype Choice for Scan, Pattern, Top SNPs
    if((ui_task() == "SNP Association" & 
        input$snp_assoc %in% c("Scan")) |
       (ui_task() == "Allele Pattern" &
        input$allele_pat %in% c("Pattern","Top SNPs"))) {
      uiOutput(ns("pheno_assoc"))
    }
  })
  output$win_choice <- renderUI({
    ## Show Window for Scan, Genes, Pattern, Alls
    if((ui_task() == "SNP Association" & 
        input$snp_assoc %in% c("Scan","Genes")) |
       (ui_task() == "Allele Pattern" &
        input$allele_pat %in% c("Pattern","All Phenos","All Patterns"))) {
      uiOutput(ns("scan_window"))
    }
  })

  observeEvent(input$pheno_assoc, {
    button_val <- c("Top SNPs","Pattern",
                    "All Phenos","All Patterns",
                    "Summary")
    if(length(input$pheno_assoc == 1))) {
      button_val <- button_val[-(3:4)]
    }
    updateRadioButtons(session, "snp_assoc", 
                       selected = button_val[1],
                       choices = button_val)
  })
  output$snp_choice <- renderUI({
    switch(req(ui_task()),
           "SNP Association" = {
             tagList(
               radioButtons(ns("snp_assoc"), "",
                            c("Scan", "Genes", "Exons", "Summary")))
           },
           "Allele Pattern" = {
             tagList(
               radioButtons(ns("allele_pat"), "",
                            c("Top SNPs","Pattern",
                              "All Phenos","All Patterns",
                              "Summary")))
#               shinyTopSNPUI(ns("top_snp")))
           })
  })
  output$extra_choice <- renderUI({
    switch(req(ui_task()),
           "SNP Association" = {
             switch(input$snp_assoc,
                    Summary = shinySNPInput(ns("best_snp")))
           },
           "Allele Pattern" = {
             switch(input$allele_pat,
                    "Top SNPs" = shinyTopFeatureUI(ns("top_feature")))
           })
  })
  output$snp_output <- renderUI({
    switch(req(ui_task()),
           "SNP Association" = {
             switch(input$snp_assoc,
                    Scan = shinyScan1SNPOutput(ns("snp_scan")),
                    Summary = shinySNPOutput(ns("best_snp")),
                    shinySNPCsqOutput(ns("snp_csq")))
             },
           "Allele Pattern" = {
              switch(input$allele_pat,
                     "Top SNPs" = shinyTopFeatureOutput(ns("top_feature")),
                     shinyTopSNPOutput(ns("top_snp")))
             })
  })
  output$title <- renderUI({
   if(snp_action() == "basic")
      h4(strong(req(ui_task())))
  })

  ## Downloads
  output$download_snp_assoc <- renderUI({
    if(input$snp_assoc == "Scan") {
      shinyScan1SNPUI(ns("snp_scan"))
    } else {
      shinySNPCsqUI(ns("snp_csq"))
    }
  })
  output$download_csv_plot <- renderUI({
    switch(req(ui_task()),
           "SNP Association" = {
             fluidRow(
               column(6, shinySNPUI(ns("best_snp"))),
               column(6, uiOutput(ns("download_snp_assoc"))))
           },
           "Allele Pattern" = {
#             fluidRow(
#               column(6, shinyTopSNPUI(ns("best_snp"))),
#               column(6, shinyScan1SNPUI(ns("snpPlot"))))
           })
  })
  ########################################
  ## Not used here. Where to put?
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_effects_", chr_pos(), ".csv")) },
    content = function(file) {
      write.csv(sum_top_pat(), file)
    }
  )
  ########################################
  
  ## Return patterns
  reactive({ # patterns
    summary(top_snps_tbl()) %>%
      filter(max_lod >= 3) %>%
      mutate(contrast = snp_action()) %>%
      arrange(desc(max_lod))
  })
}
#' @param id identifier for \code{\link{shinySNPAssoc}} use
#' @rdname shinySNPAssoc
#' @export
shinySNPAssocUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("title")),
    uiOutput(ns("snp_choice")),
    uiOutput(ns("phe_choice")),
    uiOutput(ns("win_choice")),
    uiOutput(ns("extra_choice")),
    uiOutput(ns("download_csv_plot")))
}
#' @rdname shinySNPAssoc
#' @export
shinySNPAssocOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_output"))
}
