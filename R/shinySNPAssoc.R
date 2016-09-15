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
                          snp_action = reactive({"basic"})) {
  ns <- session$ns

  ## Scan Window needed here.
  
  chr_pos <- reactive({
    make_chr_pos(req(win_par()$chr_id), range = req(input$scan_window))
  })
  
  ## Shiny Modules
  ## SNP Association Scan
  snp_scan_obj <- callModule(shinyScan1SNP, "snp_scan",
                             win_par, input, phe_df, cov_mx, 
                             pheno_anal, probs_obj, K_chr,
                             snp_action)
  ## SNP Summary
  callModule(shinySNP, "best_snp", chr_pos, top_snps_tbl)
  ## Top SNPs table (used in several places)
  top_snps_tbl <- reactive({
    req(snp_action())
    withProgress(message = 'Get Top SNPs ...', value = 0, {
      setProgress(1)
      get_top_snps_tbl(snp_scan_obj())
    })
  })
  ## SNP Consequence: Genes & Exons & Summary Table
  gene_exon_tbl <- callModule(shinySNPCsq, "snp_csq", 
             chr_pos, top_snps_tbl, snp_action)

  ## Top SNPs = Allele Patterns
  callModule(shinyTopSNP, "top_snp",
             chr_pos, snp_scan_obj, snp_action)
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
  output$gene_choice <- renderUI({
    ## Show Options for Top SNPs
    if(ui_task() == "Allele Pattern" &
       input$allele_pat %in% c("Top SNPs")) {
      shinyTopFeatureUI(ns("top_feature"))
    }
    ## Show Options for Exons
    if(ui_task() == "SNP Association" &
       input$allele_pat %in% c("Exons")) {
      shinyTopFeatureUI(ns("top_feature"))
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
#               shinyScan1SNPUI(ns("snpPlot")),
#               shinySNPCsqUI(ns("snp_csq")))
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
  output$snp_scan <- renderUI({
    switch(req(ui_task()),
           "SNP Association" = {
             if(input$snp_assoc == "Scan") {
               shinyScan1SNPOutput(ns("snpPlot"))
             } else {
               shinySNPCsqOutput(ns("snp_csq"))
             }},
           "Allele Pattern" = {
             tagList(
               shinyTopFeatureOutput(ns("top_feature")),
               shinyTopSNPOutput(ns("top_snp")))
             })
  })
  output$title <- renderUI({
    if(snp_action() == "basic")
      h4(strong("SNP Association Plots"))
  })

  ## Downloads
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_effects_", chr_pos(), ".csv")) },
    content = function(file) {
      write.csv(sum_top_pat(), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_scan_", chr_pos(), ".pdf")) },
    content = function(file) {
      scans <- req(scan_obj())
      snp_w <- req(input$scan_window)
      phenos <- req(phename())
      pdf(file)
      ## Plots over all phenotypes
      print(top_pat_plot(phenos, scans, snp_w,
                         group = "pheno", snp_action = snp_action()))
      print(top_pat_plot(phenos, scans, snp_w,
                         group = "pattern", snp_action = snp_action()))
      ## Plots by phenotype.
      for(pheno in phenos) {
        print(top_pat_plot(pheno, scans, snp_w, FALSE,
                           snp_action = snp_action()))
        top_snp_asso(pheno, scans, snp_w)
      }
      dev.off()
    }
  )
  
  ## Return patterns
  reactive({
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
    fluidRow(
       column(6, downloadButton(ns("downloadData"), "CSV")),
       column(6, downloadButton(ns("downloadPlot"), "Plots"))))
}
#' @rdname shinySNPAssoc
#' @export
shinySNPAssocOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_scan"))
}
