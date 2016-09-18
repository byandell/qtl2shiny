#' Shiny SNP and Allele analysis and plot module
#'
#' Shiny module to coordinate SNP and allele analyses and plots.
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
shinySNPAllele <- function(input, output, session,
                          win_par, phe_df, cov_mx,
                          pheno_anal, probs_obj, K_chr,
                          ui_task, snp_action = reactive({"basic"})) {
  ns <- session$ns
  
  ##### Make sure snp_scan_obj gets pheno_anal right.

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
  ## SNP Scan.
  snp_scan_obj <- reactive({
    req(pheno_anal())
    snpprobs_act <- 
      snpprob_collapse(req(snpprobs_obj()), snp_action())
    withProgress(message = "SNP Scan ...", value = 0, {
      setProgress(1)
      scan1(snpprobs_act, phe_df(), K_chr(), cov_mx())
    })
  })
  ## Top SNPs table.
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
  ass_par <- callModule(shinySNPAssoc, "snp_assoc",
             input, chr_pos, snp_scan_obj, top_snps_tbl, 
             gene_exon_tbl, snp_action)
  ## Allele Patterns
  pat_par <- callModule(shinyAllelePat, "allele_pat",
             input, chr_pos, pheno_anal, snp_scan_obj, top_snps_tbl, 
             gene_exon_tbl, snp_action)

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
        ass_par$button %in% c("Scan")) |
       (ui_task() == "Allele Pattern" &
        input$allele_pat %in% c("Pattern","Top SNPs"))) {
      uiOutput(ns("pheno_assoc"))
    }
  })
  output$win_choice <- renderUI({
    ## Show Window for Scan, Genes, Pattern, Alls
    if((ui_task() == "SNP Association" & 
        ass_par$button %in% c("Scan","Genes")) |
       (ui_task() == "Allele Pattern" &
        input$allele_pat %in% c("Pattern","All Phenos","All Patterns"))) {
      uiOutput(ns("scan_window"))
    }
  })

  ## UI Logic
  output$title <- renderUI({
    if(snp_action() == "basic")
      h4(strong(req(ui_task())))
  })
  output$snp_choice <- renderUI({
    switch(req(ui_task()),
           "SNP Association" = shinySNPAssocInput(ns("snp_assoc")),
           "Allele Pattern"  = shinyAllelePatInput(ns("allele_pat")))
  })
  output$snp_output <- renderUI({
    switch(req(ui_task()),
           "SNP Association" = shinySNPAssocOutput(ns("snp_assoc")),
           "Allele Pattern"  = shinyAllelePatOutput(ns("allele_pat")))
  })

  ## Downloads
  output$download_csv_plot <- renderUI({
    switch(req(ui_task()),
           "SNP Association" = shinySNPAssocUI(ns("snp_assoc")),
           "Allele Pattern"  = shinyAllelePatUI(ns("allele_pat")))
  })
  
  ## Return patterns
  reactive({ # patterns
    summary(top_snps_tbl()) %>%
      filter(max_lod >= 3) %>%
      mutate(contrast = snp_action()) %>%
      arrange(desc(max_lod))
  })
}
#' @param id identifier for \code{\link{shinySNPAllele}} use
#' @rdname shinySNPAllele
#' @export
shinySNPAlleleUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("title")),
    uiOutput(ns("phe_choice")),
    uiOutput(ns("win_choice")),
    uiOutput(ns("snp_choice")),
    uiOutput(ns("download_csv_plot")))
}
#' @rdname shinySNPAllele
#' @export
shinySNPAlleleOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_output"))
}
