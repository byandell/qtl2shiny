#' Shiny SNP and Allele analysis and plot module
#'
#' Shiny module to coordinate SNP and allele analyses and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,probs_obj,K_chr,snp_action reactive arguments
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
                          job_par, win_par, 
                          phe_df, cov_mx,
                          probs_obj, K_chr,
                          snp_action = reactive({"basic"})) {
  ns <- session$ns
  
  chr_pos <- reactive({
    make_chr_pos(req(win_par$chr_id), range = req(input$scan_window))
  })
  pheno_names <- reactive({
    names(phe_df())
  })
  
  ## Reactives
  ## SNP Probabilities.
  snpprobs_obj <- callModule(shinySNPProbs, "snp_probs",
                  win_par, pheno_names, probs_obj)
  ## SNP Scan.
  snp_scan_obj <- reactive({
    cat(file = stderr(), "snp_scan_obj", snp_action(), "\n")
    req(phe_df())
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
             input, chr_pos, snp_scan_obj, top_snps_tbl, 
             gene_exon_tbl, snp_action)

  # Scan Window slider
  output$scan_window <- renderUI({
    rng <- round(req(win_par$peak_Mbp) + 
                   c(-1,1) * req(win_par$window_Mbp), 
                 1)
    if(is.null(selected <- input$scan_window))
      selected <- rng
    sliderInput(ns("scan_window"), NULL, rng[1], rng[2],
                selected, step=.1)
  })

  ## Select phenotype for plots.
  output$pheno_assoc <- renderUI({
    selectInput(ns("pheno_assoc"), NULL,
                choices = req(pheno_names()),
                selected = input$pheno_assoc)
  })
  
  ## Button Options.
  output$phe_choice <- renderUI({
    ## Show Phenotype Choice for Scan, Pattern, Top SNPs
    pheno_choice <- FALSE
    switch(req(job_par$button),
           "SNP Association" = {
             if(!is.null(ass_par$button)) {
               pheno_choice <- (ass_par$button %in% c("Scan"))
             }
           },
           "Allele Pattern" = {
             if(!is.null(pat_par$button)) {
               pheno_choice <- (pat_par$button %in% c("Pattern","Top SNPs"))
             }
           })
    if(pheno_choice) {
      uiOutput(ns("pheno_assoc"))
    }
  })
  output$win_choice <- renderUI({
    ## Show Window for Scan, Genes, Pattern, Alls
    win_choice <- FALSE
    switch(req(job_par$button),
           "SNP Association" = {
             if(!is.null(ass_par$button)) {
               win_choice <- (ass_par$button %in% c("Scan","Genes"))
             }
           },
           "Allele Pattern" = {
             if(!is.null(pat_par$button)) {
               win_choice <- 
                 (pat_par$button %in%
                    c("Pattern","All Phenos","All Patterns"))
             }
           })
    if(win_choice) {
      uiOutput(ns("scan_window"))
    }
  })

  ## UI Logic
  output$title <- renderUI({
    if(snp_action() == "basic")
      h4(strong(req(job_par$button)))
  })
  output$snp_input <- renderUI({
    switch(req(job_par$button),
           "SNP Association" = shinySNPAssocInput(ns("snp_assoc")),
           "Allele Pattern"  = shinyAllelePatInput(ns("allele_pat")))
  })
  output$snp_output <- renderUI({
    switch(req(job_par$button),
           "SNP Association" = shinySNPAssocOutput(ns("snp_assoc")),
           "Allele Pattern"  = shinyAllelePatOutput(ns("allele_pat")))
  })

  ## Downloads
  output$download_csv_plot <- renderUI({
    switch(req(job_par$button),
           "SNP Association" = shinySNPAssocUI(ns("snp_assoc")),
           "Allele Pattern"  = shinyAllelePatUI(ns("allele_pat")))
  })
  
  ## Return patterns
  reactive({ # patterns
    cat(file = stderr(), "patterns\n")
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
    uiOutput(ns("snp_input")),
    uiOutput(ns("phe_choice")),
    uiOutput(ns("win_choice")),
    uiOutput(ns("download_csv_plot")))
}
#' @rdname shinySNPAllele
#' @export
shinySNPAlleleOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_output"))
}
