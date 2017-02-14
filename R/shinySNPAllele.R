#' Shiny SNP and Allele analysis and plot module
#'
#' Shiny module to coordinate SNP and allele analyses and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,probs_obj,K_chr,data_path,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom dplyr arrange desc filter mutate
#' @importFrom CCSanger get_gene_exon_snp 
#' @importFrom qtl2pattern top_snps_all snpprob_collapse
#' @importFrom qtl2scan scan1
#' @importFrom shiny callModule NS reactive req 
#'   selectInput sliderInput
#'   uiOutput
#'   renderUI
#'   tagList
#'   withProgress setProgress
shinySNPAllele <- function(input, output, session,
                          job_par, win_par, 
                          phe_df, cov_mx,
                          probs_obj, K_chr,
                          data_path,
                          snp_action = shiny::reactive({"basic"})) {
  ns <- session$ns
  
  chr_pos <- shiny::reactive({
    make_chr_pos(shiny::req(win_par$chr_id), 
                 range = shiny::req(input$scan_window))
  })
  pheno_names <- shiny::reactive({
    names(phe_df())
  })
  
  ## Reactives
  ## SNP Probabilities.
  snpprobs_obj <- shiny::callModule(shinySNPProbs, "snp_probs",
                  win_par, pheno_names, probs_obj,
                  data_path)
  ## SNP Scan.
  snp_scan_obj <- shiny::reactive({
    shiny::req(phe_df())
    snpprobs_act <- 
      qtl2pattern::snpprob_collapse(shiny::req(snpprobs_obj()), snp_action())
    shiny::withProgress(message = "SNP Scan ...", value = 0, {
      shiny::setProgress(1)
      qtl2scan::scan1(snpprobs_act, phe_df(), K_chr(), cov_mx())
    })
  })
  ## Top SNPs table.
  top_snps_tbl <- shiny::reactive({
    shiny::req(snp_action())
    shiny::withProgress(message = 'Get Top SNPs ...', value = 0, {
      shiny::setProgress(1)
      qtl2pattern::top_snps_all(shiny::req(snp_scan_obj()))
    })
  })
  ## Genes and Exons.
  gene_exon_tbl <- shiny::reactive({
    shiny::req(snp_action())
    shiny::withProgress(message = 'Gene Exon Calc ...', value = 0, {
      shiny::setProgress(1)
      tops <- shiny::req(top_snps_tbl())
      CCSanger::get_gene_exon_snp(tops,
                        file.path(data_path(), "mgi_db.sqlite"))
    })
  })
  
  ## Shiny Modules
  ## SNP Association
  ass_par <- shiny::callModule(shinySNPAssoc, "snp_assoc",
             input, chr_pos, pheno_names,
             snp_scan_obj, top_snps_tbl, 
             gene_exon_tbl, data_path,
             snp_action)
  ## Allele Patterns
  pat_par <- shiny::callModule(shinyAllelePat, "allele_pat",
             input, chr_pos, pheno_names,
             snp_scan_obj, top_snps_tbl, 
             gene_exon_tbl, snp_action)

  # Scan Window slider
  output$scan_window <- shiny::renderUI({
    rng <- round(shiny::req(win_par$peak_Mbp) + 
                   c(-1,1) * shiny::req(win_par$window_Mbp), 
                 1)
    if(is.null(selected <- input$scan_window))
      selected <- rng
    shiny::sliderInput(ns("scan_window"), NULL, rng[1], rng[2],
                selected, step=.1)
  })

  ## Select phenotype for plots.
  output$pheno_name <- shiny::renderUI({
    shiny::selectInput(ns("pheno_name"), NULL,
                choices = shiny::req(pheno_names()),
                selected = input$pheno_name)
  })
  
  ## Button Options.
  output$phe_choice <- shiny::renderUI({
    ## Show Phenotype Choice for Scan, Pattern, Top SNPs
    pheno_choice <- FALSE
    switch(shiny::req(job_par$button),
           "SNP Association" = {
             if(!is.null(ass_par$button)) {
               pheno_choice <- (ass_par$button %in% c("Genes","Exons"))
             }
           },
           "Allele Pattern" = {
             if(!is.null(pat_par$button)) {
               pheno_choice <- (pat_par$button %in% c("Pattern","Top SNPs"))
             }
           })
    if(pheno_choice) {
      shiny::uiOutput(ns("pheno_name"))
    }
  })
  output$win_choice <- shiny::renderUI({
    ## Show Window for Scan, Genes, Pattern, Alls
    win_choice <- FALSE
    switch(shiny::req(job_par$button),
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
      shiny::uiOutput(ns("scan_window"))
    }
  })

  ## UI Logic
  output$title <- shiny::renderUI({
    if(snp_action() == "basic")
      strong(shiny::req(job_par$button))
  })
  output$snp_input <- shiny::renderUI({
    switch(shiny::req(job_par$button),
           "SNP Association" = shinySNPAssocInput(ns("snp_assoc")),
           "Allele Pattern"  = shinyAllelePatInput(ns("allele_pat")))
  })
  output$snp_output <- shiny::renderUI({
    switch(shiny::req(job_par$button),
           "SNP Association" = shinySNPAssocOutput(ns("snp_assoc")),
           "Allele Pattern"  = shinyAllelePatOutput(ns("allele_pat")))
  })

  ## Downloads
  output$download_csv_plot <- shiny::renderUI({
    switch(shiny::req(job_par$button),
           "SNP Association" = shinySNPAssocUI(ns("snp_assoc")),
           "Allele Pattern"  = shinyAllelePatUI(ns("allele_pat")))
  })
  
  ## Return patterns
  shiny::reactive({ # patterns
    dplyr::arrange(
      dplyr::mutate(
        dplyr::filter(
          summary(top_snps_tbl()), 
          max_lod >= 3), 
        contrast = snp_action()), 
      dplyr::desc(max_lod))
  })
}
#' @param id identifier for \code{\link{shinySNPAllele}} use
#' @rdname shinySNPAllele
#' @export
shinySNPAlleleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("title")),
    shiny::uiOutput(ns("snp_input")),
    shiny::uiOutput(ns("phe_choice")),
    shiny::uiOutput(ns("win_choice")),
    shiny::uiOutput(ns("download_csv_plot")))
}
#' @rdname shinySNPAllele
#' @export
shinySNPAlleleOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("snp_output"))
}
