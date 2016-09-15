#' Shiny scan1 SNP analysis and plot module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,pheno_anal,probs_obj,K_chr,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyScan1SNP <- function(input, output, session,
                          win_par, snp_par, phe_df, cov_mx,
                          pheno_anal, probs_obj, K_chr,
                          snp_action = reactive({"basic"})) {
  ns <- session$ns

  chr_id <- reactive({win_par()$chr_id})

  ## SNP analyses.
  snpprobs_obj <- reactive({
    withProgress(message = 'SNP Probs ...', value = 0, {
      setProgress(1)
      get_snpprobs(chr_id(), win_par()$peak_Mbp, win_par()$window_Mbp,
                   names(phe_df()), probs_obj(),
                   datapath)
    })
  })

  snpprobs_act <- reactive({
    snpprobs <- req(snpprobs_obj())
    snpprob_collapse(snpprobs, snp_action())
  })

  ## Scan1
  snp_scan_obj <- reactive({
    req(pheno_anal())
    withProgress(message = "SNP Scan ...", value = 0, {
      setProgress(1)
      scan1(snpprobs_act(), phe_df(), K_chr(), cov_mx())
    })
  })
  
  #################################################
  ## THIS IS NOT USED RIGHT NOT. WHY NOT?
  ## Scan1 of coefficients/effects
  snp_coef <- reactive({
    req(pheno_anal())
    coef_topsnp(snp_scan_obj(),snpprobs_act(),phe_df(),K_chr(),cov_mx())
  })
  output$scanTable <- renderDataTable({
    snp_coef()
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  #################################################

  pheno_id <- reactive({
    pheno_anal <- req(snp_par$pheno_assoc)
    pheno_anal()[pheno_anal]
  })

  ## Reactives for SNP analysis.
  ## Want to have slider for snp_w
  phename <- reactive({dimnames(snp_scan_obj()$lod)[[2]]})
  output$snpPlot <- renderPlot({
    if(is.null(pheno_id()) | is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP plots ...', value = 0, {
      setProgress(1)
      top_snp_asso(pheno_id(), snp_scan_obj(), snp_par$scan_window, snp_action())
    })
  })
  

  ## Downloads.
  chr_pos <- reactive({
    make_chr_pos(chr_id(), range = snp_par$scan_window)
  })
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_scan_", chr_pos(), ".pdf")) },
    content = function(file) {
      scans <- req(snp_scan_obj())
      snp_w <- req(snp_par$scan_window)
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
  snp_scan_obj
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyScan1SNP
#' @export
shinyScan1SNPUI <- function(id) {
  ns <- NS(id)
  tagList(
    downloadButton(ns("downloadPlot"), "Plots"))
}
#' @rdname shinyScan1SNP
#' @export
shinyScan1SNPOutput <- function(id) {
  ns <- NS(id)
  plotOutput(ns("snpPlot")),
}
