#' Shiny scan1 SNP analysis and plot module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param chr_id,phe_df,cov_mx,pheno_anal,probs_obj,K_chr reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyScan1SNP <- function(input, output, session,
                          chr_id, phe_df, cov_mx, 
                          pheno_anal, probs_obj, K_chr) {
  ns <- session$ns

  ## SNP analyses.
  snpprobs_obj <- reactive({
    window_Mbp <- req(win_par$window_Mbp)
    peak_Mbp <- req(win_par$peak_Mbp)
    withProgress(message = 'SNP Scan ...', value = 0, {
      setProgress(1)
      get_snpprobs(chr_id(), peak_Mbp, window_Mbp,
                   names(phe_df()), probs_obj(),
                   pattern = "AB1NZCPW", datapath)
    })
  })
  
  # Scan Window slider
  output$choose_scan_window <- renderUI({
    req(chr_id(),pheno_anal(),snpprobs_obj())
    rng <- round(range(snpprobs_obj()$map[[chr_id()]]), 2)
    sliderInput(ns("scan_window"), NULL, rng[1], rng[2],
                rng, step=.1)
  })
  
  ## Select phenotype for plots.
  output$choose_pheno <- renderUI({
    req(pheno_anal())
    selectInput(ns("pheno_anal"), NULL,
                choices = names(pheno_anal()))
  })
  pheno_id <- reactive({
    pheno_anal <- req(input$pheno_anal)
    pheno_anal()[pheno_anal]
  })

  ## Reactives for SNP analysis.
  ## Want to have slider for snp_w
  phename <- reactive({dimnames(scan_obj()$lod)[[2]]})
  output$snpPlot <- renderPlot({
    withProgress(message = 'SNP plots ...', value = 0, {
      setProgress(1)
      top_snp_asso(pheno_id(), scan_obj(), scan_window())
    })
  })
  output$snpPatternSum <- renderTable({
    req(pheno_id(), scan_obj())
    scan_snp <- scan_obj()
    if(max(scan_snp$lod) <= 1.5)
      return(NULL)
    summary(topsnp_pattern(scan_snp, pheno_id()))
  })
  output$snpPatternPlot <- renderPlot({
    req(pheno_id(), scan_obj())
    withProgress(message = 'SNP pattern plots ...', value = 0, {
      setProgress(1)
      top_pat_plot(pheno_id(), scan_obj(), scan_window(), TRUE)
    })
  })
  output$snp_phe_pat <- renderPlot({
    withProgress(message = 'SNP Pheno patterns ...', value = 0, {
      setProgress(1)
      top_pat_plot(phename(), scan_obj(), scan_window(),
                   group = "pheno", top_type = "SNP allele pattern")
    })
  })
  output$snp_pat_phe <- renderPlot({
    withProgress(message = 'SNP Pattern phenos ...', value = 0, {
      setProgress(1)
      top_pat_plot(phename(), scan_obj(), scan_window(),
                   group = "pattern", top_type = "Phenotype")
    })
  })
  snp_coef <- reactive({
    req(plot_type(),pheno_anal())
    phename <- dimnames(scan_obj()$lod)[[2]]
    coef_topsnp(scan_obj(),probs_obj(),phe_df(),K_chr(),cov_mx())
  })
  output$scanTable <- renderDataTable({
    snp_coef()
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  
  ## Downloads.
  chr_pos <- reactive({
    make_chr_pos(chr_id(), range = scan_window())
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_effects_", chr_pos(), ".csv")) },
    content = function(file) {
      write.csv(snp_coef(), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_scan_", chr_pos(), ".pdf")) },
    content = function(file) {
      scans <- req(scan_obj())
      snp_w <- req(scan_window())
      phenos <- req(phename())
      pdf(file)
      ## Plots over all phenotypes
      print(top_pat_plot(phenos, scans, snp_w, group = "pheno",
                         top_type = "SNP allele pattern"))
      print(top_pat_plot(phenos, scans, snp_w, group = "pattern",
                         top_type = "Phenotype"))
      ## Plots by phenotype.
      for(pheno in phenos) {
        print(top_pat_plot(pheno, scans, snp_w, FALSE))
        top_snp_asso(pheno, scans, snp_w)
      }
      dev.off()
    }
  )
  scan_obj
}

#' UI for shinyScan1SNP Shiny Module
#'
#' UI for scan1 analyses and summaries to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyScan1SNP
#' @export
shinyScan1SNPOutput <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("snpPatternPlot")),
    plotOutput(ns("snpPlot"))
  )
}
#' Output for shinyScan1SNP Shiny Module
#'
#' Output for scan1 analyses and summaries to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyScan1SNP
#' @export
shinyScan1SNPUI <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      tabPanel("summary",
               downloadButton(ns("downloadPlot"),
                              "Download Plots"),
               tableOutput(ns("snpPatternSum"))),
      tabPanel("by pheno",
               plotOutput(ns("snpPatternPlot")),
               plotOutput(ns("snpPlot"))),
      tabPanel("all phenos",
               plotOutput(ns("snp_phe_pat"))),
      tabPanel("all patterns",
               plotOutput(ns("snp_pat_phe")))
    )
  )
}
