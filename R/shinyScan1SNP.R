#' Shiny scan1 SNP analysis and plot module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param win_par,phe_df,cov_mx,pheno_anal,probs_obj,K_chr,dom_type reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyScan1SNP <- function(input, output, session,
                          win_par, phe_df, cov_mx,
                          pheno_anal, probs_obj, K_chr,
                          snp_action = reactive({NULL})) {
  ns <- session$ns

  chr_id <- reactive({win_par$chr_id})

  ## SNP analyses.
  snpprobs_obj <- reactive({
    withProgress(message = 'SNP Probs ...', value = 0, {
      setProgress(1)
      get_snpprobs(chr_id(), win_par$peak_Mbp, win_par$window_Mbp,
                   names(phe_df()), probs_obj(),
                   pattern = "AB1NZCPW", datapath)
    })
  })

  snpprobs_act <- reactive({
    snpprobs <- req(snpprobs_obj())
    snpprob_collapse(snpprobs, snp_action())
  })

  ## Scan1
  scan_obj <- reactive({
    req(pheno_anal())
    withProgress(message = "SNP Scan ...", value = 0, {
      setProgress(1)
      scan1(snpprobs_act(), phe_df(), K_chr(), cov_mx())
    })
  })

  # Scan Window slider
  output$scan_window <- renderUI({
    req(chr_id(),pheno_anal(),snpprobs_obj())
    rng <- round(range(snpprobs_obj()$map[[chr_id()]]), 2)
    sliderInput(ns("scan_window"), NULL, rng[1], rng[2],
                rng, step=.1)
  })

  ## Select phenotype for plots.
  output$pheno_assoc <- renderUI({
    req(pheno_anal())
    selectInput(ns("pheno_assoc"), NULL,
                choices = names(pheno_anal()),
                selected = input$pheno_assoc)
  })
  pheno_id <- reactive({
    pheno_anal <- req(input$pheno_assoc)
    pheno_anal()[pheno_anal]
  })

  ## Reactives for SNP analysis.
  ## Want to have slider for snp_w
  phename <- reactive({dimnames(scan_obj()$lod)[[2]]})
  output$snpPlot <- renderPlot({
    withProgress(message = 'SNP plots ...', value = 0, {
      setProgress(1)
      top_snp_asso(pheno_id(), scan_obj(), input$scan_window)
    })
  })
  output$snpPatternSum <- renderDataTable({
    req(pheno_id(), scan_obj())
    scan_snp <- scan_obj()
    if(max(scan_snp$lod) <= 1.5)
      return(NULL)
    summary(topsnp_pattern(scan_snp, pheno_id()))
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  output$snpPatternPlot <- renderPlot({
    req(pheno_id(), scan_obj())
    withProgress(message = 'SNP pattern plots ...', value = 0, {
      setProgress(1)
      top_pat_plot(pheno_id(), scan_obj(), input$scan_window, TRUE)
    })
  })
  output$snp_phe_pat <- renderPlot({
    withProgress(message = 'SNP Pheno patterns ...', value = 0, {
      setProgress(1)
      top_pat_plot(phename(), scan_obj(), input$scan_window,
                   group = "pheno", top_type = "SNP allele pattern")
    })
  })
  output$snp_pat_phe <- renderPlot({
    withProgress(message = 'SNP Pattern phenos ...', value = 0, {
      setProgress(1)
      top_pat_plot(phename(), scan_obj(), input$scan_window,
                   group = "pattern", top_type = "Phenotype")
    })
  })
  snp_coef <- reactive({
    req(pheno_anal())
    phename <- dimnames(scan_obj()$lod)[[2]]
    coef_topsnp(scan_obj(),probs_obj(),phe_df(),K_chr(),cov_mx())
  })
  output$scanTable <- renderDataTable({
    snp_coef()
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))

  ## Downloads.
  chr_pos <- reactive({
    make_chr_pos(chr_id(), range = input$scan_window)
  })

  output$snp_scan <- renderUI({
    switch(req(input$snp_scan),
           Association = {
             tagList(
               uiOutput(ns("pheno_assoc")),
               plotOutput(ns("snpPlot")))
             },
           Pattern = {
             tagList(
               uiOutput(ns("pheno_assoc")),
               plotOutput(ns("snpPatternPlot")))
             },
           "All Phenos" = {
             plotOutput(ns("snp_phe_pat"))
             },
           "All Patterns" = {
            plotOutput(ns("snp_pat_phe"))
            },
          Summary = {
             dataTableOutput(ns("snpPatternSum"))
            })
  })

  ## Downloads
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
      snp_w <- req(input$scan_window)
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
    fluidRow(
      column(3,
             h4(strong("SNP Plots")),
             radioButtons(ns("snp_scan"), "",
                          c("Association","Pattern",
                            "All Phenos","All Patterns",
                            "Summary")),
             uiOutput(ns("scan_window")),
             fluidRow(
               column(6, downloadButton(ns("downloadData"), "CSV")),
               column(6, downloadButton(ns("downloadPlot"), "Plots")))),
      column(9,
             uiOutput(ns("snp_scan"))))
  )
}
