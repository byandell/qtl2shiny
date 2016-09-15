#' Shiny top SNP analysis and plot module
#'
#' Shiny module for top SNP analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param chr_pos,snp_scan_obj,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#' 
#' @export
shinyTopSNP <- function(input, output, session,
                          chr_pos, snp_scan_obj,
                          snp_action = reactive({"basic"})) {
  ns <- session$ns

  phename <- reactive({dimnames(snp_scan_obj()$lod)[[2]]})

  
  sum_top_pat <- reactive({
    req(pheno_id())
    scan_snp <- req(snp_scan_obj())
    if(max(scan_snp$lod) <= 1.5)
      return(NULL)
    summary(topsnp_pattern(scan_snp, pheno_id()))
  })
  
  output$snpPatternSum <- renderDataTable({
    sum_top_pat()
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  
  output$snpPatternPlot <- renderPlot({
    if(is.null(pheno_id()) | is.null(snp_scan_obj()) |
       is.null(input$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP pattern plots ...', value = 0, {
      setProgress(1)
      top_pat_plot(pheno_id(), snp_scan_obj(), input$scan_window, TRUE,
                   snp_action = snp_action())
    })
  })
  
  ## SNP Pheno patterns
  output$snp_phe_pat <- renderPlot({
    if(is.null(phename()) | is.null(snp_scan_obj()) |
       is.null(input$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP Pheno patterns ...', value = 0, {
      setProgress(1)
      top_pat_plot(phename(), snp_scan_obj(), input$scan_window,
                   group = "pheno", snp_action = snp_action())
    })
  })
  
  ## SNP Pattern phenos
  output$snp_pat_phe <- renderPlot({
    if(is.null(phename()) | is.null(snp_scan_obj()) |
       is.null(input$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP Pattern phenos ...', value = 0, {
      setProgress(1)
      top_pat_plot(phename(), snp_scan_obj(), input$scan_window,
                   group = "pattern", snp_action = snp_action())
    })
  })

  ## SNP output choice
  output$snp_choice <- renderUI({
    if(req(input$snp_scan) %in% c("Top SNPs","Pattern")) {
      uiOutput(ns("pheno_assoc"))
    }
    if(req(input$snp_scan) == "Top SNPs")
      shinyTopFeatureUI(ns("top_feature"))
  })
  output$snp_scan <- renderUI({
    switch(req(input$snp_scan),
           "Top SNPs" = {
             shinyTopFeatureOutput(ns("top_feature"))
           },
           Pattern = {
             plotOutput(ns("snpPatternPlot"))
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
  output$title <- renderUI({
    if(snp_action() == "basic")
      h4(strong("SNP Plots"))
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
      scans <- req(snp_scan_obj())
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
}
#' @param id identifier for \code{\link{shinyTopSNP}} use
#' @rdname shinyTopSNP
#' @export
shinyTopSNPUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("title")),
    radioButtons(ns("snp_scan"), "",
                 c("Top SNPs","Pattern",
                    "All Phenos","All Patterns",
                   "Summary")),
    uiOutput(ns("snp_choice")),
    uiOutput(ns("scan_window")),
    fluidRow(
       column(6, downloadButton(ns("downloadData"), "CSV")),
       column(6, downloadButton(ns("downloadPlot"), "Plots"))))
}
#' @rdname shinyTopSNP
#' @export
shinyTopSNPOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("snp_scan"))
}
