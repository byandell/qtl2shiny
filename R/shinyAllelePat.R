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
shinyAllelePat <- function(input, output, session,
                        snp_par, chr_pos, pheno_anal, snp_scan_obj, 
                        top_snps_tbl, gene_exon_tbl, 
                        snp_action = reactive({"basic"})) {
  ns <- session$ns

  ## Phenotype names.
  phename <- reactive({
    dimnames(req(snp_scan_obj())$lod)[[2]]
  })
  pheno_id <- reactive({
    pheno_anal <- req(phename)
    pheno_anal()[pheno_anal]
  })
  observeEvent(phename(), {
    button_val <- c("Top SNPs","Pattern",
                    "All Phenos","All Patterns",
                    "Summary")
    if(length(phename() == 1)) {
      button_val <- button_val[-(3:4)]
    }
    updateRadioButtons(session, "snp_assoc", 
                       selected = button_val[1],
                       choices = button_val)
  })
  

  ## Shiny Module
  callModule(shinyTopFeature, "top_feature",
             chr_pos, snp_scan_obj, top_snps_tbl, gene_exon_tbl)
  
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
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP pattern plots ...', value = 0, {
      setProgress(1)
      top_pat_plot(pheno_id(), snp_scan_obj(), snp_par$scan_window, TRUE,
                   snp_action = snp_action())
    })
  })
  
  ## SNP Pheno patterns
  output$snp_phe_pat <- renderPlot({
    if(is.null(phename()) | is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP Pheno patterns ...', value = 0, {
      setProgress(1)
      top_pat_plot(phename(), snp_scan_obj(), snp_par$scan_window,
                   group = "pheno", snp_action = snp_action())
    })
  })
  
  ## SNP Pattern phenos
  output$snp_pat_phe <- renderPlot({
    if(is.null(phename()) | is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP Pattern phenos ...', value = 0, {
      setProgress(1)
      top_pat_plot(phename(), snp_scan_obj(), snp_par$scan_window,
                   group = "pattern", snp_action = snp_action())
    })
  })

  output$pat_input <- renderUI({
    switch(req(input$button),
           "Top SNPs"     = shinyTopFeatureUI(ns("top_feature")))
  })
  output$pat_output <- renderUI({
    switch(req(input$button),
           "Top SNPs"     = shinyTopFeatureOutput(ns("top_feature")),
           Pattern        = plotOutput(ns("snpPatternPlot")),
           "All Phenos"   = plotOutput(ns("snp_phe_pat")),
           "All Patterns" = plotOutput(ns("snp_pat_phe")),
           Summary        = dataTableOutput(ns("snpPatternSum")))
  })
  output$title <- renderUI({
    if(snp_action() == "basic")
      h4(strong("SNP Plots"))
  })

  ## Downloads
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("pattern_", chr_pos(), ".csv")) },
    content = function(file) {
      write.csv(sum_top_pat(), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("pattern_", chr_pos(), ".pdf")) },
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
  input
}
#' @param id identifier for \code{\link{shinyAllelePat}} use
#' @rdname shinyAllelePat
#' @export
shinyAllelePatInput <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("button"), "",
                 c("Top SNPs","Pattern",
                   "All Phenos","All Patterns",
                   "Summary")),
    uiOutput(ns("pat_input"))
  )
}
#' @rdname shinyAllelePat
#' @export
shinyAllelePatUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(6, downloadButton(ns("downloadData"), "CSV")),
    column(6, downloadButton(ns("downloadPlot"), "Plots")))
}
#' @rdname shinyAllelePat
#' @export
shinyAllelePatOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("pat_output"))
}
