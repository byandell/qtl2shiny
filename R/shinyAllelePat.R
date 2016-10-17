#' Shiny top SNP analysis and plot module
#'
#' Shiny module for top SNP analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param snp_par,chr_pos,pheno_names,snp_scan_obj,top_snps_tbl,gene_exon_tbl,snp_action reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#' 
#' @export
shinyAllelePat <- function(input, output, session,
                        snp_par, chr_pos, pheno_names,
                        snp_scan_obj, top_snps_tbl, 
                        gene_exon_tbl, 
                        snp_action = reactive({"basic"})) {
  ns <- session$ns

  ## Shiny Module
  callModule(shinyTopFeature, "top_feature",
             snp_par, chr_pos, 
             snp_scan_obj, top_snps_tbl, 
             gene_exon_tbl, snp_action)
  
  sum_top_pat <- reactive({
    scan_snp <- req(snp_scan_obj())
    if(max(scan_snp$lod) <= 1.5)
      return(NULL)
#    summary(topsnp_pattern(scan_snp, pheno_names()))
    topsnp_pattern(scan_snp, pheno_names()) %>%
      mutate(pattern = sdp_to_pattern(sdp)) %>%
      arrange(desc(lod))
  })
  
  output$snpPatternSum <- renderDataTable({
    sum_top_pat()
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  
  output$snpPatternPlot <- renderPlot({
    if(is.null(snp_par$pheno_name) | is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP pattern plots ...', value = 0, {
      setProgress(1)
      top_pat_plot(snp_par$pheno_name, snp_scan_obj(), 
                   snp_par$scan_window, TRUE,
                   snp_action = snp_action())
    })
  })
  
  ## SNP Pheno patterns
  output$snp_phe_pat <- renderPlot({
    if(is.null(pheno_names()) | is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP Pheno patterns ...', value = 0, {
      setProgress(1)
      top_pat_plot(pheno_names(), snp_scan_obj(), snp_par$scan_window,
                   group = "pheno", snp_action = snp_action())
    })
  })
  ## SNP Pattern phenos
  output$snp_pat_phe <- renderPlot({
    if(is.null(pheno_names()) | is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    withProgress(message = 'SNP Pattern phenos ...', value = 0, {
      setProgress(1)
      top_pat_plot(pheno_names(), snp_scan_obj(), snp_par$scan_window,
                   group = "pattern", snp_action = snp_action())
    })
  })

  output$pat_input <- renderUI({
    switch(req(input$button),
           "Top SNPs"     = shinyTopFeatureInput(ns("top_feature")))
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
  output$download_csv_plot <- renderUI({
    switch(req(input$button),
           "Top SNPs"     = shinyTopFeatureUI(ns("top_feature")),
           tagList(fluidRow(
             column(6, downloadButton(ns("downloadData"), "CSV")),
             column(6, downloadButton(ns("downloadPlot"), "Plots")))))
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("pattern_", chr_pos(), "_", snp_action(), ".csv")) },
    content = function(file) {
      write.csv(sum_top_pat(), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("pattern_", chr_pos(), "_", snp_action(), ".pdf")) },
    content = function(file) {
      scans <- req(snp_scan_obj())
      snp_w <- req(snp_par$scan_window)
      phenos <- req(pheno_names())
      pdf(file, width = 9)
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
  output$radio <- renderUI({
    button_val <- c("Pattern",
                    "All Phenos","All Patterns",
                    "Top SNPs","Summary")
    if(length(pheno_names()) == 1) {
      button_val <- button_val[-(2:3)]
    }
    if(!is.null(selected <- input$button)) {
      if(!(selected %in% button_val))
        selected <- button_val[1]
    }
    radioButtons(ns("button"), "",
                 button_val, selected)
  })
  ## Update Radio Button if 1 or >1 Phenotype Names.
  observeEvent(pheno_names(), {
    button_val <- c("Pattern",
                    "All Phenos","All Patterns",
                    "Top SNPs","Summary")
    if(length(pheno_names()) == 1) {
      button_val <- button_val[-(2:3)]
    }
    selected <- input$button
    if(!is.null(selected)) {
      if(!(selected %in% button_val))
        selected <- button_val[1]
      updateRadioButtons(session, "button", 
                         selected = selected,
                         choices = button_val)
    }
  })

  input
}
#' @param id identifier for \code{\link{shinyAllelePat}} use
#' @rdname shinyAllelePat
#' @export
shinyAllelePatInput <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("radio")),
    uiOutput(ns("pat_input"))
  )
}
#' @rdname shinyAllelePat
#' @export
shinyAllelePatUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("download_csv_plot"))
}
#' @rdname shinyAllelePat
#' @export
shinyAllelePatOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("pat_output"))
}
