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
#' @importFrom dplyr distinct
#' @importFrom doqtl2 topsnp_pattern
#' @importFrom qtl2ggplot sdp_to_pattern
#' @importFrom shiny callModule NS reactive req 
#'   radioButtons selectInput updateRadioButtons
#'   dataTableOutput plotOutput uiOutput
#'   renderDataTable renderPlot renderUI
#'   mainPanel sidebarPanel column strong tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
shinyAllelePat <- function(input, output, session,
                        snp_par, chr_pos, pheno_names,
                        snp_scan_obj, top_snps_tbl, 
                        gene_exon_tbl, 
                        snp_action = shiny::reactive({"basic"})) {
  ns <- session$ns

  ## Shiny Module
  shiny::callModule(shinyTopFeature, "top_feature",
             snp_par, chr_pos, 
             snp_scan_obj, top_snps_tbl, 
             gene_exon_tbl, snp_action)
  
  sum_top_pat <- shiny::reactive({
    summary(doqtl2::topsnp_pattern(shiny::req(snp_scan_obj()), 
                           pheno_names()))
  })
  
  output$snpPatternSum <- shiny::renderDataTable({
    sum_top_pat()
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10))
  
  output$snpPatternPlot <- shiny::renderPlot({
    if(is.null(snp_par$pheno_name) | is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    shiny::withProgress(message = 'SNP pattern plots ...', value = 0, {
      shiny::setProgress(1)
      top_pat_plot(snp_par$pheno_name, snp_scan_obj(), 
                   snp_par$scan_window, TRUE,
                   snp_action = snp_action())
    })
  })
  
  ## SNP Pheno patterns
  output$snp_phe_pat <- shiny::renderPlot({
    if(is.null(pheno_names()) | is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    shiny::withProgress(message = 'SNP Pheno patterns ...', value = 0, {
      shiny::setProgress(1)
      top_pat_plot(pheno_names(), snp_scan_obj(), snp_par$scan_window,
                   group = "pheno", snp_action = snp_action())
    })
  })

  top_pattern <- shiny::reactive({
    doqtl2::topsnp_pattern(snp_scan_obj(), pheno_names())
  })
  output$pattern <- shiny::renderUI({
    shiny::req(snp_action())
    top_pat <- shiny::req(top_pattern())
    choices <- qtl2ggplot::sdp_to_pattern(dplyr::distinct(top_pat, sdp)$sdp)
    if(!is.null(selected <- input$pattern)) {
      if(!selected %in% choices)
        selected <- NULL
    }
    shiny::selectInput(ns("pattern"), NULL,
                choices = choices,
                selected = selected)
  })
  ## SNP Pattern phenos
  output$snp_pat_phe <- shiny::renderPlot({
    if(is.null(pheno_names()) | is.null(snp_scan_obj()) |
       is.null(snp_par$scan_window) | is.null(snp_action()))
      return(plot_null())
    shiny::req(input$pattern)
    top_pat <- shiny::req(top_pattern())
    patterns <- qtl2ggplot::sdp_to_pattern(top_pat$sdp)
    shiny::withProgress(message = 'SNP Pattern phenos ...', value = 0, {
      shiny::setProgress(1)
      top_pat_plot(pheno_names(), snp_scan_obj(), snp_par$scan_window,
                   top_pattern = top_pat[patterns == input$pattern,],
                   group = "pattern", snp_action = snp_action())
    })
  })

  output$pat_input <- shiny::renderUI({
    switch(shiny::req(input$button),
           "All Patterns" = shiny::uiOutput(ns("pattern")),
           "Top SNPs"     = shinyTopFeatureInput(ns("top_feature")))
  })
  output$pat_output <- shiny::renderUI({
    switch(shiny::req(input$button),
           "Top SNPs"     = shinyTopFeatureOutput(ns("top_feature")),
           Pattern        = shiny::plotOutput(ns("snpPatternPlot")),
           "All Phenos"   = shiny::plotOutput(ns("snp_phe_pat")),
           "All Patterns" = shiny::plotOutput(ns("snp_pat_phe")),
           Summary        = shiny::dataTableOutput(ns("snpPatternSum")))
  })
  output$title <- shiny::renderUI({
    if(snp_action() == "basic")
      shiny::strong("SNP Plots")
  })

  ## Downloads
  output$download_csv_plot <- shiny::renderUI({
    switch(shiny::req(input$button),
           "Top SNPs"     = shinyTopFeatureUI(ns("top_feature")),
           shiny::tagList(fluidRow(
             shiny::column(6, shiny::downloadButton(ns("downloadData"), "CSV")),
             shiny::column(6, shiny::downloadButton(ns("downloadPlot"), "Plots")))))
  })
  output$downloadData <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("pattern_", chr_pos(), "_", snp_action(), ".csv")) },
    content = function(file) {
      write.csv(sum_top_pat(), file)
    }
  )
  output$downloadPlot <- shiny::downloadHandler(
    filename = function() {
      file.path(paste0("pattern_", chr_pos(), "_", snp_action(), ".pdf")) },
    content = function(file) {
      scans <- shiny::req(snp_scan_obj())
      snp_w <- shiny::req(snp_par$scan_window)
      phenos <- shiny::req(pheno_names())
      pdf(file, width = 9)
      ## Plots over all phenotypes
      print(top_pat_plot(phenos, scans, snp_w,
                         top_pattern = top_pattern(),
                         group = "pheno", snp_action = snp_action()))

      top_pat <- shiny::req(top_pattern())
      patterns <- qtl2ggplot::sdp_to_pattern(top_pat$sdp)
      upat <- unique(patterns)
      for(pattern in upat) {
        print(top_pat_plot(pheno_names(), snp_scan_obj(), 
                           snp_par$scan_window,
                           top_pattern = 
                             top_pattern()[patterns == pattern,],
                           group = "pattern", 
                           snp_action = snp_action()))
      }
      
      ## Plots by phenotype.
      for(pheno in phenos) {
        print(top_pat_plot(pheno, scans, snp_w, FALSE,
                           snp_action = snp_action()))
        top_snp_asso(pheno, scans, snp_w)
      }
      dev.off()
    }
  )
  output$radio <- shiny::renderUI({
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
    shiny::radioButtons(ns("button"), "",
                 button_val, selected)
  })
  ## Update Radio Button if 1 or >1 Phenotype Names.
  shiny::observeEvent(pheno_names(), {
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
      shiny::updateRadioButtons(session, "button", 
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
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("radio")),
    shiny::uiOutput(ns("pat_input"))
  )
}
#' @rdname shinyAllelePat
#' @export
shinyAllelePatUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("download_csv_plot"))
}
#' @rdname shinyAllelePat
#' @export
shinyAllelePatOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("pat_output"))
}
