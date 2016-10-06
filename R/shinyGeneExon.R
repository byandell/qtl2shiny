#' Shiny Genes and Exons with nearby SNPs module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param snp_par,chr_pos,top_snps_tbl,gene_exon_tbl reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyGeneExon <- function(input, output, session,
                     snp_par, chr_pos, top_snps_tbl, gene_exon_tbl) {
  ns <- session$ns
  
  pheno_names <- reactive({
    sort(unique(req(top_snps_tbl()$pheno)))
  })
  summary_gene_exon <- reactive({
    summary(req(gene_exon_tbl()),,
            req(top_snps_tbl()))
  })
  gene_names <- reactive({
    pheno_name <- req(snp_par$pheno_name)
    gene_in <- summary_gene_exon()
    as.data.frame(gene_in)[gene_in[[pheno_name]] > 0, "gene"]
  })
  gene_exon_pheno <- reactive({
    pheno_name <- req(snp_par$pheno_name)
    gene_in <- gene_names()
    if(length(gene_in)) {
      subset(gene_exon_tbl(), gene_in)
    } else {
      NULL
    }
  })
  
  output$gene_sum <- renderDataTable({
    withProgress(message = 'Gene Exon Table ...', value = 0, {
      setProgress(1)
      summary_gene_exon()
    })
  }, options = list(scrollX = TRUE, pageLength = 5,
                    lengthMenu = list(c(5,10,20,-1),
                                      list("5","10","20","all"))))
  output$gene_name <- renderUI({
    selected <- input$gene_name
    choices <- gene_names()
    if(!isTruthy(selected %in% choices))
      selected <- choices[1]
    selectInput(ns("gene_name"), NULL,
                choices = choices,
                selected = selected)
  })
  observeEvent(gene_names(), {
    selected <- input$gene_name
    choices <- gene_names()
    if(!isTruthy(selected %in% choices))
      selected <- choices[1]
    updateSelectInput(session, "gene_name", NULL,
                choices = choices,
                selected = selected)
  })

  output$gene_plot <- renderPlot({
    req(top_snps_tbl(), gene_exon_tbl(), gene_names())
    gene_name <- req(input$gene_name)
    pheno_name <- req(snp_par$pheno_name)
    withProgress(message = 'Gene Exon Plot ...', value = 0, {
      setProgress(1)
      plot_gene_exon(gene_exon_pheno(), 
                     top_snps_tbl() %>%
                       filter(pheno == pheno_name),
                     gene_name, pheno_name)
    })
  })
  
  ## Outputs
  output$exon_input <- renderUI({
    switch(req(input$button),
           Plot    = uiOutput(ns("gene_name")))
  })
  output$exon_output <- renderUI({
    switch(req(input$button),
           Plot    = plotOutput(ns("gene_plot")),
           Summary = dataTableOutput(ns("gene_sum")))
  })
  
  ## Downloads.
  output$downloadData <- downloadHandler(
    filename = function() {
      file.path(paste0("gene_exon_", chr_pos(), ".csv")) },
    content = function(file) {
      write.csv(req(gene_exon_tbl()), file)
    }
  )
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("gene_exon_", chr_pos(), ".pdf")) },
    content = function(file) {
      gene_exon <- req(gene_exon_pheno())
      pheno_name <- req(snp_par$pheno_name)
      top_snps <- req(top_snps_tbl()) %>%
        filter(pheno == pheno_name)
      pdf(file)
      for(gene_name in req(gene_names())) {
        print(plot_gene_exon(gene_exon, top_snps,
                             gene_name, pheno_name))
      }
      dev.off()
    }
  )
  output$select <- renderUI({
    selectInput(ns("button"), NULL, c("Plot","Summary"),
                input$button)
  })
}
#' @param id identifier for \code{\link{shinyScan1SNP}} use
#' @rdname shinyGeneExon
#' @export
shinyGeneExonInput <- function(id) {
  ns <- NS(id)
  fluidRow(
    uiOutput(ns("select")),
    uiOutput(ns("exon_input")))
}
#' @rdname shinyGeneExon
#' @export
shinyGeneExonUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(6, downloadButton(ns("downloadData"), "CSV")),
    column(6, downloadButton(ns("downloadPlot"), "Plots")))
}
#' @rdname shinyGeneExon
#' @export
shinyGeneExonOutput <- function(id) {
  ns <- NS(id)
  uiOutput(ns("exon_output"))
}
