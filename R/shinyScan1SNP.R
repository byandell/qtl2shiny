#' Shiny scan1 SNP analysis and plot module
#'
#' Shiny module for scan1 analysis and plots.
#'
#' @param input,output,session standard shiny arguments
#' @param chr_id,scan_window,chr_pos,pheno_id,scan_obj,probs_obj reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyScan1SNP <- function(input, output, session,
                          chr_id, scan_window, chr_pos, pheno_id,
                          scan_obj, probs_obj) {
  ## Reactives for SNP analysis.
  ## Want to have slider for snp_w
  phename <- reactive({dimnames(scan_obj()$lod)[[2]]})
  top_snp_asso <- function(pheno, phename, snp_window) {
    plot_snpasso(subset(scan_obj(),
                        lodcolumn=match(pheno, phename())),
                 show_all_snps=FALSE, drop.hilit=1.5,
                 xlim=snp_window)
    title(paste(pheno, "chr", chr_id()))
  }
  output$snpPlot <- renderPlot({
    withProgress(message = 'SNP plots ...', value = 0,
                 {
                   setProgress(1)
                   top_snp_asso(pheno_id(), phename(), scan_window())
                 })
  })
  output$snpPatternSum <- renderTable({
    req(pheno_id(), scan_obj())
    scan_snp <- scan_obj()
    if(max(scan_snp$lod) <= 1.5)
      return(NULL)
    top_pattern <- topsnp_pattern(scan_snp, pheno_id())
    summary(top_pattern)
  })
  top_pat_plot <- function(pheno, fill.null=TRUE,
                           group = "pheno",
                           top_type = "SNP allele pattern") {
    snp_w <- scan_window()
    scan_snp <- scan_obj()
    drop <- max(max(scan_snp$lod) - 1.5, 1.5)
    top_pattern <- topsnp_pattern(scan_snp, pheno, drop)
    if(is.null(top_pattern)) {
      if(fill.null)
        return(ggplot(,aes(x=1,y=1,label="no SNP allele patterns to plot")) +
                 geom_text(size=10) + theme_void())
      else
        return()
    }
    plot(top_pattern, group=group) + xlim(snp_w) +
      ggtitle(paste(pheno, "on chr", chr_id(), "by", top_type))
  }
  output$snpPatternPlot <- renderPlot({
    req(pheno_id(), scan_obj())
    withProgress(message = 'SNP pattern plots ...', value = 0,
                 {
                   setProgress(1)
                   top_pat_plot(pheno_id(), TRUE)
                 })
  })
  output$snp_phe_pat <- renderPlot({
    withProgress(message = 'SNP Pheno patterns ...', value = 0,
                 {
                   setProgress(1)
                   top_pat_plot(phename(), group = "pheno",
                                top_type = "SNP allele pattern")
                 })
  })
  output$snp_pat_phe <- renderPlot({
    withProgress(message = 'SNP Pattern phenos ...', value = 0,
                 {
                   setProgress(1)
                   top_pat_plot(phename(), group = "pattern",
                                top_type = "Phenotype")
                 })
  })
  output$downloadPlot <- downloadHandler(
    filename = function() {
      file.path(paste0("snp_scan_", chr_pos(), ".pdf")) },
    content = function(file) {
      req(scan_obj())
      snp_w <- req(scan_window())
      phenos <- req(phename())
      pdf(file)
      ## Plots over all phenotypes
      print(top_pat_plot(phenos, group = "pheno",
                         top_type = "SNP allele pattern"))
      print(top_pat_plot(phenos, group = "pattern",
                         top_type = "Phenotype"))
      ## Plots by phenotype.
      for(pheno in phenos) {
        print(top_pat_plot(pheno, FALSE))
        top_snp_asso(pheno, phenos, snp_w)
      }
      dev.off()
    }
  )
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
