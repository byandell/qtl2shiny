#' Shiny phenotype selection
#'
#' Shiny module for phenotype selection.
#'
#' @param input,output,session standard shiny arguments
#' @param setup_par,peaks_tbl,analyses_tbl,chr_peak reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyPhenos <- function(input, output, session,
                        setup_par, peaks_tbl, analyses_tbl,
                        chr_peak) {
  ns <- session$ns

  ## Set up analyses data frame.
  analyses_df <- reactive({
    ## Filter by dataset.
    dataset <- req(setup_par$dataset)
    dat <- analyses_tbl()
    if(!("all" %in% dataset)) {
      dat <- dat %>%
        filter(pheno_type %in% dataset)
    }

    ## Filter by Peak Position if use_pos=TRUE
    if(isTruthy(input$use_pos)) {
      chr_id <- req(chr_peak$chr_id)
      peak_Mbp <- req(chr_peak$peak_Mbp)
      window_Mbp <- req(chr_peak$window_Mbp)
      if(window_Mbp > 0) {
        ## Filter peaks
        peaks <- peaks_tbl() %>%
          filter(chr == chr_id,
                 pos >= peak_Mbp - window_Mbp,
                 pos <= peak_Mbp + window_Mbp)
        dat <- dat[dat$output %in% peaks$output,]
      }
    }
    dat
  })

  # Choose Phenotypes for Analysis.
  output$pheno_names <- renderUI({
    phenames <- sort(analyses_df()$pheno)

    selected <- input$pheno_names
    if("all" %in% selected)
      selected <- c(selected[!(selected %in% c("all","none"))],
                    phenames)
    if("none" %in% selected)
      selected <- ""
    if(!is.null(selected))
      selected <- sort(unique(selected))

    ## Update phenames to include selected (but not "")
    phenames <- unique(c(selected, sort(phenames)))
    phenames <- phenames[phenames != ""]

    choices <- c("all","none", phenames)
    selectInput(ns("pheno_names"), "Choose phenotypes",
                choices = choices,
                selected = selected,
                multiple = TRUE)
  })
  output$filter <- renderUI({
    checkboxInput(ns("use_pos"),
                  paste("Peak on chr", chr_peak$chr_id, "in",
                        paste(chr_peak$peak_Mbp + c(-1,1) * chr_peak$window_Mbp,
                              collapse = "-")))
  })
  input
}
#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinyPhenos
#' @export
shinyPhenosUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("pheno_names")),
    uiOutput(ns("filter"))
  )
}
