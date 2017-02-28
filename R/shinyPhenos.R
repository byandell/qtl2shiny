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
#' @importFrom dplyr filter
#' @importFrom shiny NS reactive req 
#'   checkboxInput selectInput
#'   uiOutput
#'   renderUI
#'   tagList
#'   withProgress setProgress
#'   downloadButton downloadHandler
shinyPhenos <- function(input, output, session,
                        setup_par, peaks_tbl, analyses_tbl,
                        chr_peak) {
  ns <- session$ns

  ## Set up analyses data frame.
  analyses_df <- shiny::reactive({
    ## Filter by dataset.
    dataset <- shiny::req(setup_par$dataset)
    dat <- analyses_tbl()
    if(!("all" %in% dataset)) {
      dat <- dplyr::filter(dat, pheno_type %in% dataset)
    }

    ## Filter by Peak Position if use_pos=TRUE
    if(isTruthy(input$use_pos)) {
      chr_id <- shiny::req(chr_peak$chr_id)
      peak_Mbp <- shiny::req(chr_peak$peak_Mbp)
      window_Mbp <- 2 ^ shiny::req(chr_peak$window_Mbp)
      if(window_Mbp > 0) {
        ## Filter peaks
        peaks <- dplyr::filter(peaks_tbl(), chr == chr_id,
                 pos >= peak_Mbp - window_Mbp,
                 pos <= peak_Mbp + window_Mbp)
        dat <- dat[dat$output %in% peaks$output,]
      }
    }
    dat
  })

  # Choose Phenotypes for Analysis.
  output$pheno_names <- shiny::renderUI({
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
    shiny::selectInput(ns("pheno_names"), "Choose phenotypes",
                choices = choices,
                selected = selected,
                multiple = TRUE)
  })
  output$filter <- shiny::renderUI({
    shiny::checkboxInput(ns("use_pos"),
                  paste0("Peak on chr ", chr_peak$chr_id, " in ",
                        paste(chr_peak$peak_Mbp + c(-1,1) * 2 ^ chr_peak$window_Mbp,
                              collapse = "-"), "?"),
                  TRUE)
  })
  input
}
#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinyPhenos
#' @export
shinyPhenosUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("filter")),
    shiny::uiOutput(ns("pheno_names"))
  )
}
