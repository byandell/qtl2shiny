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
  analyses_set <- shiny::reactive({
    ## Filter by dataset.
    dataset <- setup_par$dataset
    data_group <- setup_par$pheno_group

    dat <- analyses_tbl()

    # Start with phenotypes in data groups.
    if(shiny::isTruthy(data_group)) {
      dat <- dplyr::filter(dat, pheno_group %in% data_group)
    }
    # For identified datasets, drop other datasets from those data groups.
    if(shiny::isTruthy(dataset)) {
      if(!("all" %in% dataset)) {
        dat_sets <- dplyr::distinct(dat, pheno_type, pheno_group)
        dat_groups <- unique(dplyr::filter(dat_sets,
                                    pheno_type %in% dataset)$pheno_group)
        dat <- dplyr::filter(dat, 
                             (pheno_type %in% dataset) |
                             !(pheno_group %in% dat_groups))
      }
    }
    dat
  })
  analyses_df <- reactive({
    dataset <- setup_par$dataset
    if(shiny::isTruthy(dataset)) {
      dat <- shiny::req(analyses_set())
      ## Filter by Peak Position if use_pos=TRUE
      if(isTruthy(input$use_pos)) {
        window_Mbp <- shiny::req(chr_peak$window_Mbp)
        ## Infinite recursion here.
        if(shiny::isTruthy(peaks_df()))
          dat <- dat[dat$pheno %in% peaks_df()$pheno,]
      }
    }
    dat
  })

  # Choose Phenotypes for Analysis.
  peaks_df <- reactive({
    shiny::req(analyses_set(), peaks_tbl())
    chr_id <- shiny::req(chr_peak$chr_id)
    peak_Mbp <- shiny::req(chr_peak$peak_Mbp)
    window_Mbp <- shiny::req(chr_peak$window_Mbp)
    peaks_in_pos(analyses_set(), peaks_tbl(),
                 shiny::isTruthy(input$use_pos),
                 chr_id, peak_Mbp, window_Mbp)
  })
  output$pheno_lod <- shiny::renderDataTable({
    shiny::req(peaks_df())
  }, escape = FALSE,
  options = list(scrollX = TRUE, pageLength = 10,
                 lengthMenu = c(5,10,25)))
  
  output$pheno_names <- shiny::renderUI({
    phenames <- selected <- input$pheno_names
    if(shiny::isTruthy(peaks_df())) {
      phenames <- peaks_df()$pheno
      # Limit to first 1000
      nphe <- length(phenames)
      phenames <- phenames[seq_len(min(1000, nphe))]
    }
    
    if("all" %in% selected)
      selected <- c(selected[!(selected %in% c("all","none"))],
                    phenames)
    if("none" %in% selected)
      selected <- ""
    if(!is.null(selected))
      selected <- sort(unique(selected))

    ## Update phenames to include selected (but not "")
    phenames <- unique(c(selected, phenames))
    phenames <- phenames[phenames != ""]

    choices <- c("all","none", phenames)
    label = ifelse(nphe <= 1000,
                   "Choose phenotypes",
                   paste("Top 1000 of", nphe))
    shiny::selectInput(ns("pheno_names"), label,
                choices = choices,
                selected = selected,
                multiple = TRUE)
  })
  output$filter <- shiny::renderUI({
    shiny::checkboxInput(ns("use_pos"),
                  paste0("Peak on chr ", chr_peak$chr_id, " in ",
                        paste(chr_peak$peak_Mbp + c(-1,1) * chr_peak$window_Mbp,
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
#' @rdname shinyPhenos
#' @export
shinyPhenosOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::dataTableOutput(ns("pheno_lod"))
}
