#' Shiny phenotype selection
#'
#' Shiny module for phenotype selection.
#'
#' @param input,output,session standard shiny arguments
#' @param pheno_type,peaks_tbl,analyses_tbl,chr_id,peak_Mbp,window_Mbp reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
testPhenos <- function(input, output, session,
                       pheno_type, peaks_tbl, analyses_tbl) {
  ns <- session$ns

  # Drop-down selection box for which phenotype set
  output$choose_dataset <- renderUI({
    cat(file=stderr(), "choose", paste(pheno_type(), collapse=","), "\n")
    checkboxGroupInput(ns("dataset"), "Phenotype Group",
                       choices = as.list(pheno_type()),
                       selected=pheno_type()[2],inline=TRUE)
  })

  ## Set up analyses data frame.
  analyses_df <- reactive({
    dataset <- req(input$dataset)
    cat(file=stderr(), "dataset", dataset, "\n")

    ## Filter by dataset.
    dat <- analyses_tbl()
    if(!("all" %in% dataset)) {
      dat <- dat %>%
        filter(pheno_type %in% dataset)
    }
    dat
  })

  ## Duplicate for now
  output$choose_phenoanal <- renderUI({
    selected <- input$columns
    if(!is.null(input$dataset))
      choices <- analyses_df()$output
    else
      choices <- NULL
    choices <- unique(c("all",selected, choices))
    checkboxGroupInput(ns("columns"), "Choose phenotypes",
                       choices = choices,
                       selected=selected,
                       inline=TRUE)
  })
  columns_re <- reactive({
    columns <- req(input$columns)
    browser()
    columns
  })
}

#' UI for testPhenos Shiny Module
#'
#' UI for phenotype selection to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname testPhenos
#' @export
testPhenosUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("choose_dataset")),
    uiOutput(ns("choose_phenoanal"))
  )
}
