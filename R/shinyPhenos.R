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

    ## Filter by Peak Position.
    chr_id <- chr_peak$chr_id
    if(input$use_pos & !is.null(setup_par$window_Mbp)) {
      if(setup_par$window_Mbp > 0) {
        ## Filter peaks
        peaks <- peaks_tbl() %>%
          filter(chr == chr_id,
                 pos >= chr_peak$peak_Mbp - setup_par$window_Mbp,
                 pos <= chr_peak$peak_Mbp + setup_par$window_Mbp)
        dat <- dat[dat$output %in% peaks$output,]
      }
    }
    dat
  })

  # Choose Phenotypes for Analysis.
  output$choose_phenoanal <- renderUI({
    phenames <- sort(analyses_df()$output)

    selected <- input$pheno_anal
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
    selectInput(ns("pheno_anal"), "Choose phenotypes",
                choices = choices,
                selected = selected,
                multiple = TRUE)
  })
  reactive({
    pheno_anal <- req(input$pheno_anal)
    dat <- analyses_tbl() %>%
      filter(output %in% pheno_anal) %>%
      arrange(output) %>%
      select(pheno,output)
    if(nrow(dat)) {
      pheno <- dat$pheno
      names(pheno) <- dat$output
      pheno
    } else { ## all, none or some non-pheno
      pheno_anal
    }
  })
}
#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinyPhenos
#' @export
shinyPhenosUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("choose_phenoanal")),
    checkboxInput(ns("use_pos"), "Filter Traits by Peak Region")
  )
}
