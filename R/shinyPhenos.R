#' Shiny phenotype selection
#'
#' Shiny module for phenotype selection.
#'
#' @param input,output,session standard shiny arguments
#' @param pheno_type,peaks_tbl,analyses_tbl,hot_peak,win_par reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyPhenos <- function(input, output, session,
                        pheno_type, peaks_tbl, analyses_tbl,
                        hot_peak, win_par) {
  ns <- session$ns

  # Drop-down selection box for which phenotype set
  output$choose_dataset <- renderUI({
    req(choices <- pheno_type())
    if(is.null(selected <- input$dataset))
      selected <- choices[2]
    selectInput(ns("dataset"), "Phenotype Group",
                choices = as.list(choices),
                selected = selected)
  })
  observeEvent(hot_peak(), {
    scan_tbl <- hot_peak() %>%
      filter(lod==max(lod))
    updateSelectInput(session, "dataset",
                      choices = pheno_type(),
                      selected = scan_tbl$pheno[1])
  })

  # Text box for selecting phenotypes
  output$choose_analysis <- renderUI({
    selectInput(ns("analysis"), "Analysis Type",
                c("default","anal1","anal2","all"),
                "all")
  })

  ## Set up analyses data frame.
  analyses_df <- reactive({
    dataset <- req(input$dataset)

    ## Filter by dataset.
    dat <- analyses_tbl()
    if(dataset != "all") {
      dat <- dat %>%
        filter(pheno_type == dataset)
    }

    ## Filter by phename (partial match to phenotype names).
#    phename <- input$phename
#    dat <- dat %>%
#      filter(grepl(tolower(phename), tolower(pheno)))

    ## Filter by analysis type.
    analysis_type <- req(input$analysis)
    if(analysis_type == "default") {
      dat <- dat %>%
        filter(!grepl(paste(paste0("_anal", 1:2), collapse="|"), output))
    } else {
      if(analysis_type != "all")
        dat <- dat %>%
          filter(grepl(paste0("_", analysis_type), output))
    }

    ## Filter by Peak Position.
    chr_id <- win_par$chr_id
    if(input$use_pos & !is.null(win_par$window_Mbp)) {
      if(win_par$window_Mbp > 0) {
        ## Filter peaks
        peaks <- peaks_tbl() %>%
          filter(chr == chr_id,
                 pos >= win_par$peak_Mbp - win_par$window_Mbp,
                 pos <= win_par$peak_Mbp + win_par$window_Mbp)
        dat <- dat[dat$output %in% peaks$output,]
      }
    }
    dat
  })

  ## Update if dataset changes.
  observeEvent(input$dataset, {
    updateTextInput(session, "phename", value="")
    updateSelectInput(session, "analysis", selected="all")
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
    uiOutput(ns("choose_dataset")),
    uiOutput(ns("choose_analysis")),
    uiOutput(ns("choose_phenoanal")),
    checkboxInput(ns("use_pos"), "Filter Traits by Peak Region")
  )
}
