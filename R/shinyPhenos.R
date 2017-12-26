#' Shiny phenotype selection
#'
#' Shiny module for phenotype selection.
#'
#' @param input,output,session standard shiny arguments
#' @param set_par,win_par,peaks_tbl,analyses_tbl,cov_df,project_info reactive arguments
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
                        set_par, win_par, peaks_tbl, analyses_tbl, cov_df,
                        project_info) {
  ns <- session$ns

  ## Set up analyses data frame.
  analyses_set <- shiny::reactive({
    shiny::req(project_info())
    ## Filter by dataset.
    dataset <- set_par$dataset
    data_group <- set_par$pheno_group

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
    dataset <- set_par$dataset
    if(shiny::isTruthy(dataset)) {
      dat <- shiny::req(analyses_set())
      ## Filter by Peak Position if use_pos=TRUE
      if(isTruthy(input$use_pos)) {
        window_Mbp <- shiny::req(win_par$window_Mbp)
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
    chr_id <- shiny::req(win_par$chr_id)
    peak_Mbp <- shiny::req(win_par$peak_Mbp)
    window_Mbp <- shiny::req(win_par$window_Mbp)
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
    shiny::req(project_info())
    cat("output$pheno_names", input$pheno_names, "\n", file = stderr())
    phenames <- selected <- input$pheno_names
    if(shiny::isTruthy(peaks_df())) {
      phenames <- unique(peaks_df()$pheno)
      # Limit to first 1000
      nphe <- length(phenames)
      phenames <- phenames[seq_len(min(1000, nphe))]
    }
    
    if("all" %in% selected)
      selected <- c(selected[!(selected %in% c("all","none"))],
                    phenames)
    if("none" %in% selected)
      selected <- ""
    if(!is.null(selected)) {
      selected <- sort(unique(selected))
      selected <- selected[selected %in% phenames]
    }
    
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
  
  # Output the peaks table
  output$peaks_tbl <- shiny::renderDataTable({
    dplyr::arrange(
      dplyr::select(
        peaks_df(), pheno, chr, pos, lod),
      dplyr::desc(lod))
  }, options = list(scrollX = TRUE, pageLength = 10))
  
  ## Density or scatter plot of phenotypes.
  analyses_plot <- shiny::reactive({
    shiny::req(analyses_tbl())
    phename <- shiny::req(input$pheno_names)
    dplyr::filter(analyses_tbl(), pheno %in% phename)
  })
  phe_mx <- shiny::reactive({
    pheno_read(pheno_data(), analyses_plot())
  })
  raw_phe_mx <- shiny::reactive({
    pheno_read(pheno_data(), analyses_plot(), FALSE)
  })
  shiny::callModule(shinyPhenoPlot, "PhenoPlotRaw", raw_phe_mx, cov_df)
  shiny::callModule(shinyPhenoPlot, "PhenoPlotTrans", phe_mx, cov_df)
  
  # Show data.
  output$radio <- renderUI({
    shiny::radioButtons(ns("radio"), NULL,
                        c("LOD Peaks","Covariates",
                          "Trans Data","Raw Data"),
                        input$radio)
  })
  output$show_data <- renderUI({
    switch(shiny::req(input$radio),
           "LOD Peaks"  = shiny::dataTableOutput(ns("peaks_tbl")),
           "Raw Data"   = shinyPhenoPlotUI(ns("PhenoPlotRaw")),
           "Trans Data" = shinyPhenoPlotUI(ns("PhenoPlotTrans")),
           "Covariates" = shiny::dataTableOutput(ns("analyses_tbl")))
  })
  
  # Output the analyses table
  output$analyses_tbl <- shiny::renderDataTable({
    collapse_covar(analyses_plot())
  }, options = list(scrollX = TRUE, pageLength = 10))
  
  output$filter <- shiny::renderUI({
    shiny::checkboxInput(ns("use_pos"),
                  paste0("Peak on chr ", win_par$chr_id, " in ",
                        paste(win_par$peak_Mbp + c(-1,1) * win_par$window_Mbp,
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
    shiny::uiOutput(ns("pheno_names")),
    shiny::uiOutput(ns("radio"))
  )
}
#' @rdname shinyPhenos
#' @export
shinyPhenosOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::dataTableOutput(ns("pheno_lod")),
    shiny::uiOutput(ns("show_data"))
  )
}
