#' Shiny peaks selection
#'
#' Shiny module for peaks selection. Does not yet have return.
#'
#' @param input,output,session standard shiny arguments
#' @param pheno_type,peaks_tbl,pmap_obj reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return list of inputs and scan summary
#'
#' @export
shinyPeaks <- function(input, output, session,
                       pheno_type, peaks_tbl, pmap_obj) {
  ns <- session$ns

  ## Need to construct simplifed scan_obj with 1 trait.
  ## For that will need probs_obj(), which takes time to load.
  ## Maybe better to just do simplified plot here.
  ## Or could create dummy and saveRDS.

  # Drop-down selection box for which phenotype set
  output$choose_peakset <- renderUI({
    checkboxGroupInput(ns("peak_set"), "Phenotype Group",
                       choices = as.list(pheno_type()),
                       selected="", inline=TRUE)
  })
  # Select chromosome
  output$peak_chr <- renderUI({
    selectInput(ns("peak_chr_id"), strong("Key Chromosome"),
                choices = as.list(c("all",names(pmap_obj()))),
                selected = input$peak_chr_id, 
                multiple = TRUE)
    })

  scan_obj_all <- reactive({
    out_peaks <- pmap_obj()
    map_chr <- names(out_peaks)
    n_chr <- length(map_chr)
    peak_window <- req(input$peak_window)
    pheno_types <- pheno_type()
    out_peaks <- list(map=out_peaks,
                      lod= matrix(0, length(unlist(out_peaks)),
                                  length(pheno_types),
                                  dimnames = list(NULL, pheno_types)))

    peaks <- peaks_tbl()
    index1 <- 1
    withProgress(message = 'Hot peak windowing ...', value = 0,
    {
      for(chri in map_chr) {
        incProgress(1/n_chr, detail = paste("chr", chri))
        mapi <- out_peaks$map[[chri]]
        index2 <- index1 + length(mapi) - 1
        for(phenoj in pheno_types) {
          if(phenoj=="all")
            posi <- peaks
          else
            posi <- peaks %>%
              filter(pheno_type == phenoj)
          posi <- (posi %>% filter(chr==chri))$pos
          out_peaks$lod[seq(index1, index2),phenoj] <-
            apply(outer(posi, mapi,
                        function(x,y,z) abs(x-y) <= z,
                        peak_window),
                  2, sum)
        }
        index1 <- index2 + 1
      }
    })
    out_peaks
  })
  scan_obj <- reactive({
    peak_chr <- req(input$peak_chr_id)
    out_peaks <- scan_obj_all()
    withProgress(message = 'Hot peak filtering ...', value = 0,
    {
      setProgress(1)
      map_chr <- names(out_peaks$map)
      if(!("all" %in% peak_chr)) {
        ## reduce to included chromosomes. Not quite there yet.
        len <- sapply(out_peaks$map, length)
        index <- split(seq_len(nrow(out_peaks$lod)),
                       ordered(rep(map_chr, times=len), map_chr))
        index <- unlist(index[map_chr %in% peak_chr])
        out_peaks$lod <- out_peaks$lod[index,]
        out_peaks$map <- out_peaks$map[map_chr %in% peak_chr]
      }
    })
    class(out_peaks) <- c("scan1", class(out_peaks))
    out_peaks
  })

  output$peak_show <- renderPlot({
    ## Here need scan1 object.
    peak_set <- req(input$peak_set)
    out_peaks <- scan_obj()
    withProgress(message = 'Hot peak plot ...',
                 value = 0,
    {
      setProgress(1)
      ## Build up count of number of peaks
      pheno_types <- pheno_type()
      lodcolumns <- match(peak_set, pheno_types)
      col <- seq_along(pheno_types)
      names(col) <- pheno_types
      plot_scan1(out_peaks, lodcolumn=lodcolumns[1],
                 col=col[lodcolumns[1]],
                 ylab="phenotype count",
                 ylim=c(0,max(out_peaks$lod[,lodcolumns])))
      if(length(lodcolumns) > 1) {
        for(lodi in lodcolumns[-1])
          plot_scan1(out_peaks, lodcolumn=lodi,
                     add=TRUE, col=col[lodi])
      }
      ## Want to add mtext for peak_set
      title(paste0("number of ", paste(peak_set, collapse=","),
                   " in ", input$peak_window, "Mbp window"))
    })
  })
  scan_tbl <- reactive({
    lodcol <- match(req(input$peak_set), pheno_type())
    scan <- scan_obj()
    withProgress(message = 'Hot peak summary ...', value = 0,
    {
      setProgress(1)
      chr <- names(scan$map)
      summary(scan, lodcol, chr)
    })
  })
  output$peak_tbl <- renderTable({scan_tbl()})
  reactive({
    out <- scan_tbl()
    if(!is.null(window <- input$peak_window))
      out <- cbind(out,window=window)
    out
  })
}

#' UI for shinyPeaks input
#'
#' UI for peak selection to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyPeaks
#' @export
shinyPeaksInput <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("choose_peakset")),
    uiOutput(ns("peak_chr")),
    sliderInput(ns("peak_window"), "Peak Window", 0, 6, 3))
}
#' UI for shinyPeaks output
#'
#' UI for peak selection to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyPeaks
#' @export
shinyPeaksOutput <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("peak_show")),
    tableOutput(ns("peak_tbl"))
  )
}
