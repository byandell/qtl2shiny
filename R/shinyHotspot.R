#' Shiny hotspot views
#'
#' Shiny module to view hotspots for peak selection.
#'
#' @param input,output,session standard shiny arguments
#' @param set_par,win_par,pheno_type,peaks_tbl,pmap_obj reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return list of inputs and scan summary
#'
#' @export
#' @importFrom dplyr add_row arrange filter
#' @importFrom ggplot2 ggtitle
#' @importFrom shiny NS reactive req 
#'   checkboxInput selectInput
#'   plotOutput dataTableOutput uiOutput
#'   renderPlot renderDataTable renderUI
#'   fluidRow column tagList
#'   withProgress incProgress setProgress
shinyHotspot <- function(input, output, session,
                         set_par, win_par, 
                         pheno_type, peaks_tbl, pmap_obj) {
  ns <- session$ns

  chr_names <- shiny::reactive({
    names(shiny::req(pmap_obj()))
  })
  # Hotspot Search (if desired)
  # Select chromosome.
  output$chr_ct <- shiny::renderUI({
    choices <- chr_names()
    if(is.null(selected <- input$chr_ct))
      selected <- "all"
    shiny::selectInput(ns("chr_ct"), strong("Chr"),
                choices = c("all", choices),
                selected = selected)
  })
  scan_obj_all <- shiny::reactive({
    shiny::req(input$hotspot)
    map_chr <- chr_names()
    n_chr <- length(map_chr)
    peak_window <- shiny::req(win_par$window_Mbp)
    if(is.null(peak_window)) {
      NULL
    } else {
      peak_window <- 2 ^ peak_window
      pheno_types <- pheno_type()
      map <- pmap_obj()
      out_peaks <- matrix(0, length(unlist(map)),
                                  length(pheno_types),
                                  dimnames = list(qtl2scan:::map2markernames(map), pheno_types))
      
      peaks <- peaks_tbl()
      index1 <- 1
      shiny::withProgress(message = 'Hot peak windowing ...', value = 0,
                   {
                     for(chri in map_chr) {
                       shiny::incProgress(1/n_chr, detail = paste("chr", chri))
                       mapi <- map[[chri]]
                       index2 <- index1 + length(mapi) - 1
                       for(phenoj in pheno_types) {
                         if(phenoj=="all")
                           posi <- peaks
                         else
                           posi <- dplyr::filter(peaks, pheno_type == phenoj)
                         posi <- dplyr::filter(posi, chr==chri)$pos
                         out_peaks[seq(index1, index2),phenoj] <-
                           apply(outer(posi, mapi,
                                       function(x,y,z) abs(x-y) <= z,
                                       peak_window),
                                 2, sum)
                       }
                       index1 <- index2 + 1
                     }
                   })
      class(out_peaks) <- c("scan1", "matrix")
      out_peaks
    }
  })
  
  scan_obj <- shiny::reactive({
    out_peaks <- scan_obj_all()
    map <- pmap_obj()
    shiny::withProgress(message = 'Hotspot search ...', value = 0,
    {
      shiny::setProgress(1)
      chr_ct <- input$chr_ct
      if(!("all" %in% chr_ct)) {
        out_peaks <- subset(out_peaks, map, chr_ct)
      }
    })
    out_peaks
  })

  output$peak_show <- shiny::renderPlot({
    peak_set <- shiny::req(set_par$dataset)
    map <- pmap_obj()
    chr_ct <- shiny::req(input$chr_ct)
    if(!("all" %in% chr_ct)) {
      map <- map[shiny::req(chr_names()) %in% chr_ct]
    }
    out_peaks <- scan_obj()
    shiny::withProgress(message = 'Hotspot show ...',
                 value = 0, {
      shiny::setProgress(1)
      ## Build up count of number of peaks
      pheno_types <- pheno_type()
      lodcolumns <- match(peak_set, pheno_types)
      col <- seq_along(pheno_types)
      names(col) <- pheno_types
      plot(out_peaks, map, lodcolumn=lodcolumns,
           col = col[lodcolumns],
           ylab = "phenotype count",
           ylim = c(0,max(out_peaks[,lodcolumns]))) +
        ## add mtext for peak_set
        ggplot2::ggtitle(paste0("number of ", paste(peak_set, collapse=","),
                   " in ", 2 ^ win_par$window_Mbp, "Mbp window"))
    })
  })
  scan_tbl <- shiny::reactive({
    peak_set <- shiny::req(set_par$dataset)
    lodcol <- match(shiny::req(set_par$dataset), pheno_type(), pmap_obj())
    scan <- scan_obj()
    shiny::withProgress(message = 'Hotspot summary ...', value = 0, {
      shiny::setProgress(1)
      chr <- sapply(pmap_obj(), 
                    function(map, nms) any(nms %in% names(map)), 
                    rownames(scan))
      chr <- names(chr)[chr]
      summary(subset(scan, lodcolumn = peak_set), pmap_obj(), chr = chr)
    })
  })
  output$peak_tbl <- shiny::renderDataTable({
    dplyr::arrange(scan_tbl(), desc(lod))
  }, escape = FALSE,
  options = list(pageLength = 5))
  ## Return.
  scan_tbl
}
#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinyHotspot
#' @export
shinyHotspotInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(6, shiny::checkboxInput(ns("hotspot"), "Search Hotspots?")),
    shiny::column(6, shiny::uiOutput(ns("chr_ct"))))
}
#' @rdname shinyHotspot
#' @export
shinyHotspotOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::plotOutput(ns("peak_show")),
    shiny::dataTableOutput(ns("peak_tbl"))
  )
}
