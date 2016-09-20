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
shinyHotspot <- function(input, output, session,
                         set_par, win_par, 
                         pheno_type, peaks_tbl, pmap_obj) {
  ns <- session$ns

  chr_names <- reactive({
    names(req(pmap_obj()))
  })
  # Hotspot Search (if desired)
  # Select chromosome.
  output$chr_ct <- renderUI({
    choices <- chr_names()
    if(is.null(selected <- input$chr_ct))
      selected <- "all"
    selectInput(ns("chr_ct"), strong("Chr"),
                choices = c("all", choices),
                selected = selected)
  })
  scan_obj_all <- reactive({
    req(input$hotspot)
    map_chr <- chr_names()
    n_chr <- length(map_chr)
    peak_window <- req(win_par$window_Mbp)
    if(is.null(peak_window)) {
      NULL
    } else {
      pheno_types <- pheno_type()
      map <- pmap_obj()
      out_peaks <- list(map=map,
                        lod= matrix(0, length(unlist(map)),
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
    }
  })
  
  scan_obj <- reactive({
    out_peaks <- scan_obj_all()
    withProgress(message = 'Hotspot search ...', value = 0,
    {
      setProgress(1)
      map_chr <- names(out_peaks$map)
      chr_ct <- input$chr_ct
      if(!("all" %in% chr_ct)) {
        if(!is.null(chr_ct)) {
          ## reduce to included chromosomes. Not quite there yet.
          len <- sapply(out_peaks$map, length)
          index <- split(seq_len(nrow(out_peaks$lod)),
                         ordered(rep(map_chr, times=len), map_chr))
          index <- unlist(index[map_chr %in% chr_ct])
          out_peaks$lod <- out_peaks$lod[index,]
          out_peaks$map <- out_peaks$map[map_chr %in% chr_ct]
        }
      }
    })
    class(out_peaks) <- c("scan1", class(out_peaks))
    out_peaks
  })

  output$peak_show <- renderPlot({
    peak_set <- req(set_par$dataset)
    out_peaks <- scan_obj()
    withProgress(message = 'Hotspot show ...',
                 value = 0, {
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
                   " in ", win_par$window_Mbp, "Mbp window"))
    })
  })
  scan_tbl <- reactive({
    lodcol <- match(req(set_par$dataset), pheno_type())
    scan <- scan_obj()
    withProgress(message = 'Hotspot summary ...', value = 0, {
      setProgress(1)
      chr <- names(scan$map)
      summary(scan, lodcol, chr)
    })
  })
  output$peak_tbl <- renderTable({
    out <- scan_tbl() %>%
      arrange(desc(lod))
    if(nrow(out) > 4) {
      out <- out %>%
      filter(row_number() < 4) %>%
      add_row(pheno="...")
    }
    out
  })
  ## Return.
  scan_tbl
}
#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinyHotspot
#' @export
shinyHotspotInput <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(6, checkboxInput(ns("hotspot"), "Search Hotspots?")),
    column(6, uiOutput(ns("chr_ct"))))
}
#' @rdname shinyHotspot
#' @export
shinyHotspotOutput <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("peak_show")),
    tableOutput(ns("peak_tbl"))
  )
}
