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
#' @importFrom dplyr add_row arrange filter rename
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
  output$hotspot <- shiny::renderUI({
    shiny::checkboxInput(ns("hotspot"), "Search Hotspots?", input$hotspot)
  })
  # Select chromosome.
  output$chr_ct <- shiny::renderUI({
    choices <- chr_names()
    if(is.null(selected <- input$chr_ct))
      selected <- "all"
    shiny::selectInput(ns("chr_ct"), strong("Chr"),
                choices = c("all", choices),
                selected = selected,
                multiple = TRUE)
  })
  shiny::observeEvent(input$chr_ct, {
    is_all <- grep("all", input$chr_ct)
    if(length(is_all)) {
      if(length(input$chr_ct) > 1) {
        selected <- input$chr_ct[-is_all]
        choices <- chr_names()
        shiny::updateSelectInput(session, "chr_ct", strong("Chr"),
                                 choices = c("all", choices),
                                 selected = selected)
      }
    }
  })
  scan_obj_all <- shiny::reactive({
    shiny::req(input$hotspot, win_par$window_Mbp, input$minLOD)
    shiny::withProgress(message = 'Hotspot scan ...', value = 0,
    {
      shiny::setProgress(1)
      hotspot(pmap_obj(), peaks_tbl(), 2^win_par$window_Mbp, input$minLOD)
    })
  })
  
  scan_obj <- shiny::reactive({
    out_peaks <- scan_obj_all()
    shiny::withProgress(message = 'Hotspot search ...', value = 0,
    {
      shiny::setProgress(1)
      chr_ct <- input$chr_ct
      if(!("all" %in% chr_ct)) {
        out_peaks <- subset(out_peaks, chr_ct)
      }
    })
    out_peaks
  })

  output$peak_show <- shiny::renderPlot({
    if(shiny::isTruthy(set_par$dataset)) {
      peak_set <- set_par$dataset
    } else {
      peak_set <- "all"
    }
    map <- scan_obj()$map
    out_peaks <- scan_obj()$scan
    shiny::withProgress(message = 'Hotspot show ...',
                 value = 0, {
      shiny::setProgress(1)
      ## Build up count of number of peaks
      pheno_types <- colnames(out_peaks)
      lodcolumns <- match(peak_set, pheno_types)
      lodcolumns <- lodcolumns[!is.na(lodcolumns)]
      col <- seq_along(pheno_types)
      names(col) <- pheno_types
      nchr <- length(map)
      xaxt <- ifelse(nchr < 5, "y", "n")
      traits <- ifelse(length(peak_set) == 1, peak_set, "traits")
      
      # Kludge 
      if(nrow(out_peaks) == length(unlist(map)) & length(lodcolumns)) {
        plot(out_peaks, map, lodcolumn=lodcolumns,
             col = col[lodcolumns],
             ylab = "phenotype count",
             ylim = c(0,max(out_peaks[,lodcolumns])),
             xaxt = xaxt,
             gap = 25 / nchr) +
          ## add mtext for peak_set
          ggplot2::ggtitle(paste0("number of ", traits,
                                  " in ", 2 ^ win_par$window_Mbp, "Mbp window"))
      } else {
        plot_null("no data")
      }
    })
  })
  scan_tbl <- shiny::reactive({
    # Used chosen datasets, or all if not chosen.
    if(shiny::isTruthy(set_par$dataset)) {
      peak_set <- make.names(set_par$dataset)
    } else {
      peak_set <- "all"
    }
    map <- scan_obj()$map
    scan <- scan_obj()$scan
    # Match lod columns to those present.
    lodcol <- match(peak_set, make.names(colnames(scan)))
    lodcol <- lodcol[!is.na(lodcol)]
    if(length(lodcol) & (nrow(scan) == length(unlist(map)))) {
      shiny::withProgress(message = 'Hotspot summary ...', value = 0, {
        shiny::setProgress(1)
        chr <- names(map)
        dplyr::select(
          dplyr::rename(
            summary(
              subset(scan, lodcolumn = lodcol),
              map, chr = chr),
            count = lod),
          -marker)
      })
    } else {
      NULL
    }
  })
  output$peak_tbl <- shiny::renderDataTable({
    shiny::req(scan_tbl(), peaks_tbl())
    peakDataTable(scan_tbl(), peaks_tbl())
  }, escape = FALSE,
  options = list(pageLength = 5))
  
  # Minimum LOD for SNP top values.
  minLOD <- reactive({
    if(shiny::isTruthy(input$minLOD)) {
      input$minLOD
    } else {
      max(3, round(min(shiny::req(peaks_tbl())$lod), 1))
    }
  })
  output$minLOD <- shiny::renderUI({
    value <- minLOD()
    shiny::numericInput(ns("minLOD"), "LOD threshold", value, min = 0, step = 0.5)
  })
  
  ## Return.
  scan_tbl
}
#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinyHotspot
#' @export
shinyHotspotInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(6, shiny::uiOutput(ns("hotspot"))),
    shiny::column(6, shiny::uiOutput(ns("chr_ct"))),
    shiny::uiOutput(ns("minLOD")))
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
