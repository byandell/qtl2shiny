#' Shiny peaks selection
#'
#' Shiny module for peaks selection. Does not yet have return.
#'
#' @param input,output,session standard shiny arguments
#' @param set_par,pheno_type,peaks_tbl,pmap_obj reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return list of inputs and scan summary
#'
#' @export
#' 
#' @importFrom dplyr distinct filter
#' @importFrom shiny callModule NS req 
#'   selectInput sliderInput updateSelectInput updateSliderInput textInput
#'   uiOutput
#'   renderUI
#'   observeEvent
#'   fluidRow column strong tagList
#'   
shinyPeaks <- function(input, output, session,
                       set_par, pheno_type, peaks_tbl, pmap_obj) {
  ns <- session$ns

  # Select chromosome. Defaults to blank.
  output$chr_id <- shiny::renderUI({
    choices <- names(pmap_obj())
    selected <- input$chr_id
    if(!isTruthy(selected))
      choices <- c("", choices)
    shiny::selectInput(ns("chr_id"), shiny::strong("Chromosome"),
                choices = choices,
                selected = selected)
  })
  
  ## Window slider
  output$window_Mbp <- shiny::renderUI({
    if(is.null(pos <- input$window_Mbp))
      pos <- 5
    shiny::sliderInput(ns("window_Mbp"), "Window Half Width",
                0, 10, pos, step=0.5)
  })
  
  # Peak position slider.
  output$peak_Mbp <- shiny::renderUI({
    chr_id <- shiny::req(input$chr_id)
    rng <- round(range(pmap_obj()[[chr_id]]), 2)
    pos <- input$peak_Mbp
    if(is.null(pos)) {
      pos <- mean(rng)
    } else {
      if(pos < rng[1] | pos > rng[2])
        pos <- mean(rng)
    }
    shiny::sliderInput(ns("peak_Mbp"), "Peak Position (Mbp)", rng[1], rng[2], pos,
                step=.5)
  })
  
  ## shorthand 
  output$chr_pos <- shiny::renderUI({
    shiny::textInput(ns("chr_pos"), "Position", input$chr_pos)
  })
  
  scan_tbl <- shiny::callModule(shinyHotspot, "hotspot",
              set_par, input, pheno_type, peaks_tbl, pmap_obj)
  
  shiny::observeEvent(scan_tbl(), update_region())
  shiny::observeEvent(input$chr_id, update_region())
  update_region <- function() {
    scan_in <- shiny::req(scan_tbl())
    scan_in <- dplyr::filter(scan_in, lod==max(lod))
    
    chr_ct <- as.character(scan_in$chr)
    if(any(chr_ct == shiny::req(input$chr_id))) {
      chr_ct <- as.character(input$chr_id)
      scan_in <- dplyr::filter(scan_in, chr == chr_ct)
      pmap <- shiny::req(pmap_obj())
      peak_Mbp <- scan_in$pos[1]
      rng <- round(range(pmap[[chr_ct]]), 2)
      shiny::updateSliderInput(session, "peak_Mbp",
                        min=rng[1], max=rng[2],
                        value=peak_Mbp)
#      shiny::updateTextInput(session, "chr_pos",
#                      value = paste(chr_ct, round(peak_Mbp, 2), sep="@"))
    }
  }
  shiny::observeEvent(input$chr_pos, {
#    chr_pos <- strsplit(input$chr_pos, ":|@|_| |,")[[1]]
#    if(length(chr_pos) == 2) {
#      chr <- chr_pos[1]
#      pmap <- pmap_obj()
#      choices <- names(pmap)
#      if(chr %in% choices) {
#        pos <- as.numeric(chr_pos[2])
    chr <- shiny::req(input$chr_id)
    pos <- as.numeric(input$chr_pos)
    if(is.numeric(pos)) {
      if(!is.na(pos)) {
        rng <- round(range(pmap[[chr]]), 2)
        if(pos >= rng[1] & pos <= rng[2]) {
          up <- is.null(input$chr_id)
          if(!up) {
            up <- (chr != input$chr_id)
          }
          if(up) {
            shiny::updateSelectInput(session, "chr_id",
                              selected=chr, choices=choices)
          }
          shiny::updateSliderInput(session, "peak_Mbp",
                            value=pos, min=rng[1], max=rng[2])
        }
      }
    }
  })
  ## Return.
  input
}
#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinyPeaks
#' @export
shinyPeaksInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::checkboxInput(ns("local"), "Local Scan in Window?", TRUE),
    shiny::fluidRow(
      shiny::column(6, shiny::uiOutput(ns("chr_id"))),
      shiny::column(6, shiny::uiOutput(ns("chr_pos")))
    ),
    shiny::uiOutput(ns("peak_Mbp")),
    shiny::uiOutput(ns("window_Mbp")),
    shinyHotspotInput(ns("hotspot")))
}
#' @rdname shinyPeaks
#' @export
shinyPeaksOutput <- function(id) {
  ns <- shiny::NS(id)
  shinyHotspotOutput(ns("hotspot"))
}
