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
shinyPeaks <- function(input, output, session,
                       set_par, pheno_type, peaks_tbl, pmap_obj) {
  ns <- session$ns

  # Select chromosome. Defaults to blank.
  output$chr_id <- renderUI({
    choices <- names(pmap_obj())
    selected <- input$chr_id
    if(!isTruthy(selected))
      choices <- c("", choices)
    selectInput(ns("chr_id"), strong("Chromosome"),
                choices = choices,
                selected = selected)
  })
  
  ## Window slider
  output$window_Mbp <- renderUI({
    if(is.null(pos <- input$window_Mbp))
      pos <- 5
    sliderInput(ns("window_Mbp"), "Window Half Width",
                0, 10, pos, step=0.5)
  })
  
  # Peak position slider.
  output$peak_Mbp <- renderUI({
    chr_id <- req(input$chr_id)
    rng <- round(range(pmap_obj()[[chr_id]]), 2)
    pos <- input$peak_Mbp
    if(is.null(pos)) {
      pos <- mean(rng)
    } else {
      if(pos < rng[1] | pos > rng[2])
        pos <- mean(rng)
    }
    sliderInput(ns("peak_Mbp"), "Peak Position (Mbp)", rng[1], rng[2], pos,
                step=.5)
  })
  
  ## shorthand 
  output$chr_pos <- renderUI({
    textInput(ns("chr_pos"), "Position", input$chr_pos)
  })
  
  scan_tbl <- callModule(shinyHotspot, "hotspot",
              set_par, input, pheno_type, peaks_tbl, pmap_obj)
  
  observeEvent(scan_tbl(), update_region())
  observeEvent(input$chr_id, update_region())
  update_region <- function() {
    scan_in <- req(scan_tbl())
    scan_in <- scan_in %>%
      filter(lod==max(lod))
    
    chr_ct <- as.character(scan_in$chr)
    if(any(chr_ct == req(input$chr_id))) {
      chr_ct <- as.character(input$chr_id)
      scan_in <- scan_in %>%
        filter(chr == chr_ct)
      pmap <- req(pmap_obj())
      peak_Mbp <- scan_in$pos[1]
      rng <- round(range(pmap[[chr_ct]]), 2)
      updateSliderInput(session, "peak_Mbp",
                        min=rng[1], max=rng[2],
                        value=peak_Mbp)
#      updateTextInput(session, "chr_pos",
#                      value = paste(chr_ct, round(peak_Mbp, 2), sep="@"))
    }
  }
  observeEvent(input$chr_pos, {
#    chr_pos <- strsplit(input$chr_pos, ":|@|_| |,")[[1]]
#    if(length(chr_pos) == 2) {
#      chr <- chr_pos[1]
#      pmap <- pmap_obj()
#      choices <- names(pmap)
#      if(chr %in% choices) {
#        pos <- as.numeric(chr_pos[2])
    chr <- req(input$chr_id)
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
            updateSelectInput(session, "chr_id",
                              selected=chr, choices=choices)
          }
          updateSliderInput(session, "peak_Mbp",
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
  ns <- NS(id)
  tagList(
    fluidRow(
      column(6, uiOutput(ns("chr_id"))),
      column(6, uiOutput(ns("chr_pos")))
    ),
    uiOutput(ns("peak_Mbp")),
    uiOutput(ns("window_Mbp")),
    shinyHotspotInput(ns("hotspot")))
}
#' @rdname shinyPeaks
#' @export
shinyPeaksOutput <- function(id) {
  ns <- NS(id)
  shinyHotspotOutput(ns("hotspot"))
}
