#' Shiny chromosome and peak selection
#'
#' Shiny module for chromosome and peak selection.
#'
#' @param input,output,session standard shiny arguments
#' @param hot_peak,pmap_obj reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyChrPeak <- function(input, output, session,
                        hot_peak=reactive({NULL}), pmap_obj) {
  ns <- session$ns

  # Select chromosome.
  output$choose_chr <- renderUI({
    choices <- names(pmap_obj())
    selected <- input$chr_id
    selectInput(ns("chr_id"), strong("Chromosome"),
                choices = choices,
                selected = selected,
                multiple = TRUE)
  })

  observeEvent(hot_peak(), {
    browser()
    scan_tbl <- hot_peak() %>%
      filter(lod==max(lod),
             row_number() ==1)

    chr_id <- as.character(scan_tbl$chr[1])
    pmap <- pmap_obj()
    updateSelectInput(session, "chr_id",
                      choices = names(pmap),
                      selected = chr_id)
    peak_Mbp <- scan_tbl$pos[1]
    rng <- round(range(pmap[[chr_id]]), 2)
    updateSliderInput(session, "peak_Mbp",
                      min=rng[1], max=rng[2],
                      value=peak_Mbp)
    updateTextInput(session, "chr_pos",
                    value = paste(chr_id, round(peak_Mbp, 2), sep="@"))
  })

  observeEvent(input$chr_pos, {
    chr_pos <- strsplit(input$chr_pos, ":|@|_| |,")[[1]]
    if(length(chr_pos) == 2) {
      chr <- chr_pos[1]
      pmap <- pmap_obj()
      choices <- names(pmap)
      if(chr %in% choices) {
        pos <- as.numeric(chr_pos[2])
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
      }
    }
  })
  input
}
#' @param id identifier for \code{\link{shinyScan1}} use
#' @rdname shinyChrPeak
#' @export
shinyChrPeakUI <- function(id) {
  ns <- NS(id)
  tagList(
    textInput(ns("chr_pos"), "chr:pos or chr@pos"),
    uiOutput(ns("choose_chr")),
    uiOutput(ns("choose_peak"))
  )
}
