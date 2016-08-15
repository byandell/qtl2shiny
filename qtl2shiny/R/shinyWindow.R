#' Shiny window selection
#'
#' Shiny module for window selection around key peak.
#'
#' @param input,output,session standard shiny arguments
#' @param pmap_obj reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyWindow <- function(input, output, session, pmap_obj,
                        hot_peak=reactive({NULL})) {
  ns <- session$ns

  # Select chromosome.
  output$choose_chr <- renderUI({
    cat(file=stderr(),"chr", input$chr_id, "\n")
    selectInput(ns("chr_id"), strong("Key Chromosome"),
                choices = names(pmap_obj()),
                selected = input$chr_id)
  })

  # Position slider
  output$choose_peak <- renderUI({
    req(input$chr_id)
    cat(file=stderr(),"peak", input$peak_Mbp, "\n")
    rng <- round(range(pmap_obj()[[input$chr_id]]), 2)
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

  ## Window slider
  output$choose_window <- renderUI({
    cat(file=stderr(),"window", input$window_Mbp, "\n")
    if(is.null(pos <- input$window_Mbp))
      pos <- 3
    sliderInput(ns("window_Mbp"), "Window over Peak if Positive",
                0, 6, pos, step=0.5)
  })
  observeEvent(hot_peak(), {
    cat(file=stderr(), "hot_peak\n")
    scan_tbl <- hot_peak() %>%
      filter(lod==max(lod))

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
    window_Mbp <- scan_tbl$window[1]
    updateSliderInput(session, "window_Mbp",
                      value = window_Mbp)
  })

  observeEvent(input$chr_pos, {
    cat(file=stderr(), "chr_pos", input$chr_pos, "\n")
    chr_pos <- strsplit(input$chr_pos, ":|@|_| |,")[[1]]
    if(length(chr_pos) == 2) {
      chr <- chr_pos[1]
      pmap <- pmap_obj()
      choices <- names(pmap)
      if(chr %in% choices) {
        pos <- as.numeric(chr_pos[2])
        if(!is.null(pos)) {
           rng <- round(range(pmap[[chr]]), 2)
          if(pos >= rng[1] & pos <= rng[2]) {
            if(chr != req(input$chr_id)) {
              updateSelectInput(session, "chr_id",
                                selected=chr, choices=choices)
            }
            updateSliderInput(session, "peak_Mbp",
                              value=pos, min=rng[1], max=rng[2])
          }
        }
      }
    }
  })
  input
}

#' UI for shinyWindow Shiny Module
#'
#' UI for window selection to use in shiny module.
#'
#' @param id identifier for \code{\link{shinyScan1}} use
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname shinyWindow
#' @export
shinyWindowUI <- function(id) {
  ns <- NS(id)
  tagList(
    textInput(ns("chr_pos"),"Quick chr:pos or chr@pos"),
    uiOutput(ns("choose_window")),
    uiOutput(ns("choose_chr")),
    uiOutput(ns("choose_peak"))
  )
}
