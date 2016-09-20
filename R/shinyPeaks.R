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
                       setup_par, pheno_type, peaks_tbl, pmap_obj) {
  ns <- session$ns

  # Select chromosome.
  output$chr_id <- renderUI({
    choices <- names(pmap_obj())
    selected <- input$chr_id
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
      pmap <- pmap_obj()
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
#      }
#    }
  })
  
  # Hotspot Search (if desired)
  # Select chromosome.
  output$chr_ct <- renderUI({
    choices <- names(pmap_obj())
    if(is.null(selected <- input$chr_ct))
      selected <- "all"
    selectInput(ns("chr_ct"), strong("Chr"),
                choices = c("all", choices),
                selected = selected)
  })
  scan_obj_all <- reactive({
    req(input$hotspot)
    out_peaks <- pmap_obj()
    map_chr <- names(out_peaks)
    n_chr <- length(map_chr)
    peak_window <- req(input$window_Mbp)
    if(is.null(peak_window)) {
      NULL
    } else {
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
    ## Here need scan1 object.
    peak_set <- req(setup_par$dataset)
    out_peaks <- scan_obj()
    withProgress(message = 'Hotspot show ...',
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
                   " in ", input$window_Mbp, "Mbp window"))
    })
  })
  scan_tbl <- reactive({
    lodcol <- match(req(setup_par$dataset), pheno_type())
    scan <- scan_obj()
    withProgress(message = 'Hotspot summary ...', value = 0,
    {
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
    fluidRow(
      column(6, checkboxInput("hotspot", "Search Hotspots?")),
      column(6, uiOutput(ns("chr_ct"))))
    )
}
#' @rdname shinyPeaks
#' @export
shinyPeaksOutput <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("peak_show")),
    tableOutput(ns("peak_tbl"))
  )
}
