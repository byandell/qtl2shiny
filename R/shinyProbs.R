#' Shiny Probability module
#'
#' Shiny genotype probability access.
#' 
#' @param input,output,session standard shiny arguments
#' @param win_par reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyProbs <- function(input, output, session,
                       win_par) {
  ns <- session$ns
  
  reactive({
    chr_id <- req(win_par()$chr_id)
    withProgress(message = 'Read probs ...', value = 0, {
      setProgress(1)
      read_probs(chr_id, datapath)
    })
  })
}
#' @rdname shinyProbs
#' @export
shinyProbs36 <- function(input, output, session,
                       win_par) {
  ns <- session$ns

  range_val <- reactive({
    req(win_par()$peak_Mbp, win_par()$window_Mbp)
    c(win_par()$peak_Mbp + c(-1,1) * win_par()$window_Mbp)
  })

  ## Probs object for 36 diplotypes.
  reactive({
    chr_id <- req(win_par()$chr_id)
    withProgress(message = 'Diplotype Probs ...', value = 0, {
      setProgress(1)
      read_probs36(chr_id, range_val()[1], range_val()[2],
                   datapath)
    })
  })
}
