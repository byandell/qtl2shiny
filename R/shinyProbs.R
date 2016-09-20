#' Shiny Probability module
#'
#' Shiny genotype probability access.
#' 
#' @param input,output,session standard shiny arguments
#' @param win_par,pheno_names,probs_obj reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
shinyProbs <- function(input, output, session,
                       win_par) {
  ns <- session$ns

  probs_obj <- reactive({
    chr_id <- req(win_par$chr_id)
    withProgress(message = 'Read probs ...', value = 0, {
      setProgress(1)
      read_probs(chr_id, datapath)
    })
  })
  
  probs_obj
}
#' @rdname shinyProbs
#' @export
shinyProbs36 <- function(input, output, session,
                         win_par) {
  ns <- session$ns

  ## Probs object for 36 diplotypes.
  probs_obj <- reactive({
    chr_id <- req(win_par$chr_id)
    range_val <- req(win_par$peak_Mbp) + 
      c(-1,1) * req(win_par$window_Mbp)
    withProgress(message = 'Diplotype Probs ...', value = 0, {
      setProgress(1)
      read_probs36(chr_id, range_val[1], range_val[2],
                   datapath)
    })
  })
  probs_obj
}
#' @rdname shinyProbs
#' @export
shinySNPProbs <- function(input, output, session,
                          win_par, pheno_names, probs_obj) {
  ns <- session$ns
  
  reactive({
    req(win_par$chr_id, win_par$peak_Mbp, win_par$window_Mbp)
    req(pheno_names(), probs_obj())
    withProgress(message = 'SNP Probs ...', value = 0, {
      setProgress(1)
      get_snpprobs(win_par$chr_id, 
                   win_par$peak_Mbp, 
                   win_par$window_Mbp,
                   pheno_names(), 
                   probs_obj(),
                   datapath)
    })
  })
}