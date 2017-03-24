#' Shiny Probability module
#'
#' Shiny genotype probability access.
#' 
#' @param input,output,session standard shiny arguments
#' @param win_par,pheno_names,probs_obj,data_path reactive arguments
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom qtl2pattern get_snpprobs 
#' @importFrom DOread read_probs
#' @importFrom shiny reactive req 
#'   withProgress setProgress
shinyProbs <- function(input, output, session,
                       win_par, data_path) {
  ns <- session$ns

  probs_obj <- shiny::reactive({
    chr_id <- shiny::req(win_par$chr_id)
    shiny::withProgress(message = 'Read probs ...', value = 0, {
      shiny::setProgress(1)
      if(shiny::isTruthy(win_par$local)) {
        mid <- req(win_par$peak_Mbp)
        win <- 2 ^ req(win_par$window_Mbp)
        start_val <- mid - win
        end_val <- mid + win
      } else {
        start_val <- end_val <- NULL
      }
      # Note probs object keeps map with it
      DOread::read_probs(chr_id, start_val, end_val, 
                         datapath = data_path())
    })
  })
  
  probs_obj
}
#' @rdname shinyProbs
#' @export
shinyProbs36 <- function(input, output, session,
                         win_par, data_path) {
  ns <- session$ns

  ## Probs object for 36 diplotypes.
  probs_obj <- shiny::reactive({
    chr_id <- shiny::req(win_par$chr_id)
    range_val <- shiny::req(win_par$peak_Mbp) + 
      c(-1,1) * 2 ^ shiny::req(win_par$window_Mbp)
    shiny::withProgress(message = 'Diplotype Probs ...', value = 0, {
      shiny::setProgress(1)
      DOread::read_probs(chr_id, range_val[1], range_val[2],
                         datapath = data_path(),
                         allele = FALSE)
    })
  })
  probs_obj
}
#' @rdname shinyProbs
#' @export
shinySNPProbs <- function(input, output, session,
                          win_par, pheno_names, probs_obj,
                          data_path) {
  ns <- session$ns
  
  shiny::reactive({
    shiny::req(win_par$chr_id, win_par$peak_Mbp, win_par$window_Mbp)
    shiny::req(pheno_names(), probs_obj())
    shiny::withProgress(message = 'SNP Probs ...', value = 0, {
      shiny::setProgress(1)
      qtl2pattern::get_snpprobs(win_par$chr_id, 
                   win_par$peak_Mbp, 
                   2 ^ win_par$window_Mbp,
                   pheno_names(), 
                   probs_obj()$probs,
                   probs_obj()$map,
                   data_path())
    })
  })
}