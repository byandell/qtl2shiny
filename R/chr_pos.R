#' Make chr_pos from chr, peak and window
#'
#' Create character string of chr:left-right with left=peak-window, right=peak+window
#'
#' @param chr_id chromosome ID
#' @param peak_Mbp position of peak in Mbp
#' @param window_Mbp half-width of window around peak in Mbp
#' @param range left and right window positions in Mbp
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{make_chr_pos(chr_id, peak_Mbp, window_Mbp)}
#' 
#' @export
make_chr_pos <- function(chr_id=NULL, peak_Mbp=NULL, window_Mbp=NULL,
                         range = c(peak_Mbp - 2 ^ window_Mbp,
                                   peak_Mbp + 2 ^ window_Mbp)) {
  if(is.null(chr_id))
    chr_id <- "?"
  if(length(range) < 2 )
    range <- rep("?", 2)
  paste(chr_id, range[1], range[2], sep = "_")
}
