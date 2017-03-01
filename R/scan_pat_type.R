#' Plot scan patterns for phenotype.
#' 
#' @importFrom dplyr arrange desc
#' 
scan_pat_type <- function(scan_pat, type, pattern, pheno) {
  if(is.null(scan_pat))
    return(plot_null())
  
  pattern_cont <- 
    unique(dplyr::filter(scan_pat$patterns,
                         CCSanger::sdp_to_pattern(sdp) %in% pattern)$founders)
  maxpos <- NULL
  maxcol <- 1
  if(type == "coef") {
    maxpos <- (dplyr::arrange(
      summary(scan_pat$scan),
      dplyr::desc(lod)))$pos[1]
  }
  plot(scan_pat, type, pattern_cont, main = pheno,
       maxpos = maxpos, maxcol = maxcol) 
}
