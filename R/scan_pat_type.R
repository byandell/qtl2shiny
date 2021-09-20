# Plot scan patterns for phenotype.
#' @importFrom rlang .data
#' @importFrom dplyr arrange desc filter
#' @importFrom qtl2pattern sdp_to_pattern
#' 
scan_pat_type <- function(scan_pat, map, type, pattern, pheno, haplos) {
  if(is.null(scan_pat))
    return(plot_null())
  
  pattern_cont <- 
    unique(
      dplyr::filter(
        scan_pat$patterns,
        qtl2pattern::sdp_to_pattern(.data$sdp, haplos) %in% pattern)$founders)
  # plot at most 8 curves
  lodcolumn <- seq_len(min(8, length(pattern_cont)))
  maxpos <- NULL
  maxcol <- 1
  if(type == "coef") {
    maxpos <- (dplyr::arrange(
      summary(scan_pat$scan, map),
      dplyr::desc(.data$lod)))$pos[1]
  }
  title <- pheno
  if(ncol(scan_pat$scan) > 1)
    title <- FALSE
  ggplot2::autoplot(
    scan_pat, map, type, pattern_cont, main = title,
    maxpos = maxpos, maxcol = maxcol, lodcolumn = lodcolumn) 
}
