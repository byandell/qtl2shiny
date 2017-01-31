#' Plot effects
#'
#' Plot of effects.
#'
#' @param pheno name of phenotype for effect scan
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param eff_obj object of class \code{\link[doqtl2]{listof_scan1coefCC}}
#' @param xlim x limits for plot
#' @param main plot title
#' @param ... parameters past to plot routine
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#' 
#' @importFrom qtl2ggplot plot_coefCC
#' @importFrom ggplot2 geom_vline
#'
plot_eff <- function(pheno, scan_obj, eff_obj, xlim = NULL,
                     main = pheno, ...) {
  if(is.null(scan_obj) | is.null(eff_obj) | is.null(pheno))
    return(NULL)
  
  lodcol <- match(pheno, names(eff_obj))
  max_pos <- max(scan_obj, lodcolumn=lodcol)$pos[1]
  qtl2ggplot::plot_coefCC(eff_obj[[lodcol]], xlim=xlim,
              main = main,
              legend.position = "right", ...) +
    ggplot2::geom_vline(xintercept=max_pos, linetype=2,
               col=lodcol)
}
