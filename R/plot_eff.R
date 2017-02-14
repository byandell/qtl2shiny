#' Plot effects
#'
#' Plot of effects.
#'
#' @param pheno name of phenotype for effect scan
#' @param eff_obj object of class \code{\link[qtl2pattern]{listof_scan1coef}}
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param xlim x limits for plot
#' @param chr_id chromosome ID
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#' 
#' @importFrom qtl2ggplot plot_coefCC
#' @importFrom ggplot2 geom_vline
#'
plot_eff <- function(pheno, eff_obj, scan_obj, xlim = NULL,
                     addlod = FALSE) {
  if(is.null(eff_obj) | is.null(pheno) | is.null(scan_obj))
    return(NULL)
  
  main <- pheno
  lodcol <- match(pheno, names(eff_obj))
  max_pos <- max(scan_obj, lodcolumn=lodcol)$pos[1]
  chr_id <- names(scan_obj$map)[1]
  if(!addlod) {
    plot(eff_obj[[lodcol]], 
         xlim=xlim,
         main = main,
         legend.position = "right") +
      ggplot2::geom_vline(xintercept=max_pos, linetype=2,
                          col=lodcol)
  } else { # coef_and_lod
    plot(eff_obj[[lodcol]], 
         scan1_output = subset(scan_obj, 
                               chr = chr_id,
                               lodcolumn = lodcol),
         xlim = xlim,
         maxcol = lodcol,
         main = main)
  }
}
