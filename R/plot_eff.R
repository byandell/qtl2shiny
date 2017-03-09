#' Plot effects
#'
#' Plot of effects.
#'
#' @param pheno name of phenotype for effect scan
#' @param eff_obj object of class \code{\link[qtl2pattern]{listof_scan1coef}}
#' @param map map object
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param xlim x limits for plot
#' @param chr_id chromosome ID
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#' 
plot_eff <- function(pheno, eff_obj, map, scan_obj, xlim = NULL,
                     addlod = FALSE) {
  if(is.null(eff_obj) | is.null(pheno) | is.null(scan_obj))
    return(NULL)
  
  main <- pheno
  lodcol <- match(pheno, names(eff_obj))
  chr_id <- names(map)[1]
  if(!addlod) {
    max_pos <- max(scan_obj, map, lodcolumn=lodcol)$pos[1]
    plot(eff_obj[[lodcol]], map,
         xlim=xlim,
         main = main,
         legend.position = "right",
         maxpos = max_pos, maxcol = lodcol)
  } else { # coef_and_lod
    plot(eff_obj[[lodcol]], map,
         scan1_output = subset(scan_obj, map,
                               chr = chr_id,
                               lodcolumn = lodcol),
         xlim = xlim,
         maxcol = lodcol,
         main = main)
  }
}
