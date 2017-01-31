#' Top SNP Association Plot Wrapper
#'
#' @param pheno name of phenotype for effect scan
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param xlim x limits for plot
#' @param snp_action character string for plot
#' @param drop.hilit hilight SNPs within this value of max
#' @param main plot title
#' @param pch,cex,col,col.hilit plot parameters
#' @param ... additional plot parameters
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#'
top_snp_asso <- function(pheno, scan_obj, xlim,
                         snp_action="basic",
                         drop.hilit = NULL,
                         main = mytitle, pch = 1, cex = 0.75,
                         col = "#7570b3", col.hilit = "#d95f02", ...) {
  if(is.null(pheno) | is.null(scan_obj) | is.null(xlim))
    return(print(plot_null()))
  
  phename <- dimnames(scan_obj$lod)[[2]]
  chr_id <- names(scan_obj$map)[1]
  scan_pheno <- subset(scan_obj,
                       lodcolumn=match(pheno, phename))
  if(is.null(drop.hilit)) {
    drop.hilit <- min(1.5,
                      max(scan_pheno$lod) - 3)
  }
  mytitle <- paste(pheno, "chr", chr_id)
  if(snp_action != "basic")
    mytitle <- paste(mytitle, snp_action)
  plot_snpasso(scan_pheno,
               show_all_snps = FALSE,
               drop.hilit = drop.hilit,
               xlim = xlim,
               col = col, col.hilit = col.hilit,
               pch = pch, cex = cex,
               main = main, ...)
}
