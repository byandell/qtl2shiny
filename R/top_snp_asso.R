#' Top SNP Association Plot Wrapper
#'
#' @param pheno name of phenotype for effect scan
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param xlim x limits for plot
#' @param snp_action character string for plot
#' @param drop.hilit hilight SNPs within this value of max
#' @param phename names of phenotypes
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' 
top_snp_asso <- function(scan_obj, xlim,
                         snp_action="basic",
                         phename = dimnames(scan_obj$lod)[[2]],
                         drop.hilit = NULL) {
  if(is.null(pheno) | is.null(scan_obj) | is.null(xlim))
    return(print(plot_null()))
  
  chr_id <- names(scan_obj$map)[1]
  if(is.null(drop.hilit)) {
    drop.hilit <- min(1.5,
                      max(scan_obj$lod) - 3)
  }
  facet <- ifelse(length(phename) > 1, "pheno", NULL)
  plot(scan_obj, seq_along(phename),
       show_all_snps = FALSE,
       drop.hilit = drop.hilit,
       xlim = xlim,
       col = col, col.hilit = col.hilit,
       cex = cex,
       facet = facet)
}
