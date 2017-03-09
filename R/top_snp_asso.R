#' Top SNP Association Plot Wrapper
#'
#' @param pheno name of phenotype for effect scan
#' @param map map object
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param xlim x limits for plot
#' @param snp_action character string for plot
#' @param drop.hilit hilight SNPs within this value of max
#' @param phename names of phenotypes
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' 
top_snp_asso <- function(scan_obj, snpinfo, xlim,
                         snp_action="basic",
                         phename = dimnames(scan_obj)[[2]],
                         drop.hilit = NULL,
                         show_all_snps = FALSE) {
  if(is.null(phename) | is.null(scan_obj) | is.null(xlim))
    return(print(plot_null()))
  
  if(is.null(drop.hilit)) {
    drop.hilit <- min(1.5,
                      max(unclass(scan_obj)) - 3)
  }
  facet <- NULL
  facet <- if(length(phename) > 1) 
    facet <- "pheno"
  
  plot(scan_obj, snpinfo, seq_along(phename),
       show_all_snps = show_all_snps,
       drop.hilit = drop.hilit,
       xlim = xlim,
       facet = facet)
}
