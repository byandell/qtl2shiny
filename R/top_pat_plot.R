#' Top Pattern Plot
#'
#' Plot of top patterns.
#'
#' @param pheno name of phenotype for effect scan
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param xlim x limits for plot
#' @param fill.null null plot if no data
#' @param facet facet by one of c("pheno","pattern")
#' @param snp_action character string for plot
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @examples
#' \dontrun{top_pat_plot(pheno, scan_obj, xlim)}
#'
top_pat_plot <- function(pheno,
                         scan_obj, xlim,
                         facet = "pheno",
                         snp_action = "basic") {
  chr_id <- names(scan_obj$map)[1]
  mytitle <- FALSE
  if(length(pheno) == 1) {
    mytitle <- paste(pheno, "chr", chr_id)
    if(snp_action != "basic")
      mytitle <- paste(mytitle, snp_action)
  }
  if(length(pheno) == 1)
    facet <- NULL
  scan_obj <- subset(scan_obj, 
                     lodcolumn = match(pheno, dimnames(scan_obj$lod)[[2]]))
  
  plot(scan_obj, xlim = xlim, main = mytitle,
       patterns = "hilit", drop.hilit = 1.5,
       facet = facet)
}
