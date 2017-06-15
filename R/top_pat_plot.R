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
                         scan_obj,
                         chr_id,
                         map,
                         xlim,
                         drop.hilit = 1.5,
                         facet = "pheno",
                         snp_action = "basic", ...) {
  mytitle <- FALSE
  if(length(pheno) == 1) {
    mytitle <- paste(pheno, "chr", chr_id)
    if(snp_action != "basic")
      mytitle <- paste(mytitle, snp_action)
  }
  legend.title <- "pattern"
  if(length(pheno) == 1 & ncol(scan_obj) > 1) {
    pheno <- colnames(scan_obj)
  }
  if(length(pheno) == 1) {
    facet <- NULL
  } else {
    if(facet == "pattern")
      legend.title <- "pheno"
  }

  lodcol <- match(pheno, colnames(scan_obj))
  if(any(is.na(lodcol))) # Probably all sex
    lodcol <- match("AddSex", colnames(scan_obj))
  if(any(is.na(lodcol)))
    return(plot_null("scan name mismatch"))
  colnames(scan_obj)[lodcol] <- pheno
  scan_obj <- subset(scan_obj, 
                     lodcolumn = lodcol)
  
  plot(scan_obj, map, seq_along(pheno),
       xlim = xlim, main = mytitle,
       patterns = "hilit", drop.hilit = drop.hilit,
       facet = facet, legend.title = legend.title, ...)
}
