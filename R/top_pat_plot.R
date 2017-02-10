#' Top Pattern Plot
#'
#' Plot of top patterns.
#'
#' @param pheno name of phenotype for effect scan
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param xlim x limits for plot
#' @param fill.null null plot if no data
#' @param group group by one of c("pheno","pattern")
#' @param snp_action character string for plot
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @examples
#' \dontrun{top_pat_plot(pheno, scan_obj, xlim)}
#'
#' @importFrom qtl2pattern topsnp_pattern
#' 
top_pat_plot <- function(pheno,
                         scan_obj, xlim,
                         fill.null=TRUE,
                         top_pattern = qtl2pattern::topsnp_pattern(scan_obj, pheno),
                         group = "pheno",
                         snp_action = "basic") {
  if(is.null(top_pattern)) {
    if(fill.null)
      return(plot_null())
    else
      return()
  }
  chr_id <- names(scan_obj$map)[1]
  p <- plot(top_pattern, group=group) + xlim(xlim)
  if(length(pheno) == 1) {
    mytitle <- paste(pheno, "chr", chr_id)
    if(snp_action != "basic")
      mytitle <- paste(mytitle, snp_action)
    p <- p + ggtitle(mytitle)
  }
  p
}
