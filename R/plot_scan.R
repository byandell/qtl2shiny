plot_scan <- function(scan_obj, map, lodcolumn, chr_id, xlim, phe_mx) {
  p <- ggplot2::autoplot(
    scan_obj, map,
    lodcolumn = lodcolumn,
    chr = chr_id,
    xlim = xlim)
  if(is.null(p)) {
    plot_null()
  } else{
    if(ncol(phe_mx) == 1 & ncol(scan_obj) >= 1)
      p <- p + ggtitle(colnames(phe_mx))
    p    
  }
}