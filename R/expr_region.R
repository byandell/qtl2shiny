expr_region <- function(chr_id, scan_window, covar, qtls, pmap) {
  # Get expression mRMNA measurements.
  out <- query_mrna(chr_id, scan_window[1], scan_window[2], qtl = TRUE)
  # Covariate matrix covar is global.
  out$cov_med <- covar[, c("sex", "DOwave")]
  
  if(qtls == 2) {
    
    out[[2]]$driver <- qtl2geno::find_marker(pmap, chr_id, out[[2]]$qtl_pos)
  }
  
  out
}
