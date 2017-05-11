expr_region <- function(chr_id, scan_window, datapath, covar, qtls, pmap) {
  # Get expression mRMNA measurements.
  out <- DOread::read_mrna(chr_id,
                           scan_window[1], scan_window[2],
                           datapath, qtl = TRUE)
  # Covariate matrix covar is global.
  out$cov_med <- covar[, c("sex", paste0("DOwave", 2:4))]
  
  if(qtls == 2) {
    
    out[[2]]$driver <- qtl2geno::find_marker(pmap, chr_id, out[[2]]$qtl_pos)
  }
  
  out
}
