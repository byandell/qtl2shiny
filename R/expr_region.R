expr_region <- function(chr_id, scan_window, datapath, covar) {
  # Get expression mRMNA measurements.
  out <- DOread::read_mrna(chr_id,
                           scan_window[1], scan_window[2],
                           datapath)
  # Covariate matrix covar is global.
  out$cov_med <- covar[, c("sex", paste0("DOwave", 2:4))]
  out
}
