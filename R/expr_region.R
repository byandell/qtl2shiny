expr_region <- function(phe_df, chr_id, scan_window, datapath, covar) {
  indID <- rownames(phe_df)
  
  # Get expression mRMNA measurements.
  out <- DOread::read_mrna(indID, chr_id,
                           scan_window[1], scan_window[2],
                           datapath)
  # Covariate matrix covar is global.
  out$cov_med <- covar[, c("sex", paste0("DOwave", 2:4))]
  out
}
