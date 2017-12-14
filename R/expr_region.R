expr_region <- function(chr_id, scan_window, covar, qtls, pmap,
                        project_info) {
  # Get expression mRMNA measurements.
  query_mrna <- read_query_rds(project_info, "query_mrna.rds")
  out <- query_mrna(chr_id, scan_window[1], scan_window[2], qtl = TRUE)

  m <- match(tolower(c("sex", "DOwave")), colnames(covar), nomatch = 0)
  if(any(m == 0))
    cat("sex and DOwave not found in data", file = stderr())
  
  out$cov_med <- covar[, m]
  
  if(qtls == 2) {
    
    out[[2]]$driver <- qtl2geno::find_marker(pmap, chr_id, out[[2]]$qtl_pos)
  }
  
  out
}
