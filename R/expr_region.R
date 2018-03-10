expr_region <- function(chr_id, scan_window, covar, qtls, pmap,
                        project_info) {
  # Get expression mRMNA measurements.
  query_mrna <- read_query_rds(project_info, "query_mrna.rds")
  
  qtl2pattern::expr_region(chr_id, scan_window[1], scan_window[2], covar, pmap, 
                           drivers = qtls,
                           query_mrna = query_mrna)
}
