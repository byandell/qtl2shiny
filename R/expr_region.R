expr_region <- function(chr_id, scan_window, covar, qtls, pmap,
                        project_info) {
  # Get expression mRMNA measurements.
  query_mrna <- read_query_rds(project_info, "query_mrna.rds")
  out <- query_mrna(chr_id, scan_window[1], scan_window[2], qtl = TRUE)
  if(is.null(out))
    return(NULL)

  # Check covariates
  expr_covars <- unique(out$annot$covar)
  if(length(expr_covars) > 1)
    warning("only using first type of covariate for expression")
  expr_covars <- stringr::str_split(expr_covars[1], ",")[[1]]
  m <- match(tolower(expr_covars), tolower(colnames(covar)), nomatch = 0)
  if(any(m == 0))
    cat(paste(paste(expr_covars, collapse = ","), "not found in data"), file = stderr())

  out$cov_med <- covar[, m]

  if(qtls == 2) {

    out$annot$driver <- qtl2::find_marker(pmap, chr_id, out$annot$qtl_pos)
  }

  out
}
