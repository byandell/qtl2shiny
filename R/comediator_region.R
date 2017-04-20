comediator_region <- function(pheno_name, chr_id, scan_window, covar, analyses_tbl, pheno_data) {
  
  out <- list(phe_group = (dplyr::filter(peaks,
                              pheno == pheno_name))$pheno_group[1],
              annot = dplyr::filter(peaks,
                                chr == chr_id,
                                pos >= scan_window[1],
                                pos <= scan_window[2],
                                pheno != pheno_name))
  out$annot <- 
    dplyr::left_join(out$annot,
                     analyses_tbl,
                     by = c("pheno", "longname", "output",
                            "pheno_group", "pheno_type"))
  # Covariates used by comediators.
  covars <- colnames(covar)
  # Replace any NA with FALSE.
  out$annot[, covars] <- apply(out$annot[, covars], 2, 
                                      function(x) ifelse(is.na(x), FALSE, x))
  # Kludge to get names of covariates that are used by comediators.
  covars <- apply(out$annot[, covars], 2, any)
  covars <- names(covars)[covars]
  
  out$cov_med <- covar[, covars]
  
  out$comediators <- DOread::get_pheno(pheno_data, out$annot)
  
  out
}
comediator_group <- function(comed) {
  not_group <- comed$annot$pheno[comed$annot$pheno_group != comed$phe_group]
  comed$comediators <- comed$comediators[, not_group]

  comed
}
