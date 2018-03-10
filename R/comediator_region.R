comediator_region <- function(pheno_name, chr_id, scan_window, 
                              covar, analyses_tbl, pheno_data, peaks, 
                              qtls = 2, pmap) {
  
  # Covariates used by comediators.
  covars <- colnames(covar)
  covars <- covars[!is.na(match(covars, colnames(analyses_tbl)))]
  # Replace any NA with FALSE.
  analyses_tbl[, covars] <- apply(analyses_tbl[, covars], 2, 
                                  function(x) ifelse(is.na(x), FALSE, x))
  # Kludge to get names of covariates that are used by comediators.
  covars <- apply(analyses_tbl[, covars], 2, any)
  covars <- names(covars)[covars]
  
  # This is specific to CCmouse.
  peaks <- dplyr::filter(peaks,
                         !(pheno_group == "Islet.mRNA"))
  
  pheno_data <- qtl2pattern::pheno_trans(
    pheno_data, analyses_tbl$pheno, analyses_tbl$transf,
    analyses_tbl$offset, analyses_tbl$winsorize)
  
  out <- qtl2pattern::pheno_region(
    chr_id, scan_window[1], scan_window[2], covar[, covars], pmap,
    peaks, analyses_tbl, pheno_data, drivers = qtls)
  out$phe_type = (dplyr::filter(peaks, pheno == pheno_name))$pheno_type[1]
  out
}
comediator_type <- function(comed, doThis) {
  if(doThis) {
    not_type <- comed$annot$id[comed$annot$biotype != comed$phe_type]
    comed$comediators <- comed$comediators[, not_type]
  }

  comed
}
