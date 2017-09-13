comediator_region <- function(pheno_name, chr_id, scan_window, 
                              covar, analyses_tbl, pheno_data, peaks, 
                              qtls = 1, pmap) {
  
  # Annotation for phenotypes.
  annot <- dplyr::filter(peaks, 
                         chr == chr_id,
                         pos >= scan_window[1],
                         pos <= scan_window[2],
                         pheno != pheno_name)
  annot <- dplyr::left_join(annot,
                            analyses_tbl,
                            by = c("pheno", "longname", "output",
                                   "pheno_group", "pheno_type"))
  
  # Covariates used by comediators.
  covars <- colnames(covar)
  covars <- covars[!is.na(match(covars, colnames(annot)))]
  
  # Replace any NA with FALSE.
  annot[, covars] <- apply(annot[, covars], 2, 
                           function(x) ifelse(is.na(x), FALSE, x))
  
  # Add QTL peaks
  annot <- 
    dplyr::inner_join(
      annot,
      dplyr::ungroup(
        dplyr::summarize(
          dplyr::group_by(peaks, pheno),
          qtl_ct = n(),
          QTL = paste0(chr_id, "@",
                       round(pos), ":",
                       round(lod), collapse = ","))),
      by = "pheno")
  
  # Kludge to get names of covariates that are used by comediators.
  covars <- apply(annot[, covars], 2, any)
  covars <- names(covars)[covars]
  
  if(qtls == 2)
    annot$driver <- qtl2geno::find_marker(pmap, chr_id, annot$pos)
  
  list(comediators = DOread::get_pheno(pheno_data, annot),
       annot = dplyr::rename(annot, 
                             id = pheno,
                             biotype = pheno_type),
       cov_med = covar[, covars],
       phe_type = (dplyr::filter(peaks, pheno == pheno_name))$pheno_type[1])
}
comediator_type <- function(comed, doThis) {
  if(doThis) {
    not_type <- comed$annot$id[comed$annot$biotype != comed$phe_type]
    comed$comediators <- comed$comediators[, not_type]
  }

  comed
}
