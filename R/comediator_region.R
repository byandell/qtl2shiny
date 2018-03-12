comediator_region <- function(pheno_name, chr_id, scan_window, 
                              covar, analyses_tbl, pheno_data, peaks, 
                              qtls = 2, pmap) {
  
  # This is specific to CCmouse.
  peaks <- dplyr::filter(peaks,
                         !(pheno_group == "Islet.mRNA"))
  
  out <- qtl2pattern::pheno_region(
    chr_id, scan_window[1], scan_window[2], covar, pmap,
    peaks, analyses_tbl, pheno_data, drivers = qtls)
  
  out
}
comediator_type <- function(comed, peaks, doThis) {
  if(doThis & !is.null(comed)) {
    phe_type <- (dplyr::filter(peaks, pheno == pheno_name))$pheno_type[1]
    not_type <- comed$annot$id[comed$annot$biotype != phe_type]
    comed$comediators <- comed$comediators[, not_type]
  }

  comed
}
