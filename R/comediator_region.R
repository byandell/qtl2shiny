#' @importFrom rlang .data
#' 
comediator_region <- function(pheno_name, chr_id, scan_window, 
                              covar, analyses_tbl, peaks, 
                              qtls = 2, pmap, project_info) {
  
  # This is specific to CCmouse.
  peaks <- dplyr::filter(peaks,
                         .data$pheno != pheno_name,
                         !(.data$pheno_group == "Islet.mRNA"))
  
  # Filter peaks and analyses to region and drop pheno_name
  peaks_local <- dplyr::filter(peaks,
                               .data$chr == chr_id,
                               .data$pos >= scan_window[1],
                               .data$pos <= scan_window[2])
  analyses_df <- dplyr::filter(analyses_tbl,
                               .data$pheno %in% peaks_local$pheno,
                               .data$pheno != pheno_name)
  
  # Read the phenos we need.
  phenos <- analyses_df$pheno
  pheno_data <- read_project(project_info, "pheno_data", phenos)
  
  out <- pheno_region(
    chr_id, scan_window[1], scan_window[2], covar, pmap,
    peaks, analyses_tbl, pheno_data, drivers = qtls)
  
  out
}
comediator_type <- function(comed, peaks, pheno_name, doThis) {
  if(doThis & !is.null(comed)) {
    phe_type <- (dplyr::filter(peaks, .data$pheno == pheno_name))$pheno_type[1]
    not_type <- comed$annot$id[comed$annot$biotype != phe_type]
    comed$comediators <- comed$comediators[, not_type]
  }

  comed
}
