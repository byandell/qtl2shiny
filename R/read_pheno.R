read_pheno <- function(pheno_data, analyses_df, transform = TRUE) {
  ## Make sure we get only one column per distinct pheno.
  DOread::get_pheno(
    pheno_data,
    dplyr::distinct(analyses_df, pheno, .keep_all=TRUE),
    transform)
}