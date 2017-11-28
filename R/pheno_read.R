#' @export
pheno_read <- function(pheno_data, analyses_df, transform = TRUE) {
  if(transform) {
    transform <- analyses_df$transf
  } else {
    transform <- NULL
  }
  qtl2pattern::pheno_trans(pheno_data, 
                           analyses_df$pheno, 
                           transform,
                           analyses_df$offset,
                           analyses_df$winsorize)
}