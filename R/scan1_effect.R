scan1_effect <- function(genoprobs, pheno, kinship, addcovar,
                   sex_type, blups) {
  addcovar <- qtl2mediate::sexcovar(addcovar, sex_type)
  addcovar <- qtl2mediate::covar_df_mx(addcovar)
  qtl2ggplot::listof_scan1coef(genoprobs, pheno, kinship, addcovar,
                                blups)
}