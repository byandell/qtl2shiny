scan1_effect <- function(genoprobs, pheno, kinship, addcovar,
                   sex_type, blups) {
  addcovar <- sexcovar(addcovar, sex_type)
  addcovar <- covar_df_mx(addcovar)
  qtl2ggplot::listof_scan1coef(genoprobs, pheno, kinship, addcovar,
                                blups)
}