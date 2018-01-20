scan1_effect <- function(genoprobs, pheno, kinship, addcovar,
                   sex_type, blups) {
  addcovar <- qtl2pattern::sexcovar(addcovar, sex_type)
  addcovar <- qtl2pattern::covar_df_mx(addcovar)
  qtl2pattern::listof_scan1coef(genoprobs, pheno, kinship, addcovar,
                                blups)
}