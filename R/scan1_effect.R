scan1_effect <- function(genoprobs, pheno, kinship, addcovar,
                   sex_type, blups) {
  switch(sex_type,
         "F" = addcovar <- addcovar[addcovar$sex == "F",, drop = FALSE],
         "M" = addcovar <- addcovar[addcovar$sex == "M",, drop = FALSE])
  addcovar <- covar_df_mx(addcovar)
  qtl2pattern::listof_scan1coef(genoprobs, pheno, kinship, addcovar,
                                blups)
}