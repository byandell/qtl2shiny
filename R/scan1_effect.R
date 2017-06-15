scan1_effect <- function(genoprobs, pheno, kinship, addcovar,
                   sex_type, blups) {
      switch(sex_type,
             "F" = addcovar <- addcovar[addcovar[,"sex"] == 0,, drop = FALSE],
             "M" = addcovar <- addcovar[addcovar[,"sex"] == 1,, drop = FALSE])
      qtl2pattern::listof_scan1coef(genoprobs, pheno, kinship, addcovar,
                                    blups)
}