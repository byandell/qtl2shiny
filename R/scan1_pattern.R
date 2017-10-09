pull_patterns <- function(patterns, pheno_names) {
  if(all(pheno_names %in% patterns$pheno)) {
    dplyr::filter(patterns, pheno %in% pheno_names)
  }
  else {
    out <- dplyr::filter(patterns, pheno == "AddSex")
    if(nrow(out))
      out <- dplyr::mutate(out, pheno = pheno_names[1])
    out
  }
}
scan1_pattern <- function(pheno, phe_df, addcovar, probs36_obj, K_chr, analyses_df,
                          pats, sex_type, blups) {
  analyses_df <- which_covar(analyses_df)
  wh <- match(pheno, names(phe_df))
  
  addcovar <- wh_covar(analyses_df, wh, addcovar)
  addcovar <- covar_df_mx(addcovar)
  
  qtl2pattern::scan_pattern(probs36_obj$probs,
                            phe_df[,, drop=FALSE],
                            K_chr, addcovar,
                            probs36_obj$map,
                            pats,
                            blups = blups)
}