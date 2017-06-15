pull_patterns <- function(patterns, pheno_name) {
  if(pheno_name %in% patterns$pheno) {
    dplyr::filter(patterns, pheno == pheno_name)
  }
  else {
    out <- dplyr::filter(patterns, pheno == "AddSex")
    if(nrow(out))
      out <- dplyr::mutate(out, pheno = pheno_name)
    out
  }
}
scan1_pattern <- function(pheno, phe_df, cov_mx, probs36_obj, K_chr, analyses_df,
                          pats, sex_type, blups) {
  analyses_df <- which_covar(analyses_df)
  wh <- match(pheno, names(phe_df))
  covars <- unlist(analyses_df[wh,])
  
  ## Names of haplos and diplos in terms of founders.
  diplos <- dimnames(probs36_obj$probs[[1]])[[2]]
  haplos <-  unique(unlist(stringr::str_split(diplos, "")))
  
  qtl2pattern::scan_pattern(probs36_obj$probs,
                            phe_df[,pheno, drop=FALSE],
                            K_chr, cov_mx[, covars, drop=FALSE],
                            probs36_obj$map,
                            pats, haplos, diplos,
                            blups = blups)
}