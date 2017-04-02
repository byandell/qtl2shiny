#' scan1 of pattern when multiple traits may have possibly different covariates
#' 
#' @importFrom qtl2pattern scan_pattern
#' 
scan1_pattern <- function(pheno, phe_df, cov_mx, probs36_obj, K_chr, analyses_df,
                          pats, blups) {
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