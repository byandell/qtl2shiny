#' scan1 of pattern when multiple traits may have possibly different covariates
#' 
#' @importFrom qtl2pattern scan_pattern
#' 
scan1_pattern <- function(pheno, phe_df, cov_mx, probs36_obj, K_chr, analyses_df,
                          pats, haplos, diplos) {
  analyses_df <- which_covar(analyses_df)
  wh <- match(pheno, names(phe_df))
  covars <- unlist(analyses_df[wh,])

  qtl2pattern::scan_pattern(probs36_obj,
                            phe_df[,pheno, drop=FALSE],
                            K_chr, cov_mx[, covars, drop=FALSE],
                            pats, haplos, diplos)
}