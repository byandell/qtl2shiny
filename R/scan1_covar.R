#' scan1 for multiple traits with possibly different covariates
#' 
#' @importFrom qtl2scan scan1
#' @importFrom dplyr select
#' 
scan1_covar <- function(probs_obj, phe_df, K_chr, cov_mx, analyses_df) {
  ## Covariate indicators follow winsorize column.
  wh <- which("winsorize" == names(analyses_df))
  is_covar <- apply(analyses_df[, -(seq_len(wh)), drop=FALSE], 2, any)
  ## Keep only covariate indicators with at least one TRUE value.
  analyses_df <- analyses_df[, names(is_covar)[is_covar], drop=FALSE]
  ## Collapse to unique identifier for each row = each phenotype.
  covarset <- apply(analyses_df, 1, function(x) paste(1 * x, collapse = ""))
  ucov <- unique(covarset)
  wh <- which(covarset == ucov[1])
  scanfn <- function(probs_obj, phe_df, K_chr, cov_mx, analyses_df, wh) {
    # scan1 for wh phenotypes using their covariates.
    covars <- unlist(analyses_df[wh[1],])
    if(length(covars)) {
      qtl2scan::scan1(probs_obj, 
                      phe_df[, wh, drop=FALSE], 
                      K_chr, 
                      cov_mx[, covars, drop=FALSE])
    } else { # no covariates (unlikely)
      qtl2scan::scan1(probs_obj, 
                      phe_df[, wh, drop=FALSE], 
                      K_chr)
    }
  }
  scans <- scanfn(probs_obj, phe_df, K_chr, cov_mx, analyses_df, wh)
  if(length(ucov) > 1) for(i in ucov[-1]) {
    wh <- which(covarset == i)
    scans <- cbind(scans, 
                   scanfn(probs_obj, phe_df, K_chr, cov_mx, analyses_df, wh))
  }
  scans
}