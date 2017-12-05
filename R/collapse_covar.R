# Collapse covariate names in analyses table
#
# Collapse covariates used to comma-separated list and return table
collapse_covar <- function(analyses_tbl) {
  if(is.null(analyses_tbl))
    return(NULL)
  
  ## Collapse covariates (past winsorize column).
  covar_names <- names(analyses_tbl)[-seq_len(match("winsorize",
                                                    names(analyses_tbl)))]
  dplyr::select(
    dplyr::mutate(
      tidyr::unite(analyses_tbl, covar, 
                   dplyr::one_of(covar_names)), 
      covar = sapply(strsplit(covar,"_"),
                     function(x) paste(covar_names[as.logical(x)],
                                       collapse=","))),
    pheno, covar, transf, offset, winsorize)
}
