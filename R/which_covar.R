#
which_covar <- function(analyses_df) {
  ## Covariate indicators follow winsorize column.
  wh <- which("winsorize" == names(analyses_df))
  is_covar <- apply(analyses_df[, -(seq_len(wh)), drop=FALSE], 2, any)
  ## Keep only covariate indicators with at least one TRUE value.
  analyses_df[, names(is_covar)[is_covar], drop=FALSE]
}
wh_covar <- function(analyses_df, wh, cov_df) {
  # Get which covariates from condensed analyses table.
  # The analyses table as T/F with names; capture names that have TRUE.
  covars <- unlist(analyses_df[wh[1],, drop = FALSE])
  covars <- names(covars)[covars]
  data.frame(cov_df[, covars, drop=FALSE],
             stringsAsFactors = FALSE)
}