#
which_covar <- function(analyses_df) {
  ## Covariate indicators follow winsorize column.
  wh <- which("winsorize" == names(analyses_df))
  is_covar <- apply(analyses_df[, -(seq_len(wh)), drop=FALSE], 2, any)
  ## Keep only covariate indicators with at least one TRUE value.
  analyses_df[, names(is_covar)[is_covar], drop=FALSE]
}