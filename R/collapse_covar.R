#' Collapse covariate names in analyses table
#'
#' Collapse covariates used to comma-separated list and return table
#'
#' @param analyses_tbl analyses table
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{collapse_covar(analyses_tbl)}
#' 
#' @export
#' @importFrom dplyr mutate one_of select 
#' @importFrom tidyr unite
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
