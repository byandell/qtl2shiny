#' Get covariates
#'
#' Get covariates using \code{analyses_tbl} filtered by \code{phename_output}.
#'
#' @param covar matrix of all covariates
#' @param analyses_tbl table of analyses setups
#'
#' @return matrix of covariates
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"
#' 
#' # Read DOex example covar from 'qtl2data'
#' covar <- read.csv(file.path(dirpath, "DOex_covar.csv"))
#' 
#' # Select Sex and Cohort columns of covariates
#' analyses_tbl <- data.frame(Sex = TRUE, Cohort = TRUE)
#' # Pull covariates
#' covar <- get_covar(covar, analyses_tbl)
#' tidyr::pivot_wider(
#'   dplyr::count(covar, Sex, Cohort),
#'   names_from = "Sex", values_from = "n")
#'
get_covar <- function(covar, analyses_tbl) {
  ## Get covariate matrix covar.
  logical.col <- sapply(analyses_tbl, function(x) any(as.logical(x)))
  covar[, match(names(analyses_tbl)[logical.col],
                dimnames(covar)[[2]], nomatch=0),
        drop=FALSE]
}
