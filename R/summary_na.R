#' Summary across traits
#'
#' Summary taking care to include NA for all traits
#'
#' @param phe data frame of phenotypes
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{summary_na(phe)}
#' 
#' @export
summary_na <- function(phe) {
  tmpfn <- function(x) {
    out <- summary(x)
    out["NA's"] <- sum(is.na(x))
    out
  }
  t(sapply(phe, tmpfn))
}
