#' Plot across traits
#'
#' Plot of density and pairwise scatterplots.
#'
#' @param phe data frame of phenotypes
#' @param cov matrix of covariates
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @examples
#' \dontrun{plot_sex(phe)}
#' 
#' @export
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes_string geom_density geom_rug
#' @importFrom GGally ggscatmat
plot_sex <- function(phe, cov) {
  phename <- names(phe)
  if(length(phename) > 10) {
    cat(file=stderr(), "\nOnly first 10 phenotypes used\n")
    phename <- phename[seq_len(10)]
    phe <- phe[,phename]
  }
  if("sex" %in% dimnames(cov)[[2]]) {
    ## Assume sex in covar. Ignore actual covariates for analyses.
    insex <- data.frame(phe, cov)
    
    if(length(phename) == 1) {
      ggplot2::ggplot(insex, 
                      ggplot2::aes_string(phename, col="sex")) +
        ggplot2::geom_density(na.rm = TRUE) + 
        ggplot2::geom_rug()
    } else {
      any.na <- apply(insex, 1, function(x) any(is.na(x)))
      GGally::ggscatmat(insex[!any.na,], seq_along(phe), color="sex")
    }
  } else {
    if(length(phename) == 1) {
      ggplot2::ggplot(phe, 
                      ggplot2::aes_string(phename)) +
        ggplot2::geom_density(na.rm = TRUE) + 
        ggplot2::geom_rug()
    } else {
      any.na <- apply(phe, 1, function(x) any(is.na(x)))
      GGally::ggscatmat(phe[!any.na,], seq_along(phe))
    }
  }
}
