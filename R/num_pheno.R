#' Number of phenotypes examined
#'
#' Character string with number of phenotypes of total.
#'
#' @param pheno names of phenotypes
#' @param analyses_tbl table of analyses run
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{num_pheno(pheno, analyses_tbl)}
#' 
#' @export
#' 
#' @importFrom gdata humanReadable
#' @importFrom dplyr distinct
#' 
num_pheno <- function(pheno, analyses_tbl) {
  if(any(c("all","none") %in% pheno))
    return(NULL)
  
  hr_num <- function(x, digits_ct=0) {
    x <- gdata::humanReadable(x, digits=digits_ct, sep="",
                              standard="SI",
                              justify=c("right","right"))
    substring(x, 1, nchar(x) - 1)
  }
  
  num_pheno <- length(unique(pheno))
  tot_pheno <- nrow(dplyr::distinct(analyses_tbl, pheno, .keep_all=TRUE))
  ## Put in human-readable format
  num_pheno <- hr_num(num_pheno, 0)
  tot_pheno <- hr_num(tot_pheno, 2)
  paste("Phenotypes:", num_pheno, "of", tot_pheno)
}
