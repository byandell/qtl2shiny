## Various utilities for qtl2shiny

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
#' \dontrun{num_Pheno(pheno, analyses_tbl)}
#' 
#' @export
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
  tot_pheno <- nrow(analyses_tbl %>%
                      distinct(pheno, .keep_all=TRUE))
  ## Put in human-readable format
  num_pheno <- hr_num(num_pheno, 0)
  tot_pheno <- hr_num(tot_pheno, 2)
  paste("Phenotypes:", num_pheno, "of", tot_pheno)
}
#' Make chr_pos from chr, peak and window
#'
#' Create character string of chr:left-right with left=peak-window, right=peak+window
#'
#' @param chr_id chromosome ID
#' @param peak_Mbp position of peak in Mbp
#' @param window_Mbp half-width of window around peak in Mbp
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{make_chr_pos(chr_id, peak_Mbp, window_Mbp)}
#' 
#' @export
make_chr_pos <- function(chr_id, peak_Mbp, window_Mbp) {
  if(is.null(chr_id))
    chr_id <- "?"
  left_Mbp <- right_Mbp <- "?"

  if(is.null(peak_Mbp))
    peak_Mbp <- "?"
  else
    peak_Mbp <- round(peak_Mbp, 2)

  if(is.null(window_Mbp))
    window_Mbp <- "?"
  else {
    window_Mbp <- round(window_Mbp, 2)
    if(peak_Mbp != "?") {
      left_Mbp <- peak_Mbp - window_Mbp
      right_Mbp <- peak_Mbp + window_Mbp
    }
  }
  paste(chr_id, left_Mbp, right_Mbp, sep = "_")
}
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
collapse_covar <- function(analyses_tbl) {
  if(is.null(analyses_tbl))
    return(NULL)
  
  ## Collapse covariates (past winsorize column).
  covar_names <- names(analyses_tbl)[-seq_len(match("winsorize",
                                           names(analyses_tbl)))]
  analyses_tbl %>%
    unite(covar, one_of(covar_names)) %>%
    mutate(covar = sapply(strsplit(covar,"_"),
                          function(x) paste(covar_names[as.logical(x)],
                                            collapse=",")))
}
