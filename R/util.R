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
#' \dontrun{num_pheno(pheno, analyses_tbl)}
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
#' @param range left and right window positions in Mbp
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{make_chr_pos(chr_id, peak_Mbp, window_Mbp)}
#' 
#' @export
make_chr_pos <- function(chr_id=NULL, peak_Mbp=NULL, window_Mbp=NULL,
                         range = c(peak_Mbp - window_Mbp,
                                    peak_Mbp + window_Mbp)) {
  if(is.null(chr_id))
    chr_id <- "?"
  if(length(range) < 2 )
    range <- rep("?", 2)
  paste(chr_id, range[1], range[2], sep = "_")
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
                                            collapse=","))) %>%
    select(pheno,covar,transf,offset,winsorize)
}
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
plot_sex <- function(phe, cov) {
  phename <- names(phe)
  if(length(phename) > 10) {
    cat(file=stderr(), "\nOnly first 10 phenotypes used\n")
    phename <- phename[seq_len(10)]
    phe <- phe[,phename]
  }
  if("sex" %in% dimnames(cov)[[2]]) {
    ## Assume sex in covar. Ignore actual covariates for analyses.
    insex <- data.frame(phe,cov) %>%
      mutate(sex=c("female","male")[1+sex])
    
    if(length(phename) == 1) {
      ggplot(insex, aes_string(phename, col="sex")) +
        geom_density(na.rm = TRUE) + geom_rug()
    } else {
      any.na <- apply(insex, 1, function(x) any(is.na(x)))
      ggscatmat(insex[!any.na,], seq_along(phe), color="sex")
    }
  } else {
    if(length(phename) == 1) {
      ggplot(phe, aes_string(phename)) +
        geom_density(na.rm = TRUE) + geom_rug()
    } else {
      any.na <- apply(phe, 1, function(x) any(is.na(x)))
      ggscatmat(phe[!any.na,], seq_along(phe))
    }
  }
}
#' @export
plot_null <- function() {
  ggplot(data.frame(x=1,y=1), aes(x,y,label="no data")) +
    geom_text(size=10) + theme_void()
}
