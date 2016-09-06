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
#' Plot effects
#'
#' Plot of effects.
#'
#' @param pheno name of phenotype for effect scan
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param eff_obj object of class \code{\link[doqtl2]{listof_scan1coefCC}}
#' @param xlim x limits for plot
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @examples
#' \dontrun{plot_eff(pheno, scan_obj, eff_obj, xlim)}
#' 
#' @export
plot_eff <- function(pheno, scan_obj, eff_obj, xlim) {
  lodcol <- match(pheno, names(eff_obj))
  plot_coefCC(eff_obj[[lodcol]], xlim=xlim)
  max_pos <- max(scan_obj, lodcolumn=lodcol)$pos[1]
  abline(v=max_pos, lty=2, lwd=2, col=lodcol)
  at <- par("usr")[1:2]
  at <- seq(at[1],at[2],length.out=length(CCcolors))
  mtext(names(CCcolors),col=CCcolors,at=at)
  title(pheno)
}
#' Top Pattern Plot
#'
#' Plot of top patterns.
#'
#' @param pheno name of phenotype for effect scan
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param xlim x limits for plot
#' @param fill.null null plot if no data
#' @param group group by one of c("pheno","pattern")
#' @param snp_action character string for plot
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @examples
#' \dontrun{top_pat_plot(pheno, scan_obj, xlim)}
#' 
#' @export
top_pat_plot <- function(pheno, 
                         scan_obj, xlim,
                         fill.null=TRUE, group = "pheno",
                         snp_action = "basic") {
  drop <- max(max(scan_obj$lod) - 1.5, 1.5)
  top_pattern <- topsnp_pattern(scan_obj, pheno, drop)
  if(is.null(top_pattern)) {
    if(fill.null)
      return(plot_null())
    else
      return()
  }
  chr_id <- names(scan_obj$map)[1]
  p <- plot(top_pattern, group=group) + xlim(xlim)
  if(length(pheno) == 1) {
    mytitle <- paste(pheno, "chr", chr_id)
    if(snp_action != "basic")
      mytitle <- paste(mytitle, snp_action)
    p <- p + ggtitle(mytitle)
  }
  p
}
#' @export
plot_null <- function() {
  ggplot(data.frame(x=1,y=1), aes(x,y,label="no data")) +
    geom_text(size=10) + theme_void()
}
#' Top SNP Association Plot Wrapper
#'
#' @param pheno name of phenotype for effect scan
#' @param scan_obj object of class \code{\link[qtl2scan]{scan1}}
#' @param xlim x limits for plot
#' @param snp_action character string for plot
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @examples
#' \dontrun{top_snp_asso(pheno, scan_obj, xlim)}
#' 
#' @export
top_snp_asso <- function(pheno, scan_obj, xlim, snp_action="basic") {
  if(is.null(pheno) | is.null(scan_obj) | is.null(xlim))
    return(print(plot_null()))
  
  phename <- dimnames(scan_obj$lod)[[2]]
  chr_id <- names(scan_obj$map)[1]
  plot_snpasso(subset(scan_obj,
                      lodcolumn=match(pheno, phename)),
               show_all_snps=FALSE, drop.hilit=1.5,
               xlim=xlim)
  mytitle <- paste(pheno, "chr", chr_id)
  if(snp_action != "basic")
    mytitle <- paste(mytitle, snp_action)
  title(mytitle)
}
