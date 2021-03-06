% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hotspot.R
\name{hotspot}
\alias{hotspot}
\title{Hotspots for phenotypes}
\usage{
hotspot(map, peaks, peak_window = 1, minLOD = 5.5)
}
\arguments{
\item{map}{list of genetic maps}

\item{peaks}{data frame of peak information}

\item{peak_window}{half-width of peak window in Mbp}

\item{minLOD}{minimum LOD to include in count}
}
\value{
object of class hotspot as list of \code{\link[qtl2]{scan1}} and \code{map} objects.
}
\description{
Count hotspots by pheno_group and pheno_type.
}
\examples{
dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"

# Read DOex example cross from 'qtl2data'
DOex <- qtl2::read_cross2(file.path(dirpath, "DOex.zip"))
DOex <- subset(DOex, chr = "2")

# Calculate genotype and allele probabilities
pr <- qtl2::calc_genoprob(DOex, error_prob=0.002)

# Summary of coefficients at scan peak
scan_pr <- qtl2::scan1(pr, DOex$pheno)
peaks <- summary(scan_pr, DOex$pmap)

hotspot(DOex$pmap, peaks)

# Select Sex and Cohort columns of covariates
analyses_tbl <- data.frame(pheno = "OF_immobile_pct", Sex = TRUE, Cohort = TRUE)

# Get hotspot (only one phenotype here).
out <- hotspot(DOex$pmap, peaks)
summary(out)

}
\author{
Brian S Yandell, \email{brian.yandell@wisc.edu}
}
\keyword{utilities}
