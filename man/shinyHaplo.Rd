% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinyHaplo.R
\name{shinyHaplo}
\alias{shinyHaplo}
\title{Shiny haplotype analysis}
\usage{
shinyHaplo(
  id,
  win_par,
  pmap_obj,
  phe_mx,
  cov_df,
  K_chr,
  analyses_df,
  covar,
  analyses_tbl,
  peaks,
  project_info,
  allele_info
)
}
\arguments{
\item{id}{identifier for shiny reactive}

\item{win_par, pmap_obj, phe_mx, cov_df, K_chr, analyses_df, covar, analyses_tbl, peaks, project_info, allele_info}{reactive arguments}
}
\value{
No return value; called for side effects.
}
\description{
Shiny module for analysis based on haplotype alleles, with interface \code{shinyHaploUI}.
}
\author{
Brian S Yandell, \email{brian.yandell@wisc.edu}
}
\keyword{utilities}
