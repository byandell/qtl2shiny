% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinySNPPlot.R
\name{shinySNPPlot}
\alias{shinySNPPlot}
\title{Shiny scan1 SNP analysis and plot module}
\usage{
shinySNPPlot(
  id,
  snp_par,
  chr_pos,
  pheno_names,
  snp_scan_obj,
  snpinfo,
  snp_action = shiny::reactive({
     "basic"
 })
)
}
\arguments{
\item{id}{identifier for shiny reactive}

\item{snp_par, chr_pos, pheno_names, snp_scan_obj, snpinfo, snp_action}{reactive arguments}
}
\description{
Shiny module for scan1 analysis and plots, with interfaces \code{shinySNPPlotUI} and  \code{shinySNPPlotOutput}.
}
\author{
Brian S Yandell, \email{brian.yandell@wisc.edu}
}
\keyword{utilities}
