% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinyPeaks.R
\name{shinyPeaks}
\alias{shinyPeaks}
\title{Shiny peaks selection}
\usage{
shinyPeaks(id, set_par, pheno_type, peaks_tbl, pmap_obj, project_info)
}
\arguments{
\item{id}{identifier for shiny reactive}

\item{set_par, pheno_type, peaks_tbl, pmap_obj, project_info}{reactive arguments}
}
\value{
list of inputs and scan summary
}
\description{
Shiny module for peaks selection, with interfaces \code{shinyPeaksInput}, \code{shinyPeaksUI} and  \code{shinyPeaksOutput}.
}
\author{
Brian S Yandell, \email{brian.yandell@wisc.edu}
}
\keyword{utilities}
