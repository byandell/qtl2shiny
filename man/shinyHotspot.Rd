% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinyHotspot.R
\name{shinyHotspot}
\alias{shinyHotspot}
\title{Shiny hotspot views}
\usage{
shinyHotspot(id, set_par, pheno_type, peaks_tbl, pmap_obj, project_info)
}
\arguments{
\item{id}{identifier for shiny reactive}

\item{set_par, pheno_type, peaks_tbl, pmap_obj, project_info}{reactive arguments}
}
\value{
list of inputs and scan summary

No return value; called for side effects.
}
\description{
Shiny module to view hotspots for peak selection, with interfaces \code{shinyHotspotInput} and  \code{shinyHotspotOutput}.
}
\author{
Brian S Yandell, \email{brian.yandell@wisc.edu}
}
\keyword{utilities}
