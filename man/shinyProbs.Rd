% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinyProbs.R
\name{shinyProbs}
\alias{shinyProbs}
\alias{shinyPairProbs}
\alias{shinySNPProbs}
\title{Shiny Probability module}
\usage{
shinyProbs(id, win_par, project_info)

shinyPairProbs(id, win_par, project_info)

shinySNPProbs(id, win_par, pheno_names, project_info)
}
\arguments{
\item{id}{identifier for shiny reactive}

\item{win_par, pheno_names, project_info}{reactive arguments}
}
\value{
Object of class \code{probs}.
}
\description{
Shiny genotype probability access.
}
\author{
Brian S Yandell, \email{brian.yandell@wisc.edu}
}
\keyword{utilities}
