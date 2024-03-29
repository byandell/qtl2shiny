% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinySNPPattern.R
\name{shinySNPPattern}
\alias{shinySNPPattern}
\title{Shiny top SNP analysis and plot module}
\usage{
shinySNPPattern(
  id,
  snp_par,
  chr_pos,
  pheno_names,
  snp_scan_obj,
  snpinfo,
  top_snps_tbl,
  gene_exon_tbl,
  allele_info,
  snp_action = shiny::reactive({
     "basic"
 })
)
}
\arguments{
\item{id}{identifier for shiny reactive}

\item{snp_par, chr_pos, pheno_names, snp_scan_obj, snpinfo, top_snps_tbl, gene_exon_tbl, allele_info, snp_action}{reactive arguments}
}
\description{
Shiny module for top SNP analysis and plots, with interfaces \code{shinySNPPatternInput}, \code{shinySNPPatternUI} and  \code{shinySNPPatternOutput}.
}
\author{
Brian S Yandell, \email{brian.yandell@wisc.edu}
}
\keyword{utilities}
