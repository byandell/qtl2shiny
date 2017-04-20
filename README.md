# README #

Yandell R/qtl2shiny project.

### What is this repository for? ###

* Code to create shiny interface for [R/doqtl2](https://bitbucket.org/byandell/doqtl2).
* Version 0.3.2

### What has been done ###

* Shiny interface for multiple phenotypes with modules and dashboard
  + invoke shiny dashboard as doqtl2_app()
  + download RDS data
  + create plots and scans on the fly
* organized in package R/qtl2shiny
* depends on R/doqtl2

### What are open development issues ###

Major issues

* Navigation across panels improved
* speed up access
* save intermediate calculations that are reused

Minor issues

* shiny
    + reveal plots more (settings on sidebar? tabs?)
    + user save settings for quick replay of shiny
* markdown
    + flexdashboard for Rmd to create dynamic reports
    + link Rmd and shiny
* genes
    + GeneRegion and GeneExon for multiple traits giving different results
    + make sure GeneExon has null plot if no genes
    + think about many traits, BLUPs
* scans
    + multiple traits now use superset of covariates
    + needs redesign of Scan1Plot and use of cov_mx() from DOread::get_covar
    + fixed for Scan1Plot; need to do same for Scan1SNP and possibly elsewhere

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* [Brian Yandell](http://bitbucket.org/byandell)
* [Karl Broman](http://bitbucket.org/kbroman)
  + [R/qtl2](http://kbroman.org/qtl2/)

### Installation

R/qtl2 is early in development and so is not yet available on
[CRAN](http://cran.r-project.org).

You can install R/qtl2 from [GitHub](https://github.com/rqtl).

You first need to install the
[devtools](https://github.com/hadley/devtools) package, plus a set of
package dependencies: [yaml](https://cran.r-project.org/package=yaml),
[jsonlite](https://cran.r-project.org/package=jsonlite),
[data.table](https://cran.r-project.org/package=data.table),
and [RcppEigen](https://github.com/RcppCore/RcppEigen).
(Additional, secondary dependencies will also be installed)

    install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen"))

You will also need the following packages for qtl2shiny dependencies:

    install.packages(c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "stringr",
      "shiny", "shinydashboard", "grid", "gdata", "GGally", "Rcpp", "mnormt", "corpcor"))

Then, install R/qtl2 using `devtools::install_github()`.

    library(devtools)
    install_github(paste0("rqtl/qtl2", c("geno", "scan")))
    install_github("ropensci/plotly")

Once you have installed these, install qtl2shiny as

    install_github(paste0("byandell/",
      c("DOread","CCSanger","qtl2pattern","qtl2ggplot","qtl2feather","CausalMST","qtl2shiny")))

