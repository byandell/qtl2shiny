# README #

Yandell R/qtl2shiny project.

### What is this repository for? ###

* Code to create shiny interface for [R/qtl2](https://kbroman.org/qtl2).
* Version 0.9.18
* See following documents:
    + [~yandell/software/qtl2shiny](http://www.stat.wisc.edu/~yandell/software/qtl2shiny)
    + [R/qtl2shiny User Guide](https://github.com/byandell/qtl2shiny/blob/master/vignettes/UserGuide.Rmd)
    + [R/qtl2shiny Developer Guide](https://github.com/byandell/qtl2shiny/blob/master/vignettes/DeveloperGuide.Rmd)

### What has been done ###

- created [Shiny](http://shiny.rstudio.org) interface for [R/qtl2](https://kbroman.org/qtl2) data
    + handles multiple phenotypes and multiple projects
    + creates plots and scans on the fly
    + uses Shiny modules and dashboard
- organized in package R/qtl2shiny
- depends on related packages
    + [R/qtl2](https://kbroman.org/qtl2)
    + [R/qtl2ggplot2](https://github.com/byandell/qtl2ggplot2)
    + [R/qtl2feather](https://github.com/byandell/qtl2feather)
    + [R/qtl2pattern](https://github.com/byandell/qtl2pattern)
    + [R/DOread](https://github.com/byandell/DOread)
    + [R/CausalMST](https://github.com/byandell/CausalMST)

### What are open development issues ###

Major issues

- speed up access
- save intermediate calculations that are reused
- make sure multiple taxa work smoothly

Minor issues

* shiny
    + reveal plots more (settings on sidebar? tabs?)
    + user save settings for quick replay of shiny
* markdown
    + flexdashboard for Rmd to create dynamic reports
    + link Rmd and shiny documents
* genes
    + GeneRegion and GeneExon for multiple traits giving different results
    + upgrade these to use `plot_genes`
* scans
    + multiple traits now use superset of covariates
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

    install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "RSQLite"))

You will also need the following packages for qtl2shiny dependencies:

    install.packages(c("tidyverse", "RColorBrewer", "feather",
      "shiny", "shinydashboard", "grid", "gridBase", "gdata", "GGally", "Rcpp", "mnormt", "corpcor"))

Then, install R/qtl2 using `devtools::install_github()`.

    library(devtools)
    install_github("rqtl/qtl2")
    install_github("ropensci/plotly")

Once you have installed these, install qtl2shiny and related packages as

    install_github(paste0("byandell/",
      c("DOread","qtl2pattern","qtl2ggplot","qtl2feather","CausalMST","qtl2shiny")))

To install `qtl2shiny` with vignettes (takes a bit longer):

    install_github("byandell/qtl2shiny", build_vignettes=TRUE)
