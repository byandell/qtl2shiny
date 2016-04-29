# README #

Yandell contribution to doqtl2 project.

### What is this repository for? ###

* Code to liaise between Broman work on R/qtl2 and Attie Lab.
* Version 0.4

### What has been done ###

* built SQLite database for MGI genes
* rebuilt gene_plot routine for association map SNP plot
* Rmd pipelines for single pheno and single peaks
* Rmd pipeline for multiple phenotypes (say <10)
* incorporated shiny parameter setting for users
* intermediate calculations saved in rds files
* Dependencies: readr, dplyr, tidyr, ggplot2, stats, graphics
* Suggests: shiny (for params), ggbiplot (for PCA of patterns)
* organized in package doqtl2

### What are open development issues ###

Major issues

* update as [R/qtl2](http://kbroman.org/qtl2/) evolves
* incorporate Ensembl ID as well as Gene Name
* Get Csq and other gene features from SNP database rather than Sanger site
  + got Csq, but need gene name and/or other features in CSQ part of VCF
* threshold calculations
* probs will only grow, but SQLite does not seem easy answer
  + [NetCDF and HDF5](http://www.unidata.ucar.edu/software/netcdf/docs/interoperability_hdf5.html)? (used at Jax; see [RNetCDF](https://cran.r-project.org/web/packages/RNetCDF/index.html))

Minor issues

* PCA of SNP patterns
  + ggbiplot not working sometimes in plot.snp_pc (batch Rmd mode?)
  + summary.snp_pc: include loadings somehow
* revisit naming of routines
* remove AB1NZCPW dependency as this is specific to DO cross
  + have as parameter setting, etc.
  + but depends on what Karl does with genotype classes
* streamline Rmd output to focus on key findings
  + fewer plots
  + better story around figures

Redesign for use by scientists

* now down to one command to get to shiny parameter check
  + qtl2_rmd() generic interface
  + qtl2_onepheno() more customized for onepheno & onepeak run
* connect to Karl's new server and d3 tools

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* [Brian Yandell](http://bitbucket.org/byandell)
* Working with [Karl Broman](http://bitbucket.org/kbroman)
