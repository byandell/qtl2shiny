# README #

Yandell contribution to doqtl2 project.

### What is this repository for? ###

* Code to liaise between Broman work on R/qtl2 and Attie Lab.
* Version 0.4.12

### What has been done ###

* built SQLite database for MGI genes
* rebuilt gene_plot routine for association map SNP plot
* incorporated Ensembl ID
  + still need Gene Name
* Shiny interface for multiple phenotypes with modules and dashboard
  + invoke shiny dashboard as doqtl2_app()
  + download RDS data
  + create plots and scans on the fly
* markdown scripts
  + Rmd pipelines for single pheno and single peaks
  + Rmd pipeline for multiple phenotypes (say <10)
  + incorporated shiny parameter setting into rmd files
  + intermediate calculations saved in rds files
  + now down to one command to get to shiny parameter check
    + qtl2_rmd() generic interface
    + qtl2_onepheno() more customized for onepheno & onepeak run
* Dependencies: readr, dplyr, tidyr, ggplot2, stats, graphics
* Suggests: shiny (for params), ggbiplot (for PCA of patterns)
* organized in package R/doqtl2

### What are open development issues ###

Major issues

* update as [R/qtl2](http://kbroman.org/qtl2/) evolves
* threshold calculations in various places
* derived and intermediate data storage
  + using SQLite for now
  + [NetCDF and HDF5](http://www.unidata.ucar.edu/software/netcdf/docs/interoperability_hdf5.html)? (used at Jax; see [RNetCDF](https://cran.r-project.org/web/packages/RNetCDF/index.html))
  + [PostgreSQL](https://www.postgresql.org) (available, used by Ed Buckler for [TASSLE](http://www.maizegenetics.net))
  + [Arvados](https://arvados.org) (used for redesign of [GeneNetwork](http://www.genenetwork.org/))
  + for now use save/readRDS for individual chromosomes
* shiny SNP pattern plot across traits by pattern
* SNPs as covariates as function and shiny module (see inst/tests/best_snp.R)
* connect to Karl's new server and d3 tools
* dynamic plots with [d3panels](http://kbroman.org/d3panels/) and [qtlcharts](http://kbroman.org/qtlcharts/)

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
* shiny
  + reveal plots more (settings on sidebar? tabs?)
  + user save settings for quick replay of shiny
* snpinfo varying names for `pos_Mbp` as `pos` and `snp_id` as `snp`
  + `R/exon_snp.R`, `R/check_interval.R`, `R/snpinfo.R`, others?
* markdown
  + flexdashboard for Rmd to create dynamic reports
  + link Rmd and shiny

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* [Brian Yandell](http://bitbucket.org/byandell)
* [Karl Broman](http://bitbucket.org/kbroman)
  + [R/qtl2](http://kbroman.org/qtl2/)
