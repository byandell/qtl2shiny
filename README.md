# README #

Yandell contribution to doqtl2 project.

### What is this repository for? ###

* Code to liaise between Broman work on R/qtl2 and Attie Lab.
* Version 0.3

### What has been done ###

* built SQLite database for MGI genes
* rebuilt gene_plot routine for association map SNP plot
* Dependencies: readr, dplyr, tidyr, ggplot2, stats, graphics
* organized in package doqtl2

### What are open development issues ###

Major issues

* update as [R/qtl2](http://kbroman.org/qtl2/) evolves
* incorporate Ensembl ID as well as Gene Name
* Get Csq from SNP database rather than Sanger site

Minor issues

* ggbiplot not working sometimes in plot.snp_pc (batch Rmd mode?)
* summary.snp_pc: include loadings somehow
* use output rather than phename for folder
* organize DerivedData in its own folder

Redesign for use by scientists

* save `analysis_*.Rmd` as development version
* make `analysis_*.Rmd` production version not echoing commands
* could have echo set globally to turn on/off?

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Brian Yandell
* Working with [Karl Broman](http://bitbucket.org/kbroman)
