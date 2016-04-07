# README #

Yandell contribution to doqtl2 project.

### What is this repository for? ###

* Code to liaise between Broman work on R/qtl2 and Attie Lab.
* Version 0.2.1

### What has been done ###

* built SQLite database for MGI genes
* rebuilt gene_plot routine for association map SNP plot
* Dependencies: readr, dplyr, tidyr, ggplot2, stats, graphics
* organized in package doqtl2

### What are open development issues ###

Major issues

* using covariates from analysis.csv as done for pipeline
* update as [R/qtl2](http://kbroman.org/qtl2/) evolves

Minor issues

* chr 3 HOMA broken again -- no snps
* ggbiplot not working sometimes in plot.snp_pc (batch Rmd mode?)
* summary.snp_pc: include loadings somehow
* put `snp_csq` analysis before `gene_exon` plots
* try snp_csq plot by both Csq and pattern
* rename `qtl2_pipeline` as `qtl2_onepheno`

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
