---
title: "R/doqtl2 migration"
author: "Brian S Yandell, <http://www.stat.wisc.edu/~yandell>"
date: "January 2017"
output: html_document
---

migrated to R/qtl2pattern

file              | function
----------------- | --------
shinyGeneExon.R   | plot_gene_exon
shinyGeneRegion.R | get_mgi_features
shinyPattern.R    | scan_pattern
shinyProbs.R      | get_snpprobs read_probs read_probs36
shinySNPAllele.R  | get_gene_exon_snp get_top_snps_tbl snpprob_collapse
shinyScan1Plot.R  | listof_scan1coefCC
shinySetup.R      | get_pheno
shinyTopFeature.R | merge_feature

function          | file | use
----------------- | -------- | --------
get_gene_exon_snp | gene_exon.R | call get_mgi_features
get_mgi_features  | get_mgi_features.R | extract from SQLite
get_pheno         | get_traits.R | get selected phenotypes
get_snpprobs      | snpinfo.R   | snpprobs for SNPs, InDels, SVs
get_top_snps_tbl  | top_snps_tbl.R | get top SNP info based on LMMs
listof_scan1coefCC | listof_scan1coefCC.R | create list of scan1coefCC objects
merge_feature     | merge_feature.R | merge SNP LOD and other information
plot_gene_exon    | gene_exon.R | plot genes and exons
read_probs        | read_probs.R | read genoprob object for RDS
read_probs36      | read_probs.R | read genoprob object for RDS
scan_pattern      | scan_pattern.R | genome scan by pattern set
snpprob_collapse  | genoprob_to_patternprob.R | collapse from alleles to SNPs

Want `scan_pattern` to look like `plot_snpasso` with `pattern="all"`.

Need to work on geno and exon stuff to meld with `plot_genes`.
