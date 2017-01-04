## Aim of this setup file is to create (ultimately) and standalone qtl2shiny
## Some things are not in place yet, so it will not work for now.

## Data: <https://github.com/rqtl/qtl2data/>  
## Reference: Recla JM, Robledo RF, Gatti DM, Bult CJ, Churchill GA, Chesler EJ (2014)
## Precise genetic mapping and integrative bioinformatics
## in Diversity Outbred mice reveals Hydin as a novel pain gene.
## Mamm Genome 25:211-222. doi:10.1007/s00335-014-9508-0

## web_address <- file.path("https://raw.githubusercontent.com/rqtl/",
##                          "qtl2data/master/DOex")
## if(!file.exists("DOex.zip")) {
##   download.file(file.path(web_address, "DOex.zip"), do_file, quiet=TRUE)
## }
## if(!file.exists("c2_snpinfo.rds")) {
##   download.file(file.path(web_address, "c2_snpinfo.rds"), snp_file, quiet=TRUE)
## }

DOex <- read_cross2("DOex.zip")
snpinfo <- readRDS("c2_snpinfo.rds")

pmap <- DOex$pmap
covar <- DOex$covar
pheno_data <- DOex$pheno

## Filter peaks and analyses_tbl to "best" analysis.
## Need to create peaks file
peaks <- NULL

## Need to create peaks file
analyses_tbl <- NULL
pheno_type <- "all"

## Calculate genoprobs and kinship
pr <- calc_genoprob(DOex, error_prob=0.002)
apr <- genoprob_to_alleleprob(pr)
apr$map <- DOex$pmap
K <- calc_kinship(apr, "loco")

## Get peaks
pheno_scan <- scan1(apr, DOex$pheno)

## SNP info in code is drawn from files. Need to refactor.
