## ----echo=FALSE----------------------------------------------------------
## Attach packages
suppressPackageStartupMessages({
  library(qtl2geno)
  library(qtl2scan)
  library(qtl2plot)
  library(doqtl2)
})

## Directories -- required.
datapath <- params$datapath
if(is.null(datapath))
  stop("must supply datapath")
resultpath <- params$resultpath
if(is.null(resultpath))
  stop("must supply resultpath")

## Phenotype selection and chromosome ID -- required.
phename_select <- params$phename_select
if(is.null(phename_select))
  stop("must supply phename_select")
if(phename_select=="")
  stop("must supply phename_select")

chr_id <- params$chr_id
if(is.null(chr_id))
  stop("must supply chr_id")
chr_id <- as.character(chr_id)
if(chr_id=="")
  stop("must supply chr_id")

## All other arguments are optional
set_null <- function(x, params) {
  x <- params[[x]]
  if(all(x==""))
    x <- NULL
  x
}
phename_peak <- as.numeric(params$phename_peak)
if(!length(phename_peak)) {
  phename_peak <- 0
}

phename_window <- params$phename_window

## Phenotypes to drop. Convert from space delimited to vector.
phename_drop <- set_null("phename_drop", params)
if(!is.null(phename_drop))
  phename_drop <- strsplit(phename_drop, " ")[[1]]

analysis_type <- set_null("analysis_type", params)
proband <- set_null("proband", params)
best_snp <- set_null("best_snp", params)

if(is.null(params$shorten_char)) {
  shorten_char <- 0
} else {
  shorten_char <- round(as.numeric(params$shorten_char))
}

## Redo intermediate calculations (including reload of probs).
## Redo automatically happens if resultpath is new.
redo <- params$redo
if(is.null(redo))
  redo <- 0 ## 0=no redo; 1=reload probs; 2=recalc everything
if(is.character(redo))
  redo <- match(redo, c("no", "some", "all"), nomatch=1) - 1

## ------------------------------------------------------------------------
## Use analyses.csv table to determine phename and pheno_set
analysis_tbl <- read_csv(file.path(datapath, "analyses.csv"))
if(phename_select != "") {
  analysis_tbl <- analysis_tbl %>%
    filter(grepl(phename_select, pheno))
  if(length(phename_drop)) for(type in phename_drop) {
    analysis_tbl <- analysis_tbl %>%
      filter(!grepl(type, pheno))
  }
}

if(phename_peak > 0) {
  load(file.path(datapath, "peaks.RData"))
  peaks <- peaks %>%
    filter(chr==chr_id,
           pos >= phename_peak - phename_window,
           pos <= phename_peak + phename_window)
  analysis_tbl <- analysis_tbl[match(peaks$pheno, analysis_tbl$pheno),]
}
phename <- unique(analysis_tbl$pheno)
shortname <- shorten_phename(phename, shorten_char)
pheno_set <- paste0("pheno_", unique(analysis_tbl$pheno_group))

## Get traits: phenotpyes phe and covariates phe
covar <- get_traits(phename, pheno_set, TRUE, analysis_type, TRUE, datapath)
phename_output <- covar$output
phe <- covar$phe
covar <- covar$covar

## ------------------------------------------------------------------------
phename

## ------------------------------------------------------------------------
dimnames(covar)[[2]]

