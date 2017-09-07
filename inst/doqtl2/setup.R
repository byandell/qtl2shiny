dirpath <- "~/Documents/Research/attie_alan/DO/data"
datapath <- file.path(dirpath, "DerivedData")

pmap <- readRDS(file.path(datapath, "pmap.rds"))
covar <- readRDS(file.path(datapath, "covar.rds"))

peaks        <- DOread::setup_peaks(datapath) 
analyses_tbl <- DOread::setup_analyses(peaks, datapath)
pheno_data   <- DOread::setup_data(analyses_tbl, peaks, datapath)
pheno_type   <- DOread::setup_type(analyses_tbl)

K <- readRDS(file.path(datapath, "kinship.rds"))

