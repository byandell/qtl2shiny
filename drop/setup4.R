dirpath <- "~/Documents/Research/attie_alan/DO/data"
datapath <- file.path(dirpath, "DerivedData/wave4")

pmap <- readRDS(file.path(datapath, "pmap.rds"))
peaks <- readRDS(file.path(datapath, "peaks.rds"))
covar <- readRDS(file.path(datapath, "covar.rds"))

peak_info <- peaks$output
analyses_tbl <- readRDS(file.path(datapath, "analyses.rds"))
peak_info <- peaks$pheno
pheno_data <- read_pheno_tbl(analyses_tbl, datapath)
pheno_data <- pheno_data %>%
  select(which(names(pheno_data) %in% peak_info))
rm(peak_info)

pheno_type <- c("all", sort(unique(analyses_tbl$pheno_type)))

K <- readRDS(file.path(datapath, "kinship.rds"))
