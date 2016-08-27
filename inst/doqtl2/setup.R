dirpath <- "~/Documents/Research/attie_alan/DO/data"
datapath <- file.path(dirpath, "DerivedData")

pmap <- readRDS(file.path(datapath, "pmap.rds"))
covar <- readRDS(file.path(datapath, "covar.rds"))

## Filter peaks and analyses_tbl to "best" analysis.
peaks <- readRDS(file.path(datapath, "peaks.rds"))
peak_info <- (peaks %>%
  group_by(pheno) %>%
  distinct(output) %>%
  arrange(desc(output)) %>%
  summarize(output = output[1]) %>%
  ungroup)$output
peaks <- peaks %>%
  filter(output %in% peak_info)

analyses_tbl <- readRDS(file.path(datapath, "analyses.rds")) %>%
  filter(output %in% peak_info)

## Want to filter to "best" analysis based -- anal1 or anal2
peak_info <- peaks$pheno
pheno_data <- read_pheno_tbl(analyses_tbl, datapath)
pheno_data <- pheno_data %>%
  select(which(names(pheno_data) %in% peak_info))
rm(peak_info)

pheno_type <- c("all", sort(unique(analyses_tbl$pheno_type)))

K <- readRDS(file.path(datapath, "kinship.rds"))
