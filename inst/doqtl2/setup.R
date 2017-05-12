dirpath <- "~/Documents/Research/attie_alan/DO/data"
datapath <- file.path(dirpath, "DerivedData")

pmap <- readRDS(file.path(datapath, "pmap.rds"))
covar <- readRDS(file.path(datapath, "covar.rds"))

## Filter peaks and analyses_tbl to "best" analysis.
peaks <- readRDS(file.path(datapath, "peaks.rds"))

peak_info <- dplyr::ungroup(
  dplyr::summarize(
    dplyr::arrange(
      dplyr::distinct(
        dplyr::group_by(peaks, pheno), 
        output), 
      dplyr::desc(output)), 
    output = output[1]))$output

peaks <- dplyr::filter(peaks, output %in% peak_info)
peaks <- dplyr::mutate(peaks, chr = factor(chr, c(1:19,"X")))

analyses_tbl <- dplyr::filter(readRDS(file.path(datapath, "analyses.rds")), 
                              output %in% peak_info)

## Want to filter to "best" analysis based -- anal1 or anal2
peak_info <- peaks$pheno
pheno_data <- DOread::read_pheno_tbl(analyses_tbl, datapath)
pheno_data <- dplyr::select(pheno_data, 
                            which(names(pheno_data) %in% peak_info))
rm(peak_info)

## Add Closed Reference OTUs
peaks_otu <- readRDS(file.path(datapath, "otu", "peaks_OTU_CR.rds"))
peaks <- dplyr::bind_rows(peaks, 
                          peaks_otu[, names(peaks)])

analyses_otu <- readRDS(file.path(datapath, "otu", "analyses_OTU_CR.rds"))
analyses_tbl <- dplyr::bind_rows(analyses_tbl,
                                 analyses_otu[, names(analyses_tbl)])

pheno_otu <- readRDS(file.path(datapath, "otu", "pheno_OTU_CR.rds"))
tmp <- matrix(NA, nrow(pheno_data), ncol(pheno_otu))
dimnames(tmp) <- list(rownames(pheno_data), colnames(pheno_otu))
m <- match(rownames(pheno_otu), rownames(tmp))
tmp[m,] <- pheno_otu
pheno_data <- cbind(pheno_data, tmp)

##
pheno_type <- c("all", sort(unique(analyses_tbl$pheno_type)))

K <- readRDS(file.path(datapath, "kinship.rds"))
