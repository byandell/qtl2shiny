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
pheno_data <- DOread::read_pheno_tbl(analyses_tbl, datapath)
pheno_data <- dplyr::select(pheno_data, 
                            which(names(pheno_data) %in% peaks$pheno))

## Add Closed Reference OTUs
if(dir.exists(file.path(datapath, "otu"))) {
  peaks_otu <- readRDS(file.path(datapath, "otu", "peaks_OTU_CR.rds"))
  peaks_otu <- 
    dplyr::filter(
      peaks_otu,
      output %in% peaks_otu$output)
  peaks <- dplyr::bind_rows(peaks,
                            peaks_otu[, names(peaks)])
  
  analyses_otu <- 
    dplyr::filter(
      readRDS(file.path(datapath, "otu", "analyses_OTU_CR.rds")),
      output %in% peaks_otu$output)
  analyses_tbl <- 
    dplyr::bind_rows(analyses_tbl,
                     analyses_otu[, names(analyses_tbl)])
  
  pheno_otu <- readRDS(file.path(datapath, "otu", "pheno_OTU_CR.rds"))
  pheno_otu <- dplyr::select(
      as.data.frame(pheno_otu), 
      which(colnames(pheno_otu) %in% peaks_otu$pheno))
  tmp <- matrix(NA, nrow(pheno_data), ncol(pheno_otu))
  dimnames(tmp) <- list(rownames(pheno_data), colnames(pheno_otu))
  m <- match(rownames(pheno_otu), rownames(tmp))
  tmp[m,] <- as.matrix(pheno_otu)
  pheno_data <- cbind(pheno_data, 
                      as.data.frame(tmp))
}

# Add model column to analyses_tbl
analyses_tbl <- 
  dplyr::mutate(
    analyses_tbl,
    model = ifelse(pheno_type == "OTU_Bin",
                   "binary",
                   "normal"))
analyses_tbl <- 
  dplyr::select(
    analyses_tbl,
    pheno:pheno_type, model, transf:ncol(analyses_tbl))

##
pheno_type <- c("all", sort(unique(analyses_tbl$pheno_type)))

K <- readRDS(file.path(datapath, "kinship.rds"))

