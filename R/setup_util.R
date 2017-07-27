#' @export
setup_peaks <- function(datapath) {
  
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

  if(dir.exists(file.path(datapath, "otu"))) {
    ## Add Closed Reference OTUs, OTU Modules, and Bile Acids
    read_otu <- function(peaks, datapath, filename) {
      if(file.exists(file.path(datapath, "otu", filename))) {
        peaks_otu <- readRDS(file.path(datapath, "otu", filename))
        peaks_otu <- 
          dplyr::filter(
            peaks_otu,
            output %in% peaks_otu$output)
        peaks <- dplyr::bind_rows(peaks,
                                  peaks_otu[, names(peaks)])
      }
      peaks
    }
    peaks <- read_otu(peaks, datapath, "peaks_OTU_CR.rds")
    peaks <- read_otu(peaks, datapath, "peaks_OTU_Module.rds")
    ## Replace bile acid with new BileAcid
    if(file.exists(file.path(datapath, "otu", "peaks_BileAcid.rds")))
      peaks <- dplyr::filter(peaks, !(pheno_type == "bile acid"))
    peaks <- read_otu(peaks, datapath, "peaks_BileAcid.rds")
  }
  peaks
}

#' @export
setup_analyses <- function(peaks, datapath) {
  analyses_tbl <- dplyr::filter(readRDS(file.path(datapath, "analyses.rds")), 
                                output %in% peaks$output)
  if(dir.exists(file.path(datapath, "otu"))) {
    ## Add Closed Reference OTUs, OTU Modules, and Bile Acids
    read_otua <- function(analyses_tbl, peaks, datapath, filename) {
      if(file.exists(file.path(datapath, "otu", filename))) {
        analyses_otu <- 
          dplyr::filter(
            readRDS(file.path(datapath, "otu", filename)),
            output %in% peaks$output)
        analyses_tbl <- 
          dplyr::bind_rows(analyses_tbl,
                           analyses_otu[, names(analyses_tbl)])
      }
      analyses_tbl
    }
    analyses_tbl <- read_otua(analyses_tbl, peaks, datapath, "analyses_OTU_CR.rds")
    analyses_tbl <- read_otua(analyses_tbl, peaks, datapath, "analyses_OTU_Module.rds")
    ## Replace bile acid with new BileAcid
    analyses_tbl <- dplyr::filter(analyses_tbl, !(pheno_type == "bile acid"))
    analyses_tbl <- read_otua(analyses_tbl, peaks, datapath, "analyses_BileAcid.rds")
    
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
  }
  analyses_tbl
}
  
#' @export
setup_data <- function(analyses_tbl, peaks, datapath) {
  ## Want to filter to "best" analysis based -- anal1 or anal2
  analyses_std <- dplyr::filter(analyses_tbl,
                                pheno_group %in% c("clin","gutMB","otu","otufam"))
  pheno_data <- DOread::read_pheno_tbl(analyses_std, datapath)
  pheno_data <- dplyr::select(pheno_data, 
                              which(names(pheno_data) %in% peaks$pheno))
  
  ## Add Closed Reference OTUs
  if(dir.exists(file.path(datapath, "otu"))) {
    peaks_otu <- dplyr::filter(peaks,
                               !(pheno_group %in% c("clin","gutMB","otu","otufam")))
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
  pheno_data
}
  
#' @export
setup_type <- function(analyses_tbl) {
  c("all", sort(unique(analyses_tbl$pheno_type)))
}

