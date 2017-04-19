med_test <- function(chr_id, pos_Mbp, window_Mbp,
                     phe_df, cov_tar, aprobs, kinship, 
                     map, datapath) {
  indID <- rownames(phe_df)
  scan_window <- pos_Mbp + c(-1,1) * window_Mbp
  
  # Get expression mRMNA measurements.
  expr.mrna <- DOread::read_mrna(indID, chr_id, 
                                 scan_window[1], scan_window[2], 
                                 datapath)
  annot.mrna <- expr.mrna$annot
  expr.mrna <- expr.mrna$expr
  
  # Covariate matrix covar is global, but reget it here to be sure.
  cov_med <- readRDS(file.path(datapath, "covar.rds"))[, c("sex", paste0("DOwave", 2:4))]
  
  # Get genotype matrix and map at 
  peak_mar <- qtl2geno::find_marker(pmap, chr_id, pos_Mbp)
  probs_max <- subset(aprobs, chr = chr_id, mar = peak_mar)
  map[[chr_id]] <- map[[chr_id]][peak_mar]
  geno_max <- probs_max[[1]][,,1]
  
  CausalMST::mediate1_test(geno_max, phe_df, expr.mrna, qtl2scan::fit1,
                                       kinship[[chr_id]], cov_tar, cov_med,
                                       annot.mrna, "wilc", pos = pos_Mbp)
}