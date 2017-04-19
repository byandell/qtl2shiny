med_test <- function(chr_id, pos_Mbp, expr_ls,
                     phe_df, cov_tar, aprobs, kinship, map) {

  # Get expression mRMNA measurements.
  annot.mrna <- expr_ls$annot
  expr.mrna <- expr_ls$expr
  # Covariate matrix covar is global, but reget it here to be sure.
  cov_med <- expr_ls$cov_med
  
  # Get genotype matrix and map at 
  peak_mar <- qtl2geno::find_marker(map, chr_id, pos_Mbp)
  probs_max <- subset(aprobs, chr = chr_id, mar = peak_mar)
  map[[chr_id]] <- map[[chr_id]][peak_mar]
  geno_max <- probs_max[[1]][,,1]
  
  CausalMST::mediate1_test(geno_max, phe_df, expr.mrna, qtl2scan::fit1,
                                       kinship[[chr_id]], cov_tar, cov_med,
                                       annot.mrna, "wilc", pos = pos_Mbp)
}