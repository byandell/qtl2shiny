med_test <- function(chr_id, pos_Mbp, med_ls,
                     phe_df, cov_tar, aprobs, kinship, map,
                     data_type) {

  cov_med <- med_ls$cov_med
  
  # Get genotype matrix and map at 
  peak_mar <- qtl2geno::find_marker(map, chr_id, pos_Mbp)
  probs_max <- subset(aprobs, chr = chr_id, mar = peak_mar)
  map[[chr_id]] <- map[[chr_id]][peak_mar]
  geno_max <- probs_max[[1]][,,1]
  
  CausalMST::mediate1_test(med_ls, geno_max, phe_df,
                           kinship[[chr_id]], cov_tar, cov_med,
                           "wilc", pos = pos_Mbp,
                           data_type = data_type)
}