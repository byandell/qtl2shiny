med_test <- function(chr_id, pos_Mbp, med_ls,
                     phe_df, cov_tar, aprobs, kinship, map) {

  if("expr" %in% names(med_ls)) {
    mediators <- med_ls$expr
    testfn <- CausalMST::mediate1_test
  } else {
    mediators <- med_ls$comediators
    testfn <- CausalMST::comediate1_test
  }
  annotation <- med_ls$annot
  cov_med <- med_ls$cov_med
  
  # Get genotype matrix and map at 
  peak_mar <- qtl2geno::find_marker(map, chr_id, pos_Mbp)
  probs_max <- subset(aprobs, chr = chr_id, mar = peak_mar)
  map[[chr_id]] <- map[[chr_id]][peak_mar]
  geno_max <- probs_max[[1]][,,1]
  
  testfn(geno_max, phe_df, mediators, qtl2scan::fit1,
                                       kinship[[chr_id]], cov_tar, cov_med,
                                       annotation, "wilc", pos = pos_Mbp)
}