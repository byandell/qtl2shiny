allele_scan <- function(phe_mx, cov_df, probs_obj, K_chr,
                    patterns, scan_pat, blups) {
  addcovar <- covar_df_mx(cov_df)
  qtl2pattern::allele1(probs_obj$probs, phe_mx, addcovar, 
                       probs_obj$map, K_chr,
                       patterns = patterns,
                       scan_pat = scan_pat,
                       blups = blups)
}
