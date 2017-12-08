allele_scan <- function(phe_df, cov_df, probs_obj, K_chr,
                    patterns, scan_pat, blups) {
  addcovar <- covar_df_mx(cov_df)
  qtl2pattern::allele1(phe_df, addcovar, probs_obj$probs, 
                       probs_obj$map, K_chr,
                       patterns = patterns,
                       scan_pat = scan_pat,
                       blups = blups)
}