med_test <- function(med_ls, geno_max, phe_mx, kinship, cov_tar,
                     pos_Mbp, data_type, driver_med = NULL) {

  if(is.list(kinship))
    kinship <- kinship[[1]]
  
  # Somewhere in this, cov_med ends up as vector rather than column matrix.
  # Seems mediate1_test wants cov_med as dataframe, then later convert?
  # see CausalMST:::cmst_pheno is where problem lies.
  CausalMST::mediate1_test(phe_mx, med_ls[[1]], med_ls[[2]], geno_max,
                           cov_tar, med_ls$covar, kinship,
                           driver_med,
                           test = "wilc", pos = pos_Mbp,
                           data_type = data_type)
}

med_scat <- function(med_ls, geno_max, phe_mx, kinship, cov_tar, sdps, 
                     pattern, med_name, medID, haplos) {
  
  if(is.list(kinship))
    kinship <- kinship[[1]]
  
  # Could use qtl2pattern::pull_mediator here.
  
  sdp <- sdps[qtl2pattern::sdp_to_pattern(sdps, haplos) == pattern]
  id <- med_ls[[2]]$id[med_ls[[2]][[medID]] == med_name]
  if(length(id) != 1)
    return(NULL)
  
  cov_tar <- qtl2pattern::covar_df_mx(cov_tar)
  cov_med <- qtl2pattern::covar_df_mx(med_ls$covar)
  
  CausalMST:::med_scatter(geno_max, phe_mx, med_ls[[1]][, id, drop = FALSE],
                          kinship, cov_tar, cov_med,
                          qtl2::fit1,
                          sdp = sdp, allele = TRUE)
}
