med_test <- function(med_ls, geno_max, phe_mx, kinship, cov_tar,
                     pos_Mbp, data_type, driver_med = NULL) {

  if(is.list(kinship))
    kinship <- kinship[[1]]
  
  CausalMST::mediation_test(phe_mx, med_ls[[1]], med_ls[[2]], geno_max,
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
