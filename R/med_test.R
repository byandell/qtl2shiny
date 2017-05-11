med_test <- function(med_ls, geno_max, phe_df, kinship, cov_tar,
                     pos_Mbp, data_type, driver_med = NULL) {

  if(is.list(kinship))
    kinship <- kinship[[1]]
  
  CausalMST::mediate1_test(med_ls, geno_max, phe_df,
                           kinship, cov_tar, med_ls$cov_med,
                           driver_med,
                           test = "wilc", pos = pos_Mbp,
                           data_type = data_type)
}