filepath <- file.path(datapath, "mouse", "AttieDO")
query_probs <- 
  DOread::create_probs_query_func_do(filepath)
saveRDS(query_probs, "query_probs.rds")
query_mrna <- 
  DOread::create_mrna_query_func_do(filepath)
saveRDS(query_mrna, "query_mrna.rds")
