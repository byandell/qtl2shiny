query_probs <- 
  DOread::create_probs_query_func_do(
    file.path("qtl2shinyData", "mouse", "AttieDO"))
saveRDS(query_probs, "query_probs.rds")
query_mrna <- 
  DOread::create_mrna_query_func_do(
    file.path("qtl2shinyData", "mouse", "AttieDO"))
saveRDS(query_mrna, "query_mrna.rds")
