query_genes <- 
  qtl2db::create_gene_query_func(
    file.path("DerivedData", "mouse", "mouse_genes.sqlite"))
saveRDS(query_genes, "query_genes.rds")
query_variants <- 
  qtl2db::create_variant_query_func(
    file.path("DerivedData", "mouse", "cc_variants.sqlite"))
saveRDS(query_variants, "query_variants.rds")
