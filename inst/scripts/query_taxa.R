query_genes <-
  qtl2::create_gene_query_func(
    file.path("qtl2shinyData", "mouse", "mouse_genes.sqlite"))
saveRDS(query_genes, "query_genes.rds")
query_variants <-
  qtl2::create_variant_query_func(
    file.path("qtl2shinyData", "mouse", "cc_variants.sqlite"))
saveRDS(query_variants, "query_variants.rds")
