gene_region <- function(chr_id, wrng, project_info) {
  
  query_genes <- read_query_rds(project_info, "query_genes.rds")
  
  gene_tbl <- query_genes(chr_id, wrng[1], wrng[2])
  qtl2pattern::get_genes(chr_id, wrng[1], wrng[2], gene_tbl)
}