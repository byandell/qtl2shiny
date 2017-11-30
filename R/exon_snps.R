exon_snps <- function(tops, project_info) {
  
  query_genes <- read_query_rds(project_info, "query_genes.rds")
  
  CCSanger::get_gene_exon_snp(tops)
}