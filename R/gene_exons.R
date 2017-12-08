gene_exons <- function(tops, project_info) {
  
  query_genes <- read_query_rds(project_info, "query_genes.rds")
  
  qtl2pattern::get_gene_exon_snp(tops)
}