#' @importFrom ggplot2 ggtitle
#' 
plot_gene_exons <- function(gene_exon_tbl, top_snps_tbl, gene_name, pheno) {
  if(nrow(top_snps_tbl)) {
    p <- plot(gene_exon_tbl, top_snps_tbl,
              FALSE, genes = gene_name)
    p[[1]] + ggplot2::ggtitle(paste(gene_name, "SNPs for", pheno))
  } else
    plot_null()
}
