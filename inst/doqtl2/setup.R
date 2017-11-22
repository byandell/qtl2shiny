dirpath <- "~/Documents/Research/attie_alan/DO/data"
datapath <- file.path(dirpath, "DerivedData")
v2path <- "~/Documents/Research/attie_alan/DO/AttieDOv2"

pmap <- readRDS(file.path(datapath, "pmap.rds"))

filtered <- file.path(datapath, "filtered")
covar        <- readRDS(file.path(filtered, "covar.rds"))
peaks        <- readRDS(file.path(filtered, "peaks.rds"))
analyses_tbl <- readRDS(file.path(filtered, "analyses.rds"))
pheno_data   <- readRDS(file.path(filtered, "pheno.rds"))
pheno_type   <- readRDS(file.path(filtered, "pheno_type.rds"))

K <- readRDS(file.path(datapath, "kinship.rds"))

query_genes <- 
  qtl2db::create_gene_query_func(
    file.path(v2path, "qtl2db", "mouse_genes.sqlite"))
query_variants <- 
  qtl2db::create_variant_query_func(
    file.path(v2path, "qtl2db", "cc_variants.sqlite"))

