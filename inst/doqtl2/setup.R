projects <- read.csv("projects.csv", stringsAsFactors = FALSE)
project <- projects$project[1]
datapath <- file.path(projects$directory[1], project)

# Need to think about how to hand off from one project to next.
# Not so easy, as must reside within app rather than global...

pmap         <- readRDS(file.path(datapath, "pmap.rds"))
K            <- readRDS(file.path(datapath, "kinship.rds"))

filtered <- file.path(datapath, "filtered")
covar        <- readRDS(file.path(filtered, "covar.rds"))
peaks        <- readRDS(file.path(filtered, "peaks.rds"))
analyses_tbl <- readRDS(file.path(filtered, "analyses.rds"))
pheno_data   <- readRDS(file.path(filtered, "pheno.rds"))
pheno_type   <- readRDS(file.path(filtered, "pheno_type.rds"))
hotspots     <- readRDS(file.path(filtered, "hotspot.rds"))

# Set up queries
query_genes <- 
  qtl2db::create_gene_query_func(
    file.path(datapath, "qtl2db", "mouse_genes.sqlite"))
query_variants <- 
  qtl2db::create_variant_query_func(
    file.path(datapath, "qtl2db", "cc_variants.sqlite"))
query_probs <- 
  DOread::create_probs_query_func_do(datapath)
query_mrna <- 
  DOread::create_mrna_query_func_do(datapath)


