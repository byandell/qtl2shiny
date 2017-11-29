projects <- read.csv("projects.csv", stringsAsFactors = FALSE)
project <- projects$project[1]
id <- 1
datapath <- file.path(projects$directory[id], project)

# Need to think about how to hand off from one project to next.
# Not so easy, as must reside within app rather than global...

#pmap         <- qtl2shiny::qtl2shiny_read(project, "pmap")
#K            <- qtl2shiny::qtl2shiny_read(project, "kinship")

#covar        <- qtl2shiny::qtl2shiny_read(project, "covar")
#peaks        <- qtl2shiny::qtl2shiny_read(project, "peaks")
#analyses_tbl <- qtl2shiny::qtl2shiny_read(project, "analyses")
#pheno_data   <- qtl2shiny::qtl2shiny_read(project, "pheno")
#pheno_type   <- qtl2shiny::qtl2shiny_read(project, "pheno_type")
#hotspots     <- qtl2shiny::qtl2shiny_read(project, "hotspot")

# Set up queries
query_genes <- 
  qtl2db::create_gene_query_func(
    file.path(datapath, projects$genes[id]))
query_variants <- 
  qtl2db::create_variant_query_func(
    file.path(datapath, projects$variants[id]))
query_probs <- 
  DOread::create_probs_query_func_do(datapath)
query_mrna <- 
  DOread::create_mrna_query_func_do(datapath)


