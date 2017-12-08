projects <- read.csv("projects.csv", stringsAsFactors = FALSE)
read_project_data <- qtl2shiny::read_project_rds

# fixed directory approach
#source("inst/doqtl2/query_taxa.R")
#source("inst/doqtl2/query_project.R")

# Set up queries
# Want to have global queries that can be reset, or that take project info as arguments
# queries to get data
# queries to get genes, variants (for taxa)
# queries to get probs, mrna (for project)


