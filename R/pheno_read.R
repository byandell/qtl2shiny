pheno_read <- function(project_info, analyses_df, transform = TRUE) {
  # Read the phenos we need.
  phenos <- analyses_df$pheno
  pheno_data <- read_project(project_info, "pheno_data", phenos)
  
  if(transform) {
    transform <- analyses_df$transf
  } else {
    transform <- NULL
  }
  pheno_trans(pheno_data, 
              phenos, 
              transform,
              analyses_df$offset,
              analyses_df$winsorize)
}