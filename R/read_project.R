read_project <- function(project_info, dataname, filetype = c("rds","csv","feather")) {
  filetype <- match.arg(filetype)
  
  if(!nrow(project_info))
    return(NULL)
  
  project <- project_info$project
  taxa <- project_info$taxa
  
  directory <- match(project, project_info$project)
  assertthat::assert_that(!is.na(directory))
  directory <- project_info$directory[directory]
  if(is.na(directory) || is.null(directory))
    directory <- "."
  
  filepath <- file.path(directory, taxa, project, paste0(dataname, ".rds"))
  if(!file.exists(filepath)) {
    filepath <- file.path(directory, taxa, paste0(dataname, ".rds"))
  }

  # Only RDS works for now. CSV will be easy, but feather ...  
  switch(filetype,
         rds = read_project_rds(project_info, dataname))
}