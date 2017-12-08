#' @export
read_project_rds <- function(project_info, dataname) {
  project <- project_info$project
  taxa <- project_info$taxa
  
  directory <- match(project, project_info$project)
  assertthat::assert_that(!is.na(directory))
  directory <- project_info$directory[directory]
  
  filepath <- file.path(directory, taxa, project, paste0(dataname, ".rds"))
  if(file.exists(filepath))
    readRDS(filepath)
  else
    NULL
}
