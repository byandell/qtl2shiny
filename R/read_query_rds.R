#' Read Query for RDS Object
#' 
#' @param project_info table of project information
#' @param filename name of RDS file
#' 
#' @export
read_query_rds <- function(project_info, filename) {
  if(is.null(project_info))
    return(NULL)
  datapath <- project_info$directory
  projectfile <- file.path(datapath,
                           project_info$taxa,
                           project_info$project,
                           filename)
  if(!file.exists(projectfile))
    projectfile <- file.path(datapath,
                             project_info$taxa,
                             filename)
  assertthat::assert_that(file.exists(projectfile))
  readRDS(projectfile)
}