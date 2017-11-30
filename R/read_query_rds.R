#' @export
#' @importFrom assertthat assert_that
read_query_rds <- function(project_info, filename) {
  projectfile <- file.path(project_info$directory,
                           project_info$taxa,
                           project_info$project,
                           filename)
  if(!file.exists(projectfile))
    projectfile <- file.path(project_info$directory,
                             project_info$taxa,
                             filename)
  assertthat::assert_that(file.exists(projectfile))
  readRDS(projectfile)
}