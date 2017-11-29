#' @export
#' @importFrom assertthat assert_that
qtl2shiny_read <- function(project, dataname) {
  directory <- match(project, projects$project)
  assertthat::assert_that(!is.na(directory))
  directory <- projects$directory[directory]
  readRDS(file.path(directory, project, paste0(dataname, ".rds")))
}
