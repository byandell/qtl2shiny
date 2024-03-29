#' @importFrom utils packageVersion
#' 
versions <- function(packages = c(qtl2 = "qtl2",
                                  shiny = "qtl2shiny",
                                  pattern = "qtl2pattern",
                                  ggplot = "qtl2ggplot",
                                  fst = "qtl2fst",
                                  mediate = "qtl2mediate",
                                  intermediate = "intermediate")) {
  pv <- unlist(apply(as.matrix(packages), 1, function(x) as.character(utils::packageVersion(x))))
  paste("Version", paste(names(packages), pv, sep = "=", collapse = ", "))
}
