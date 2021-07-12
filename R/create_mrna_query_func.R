#' Create a function to query mRNA data
#'
#' Create a function that will connect to a database of mRNA information
#' and return a list with `probs` object and a `map` object.
#'
#' @param dbfile Name of database file
#' @param mrnadir_val name of mRNA data directory (default \code{"RNAseq"})
#'
#' @return Function with seven arguments, `chr`, `start`,
#'     `end`, `local`, `qtl`, `fast` and `mrnadir`. It returns a list with `expr`, `annot` and `peaks` objects
#'     spanning the region specified by the first three arguments.
#'
#' @details Note that this function assumes positions are in Mbp.
#'     There are required columns for each element, to be detailed in time.
#'     See \code{\link{read_mrna}} for details on how mRNA data are read.
#'     See \code{\link[qtl2]{create_variant_query_func}} for original idea.
#'
#' @export
create_mrna_query_func <- function(dbfile,
                                   mrnadir_val = "RNAseq") {
  if(missing(dbfile) || is.null(dbfile)) {
    # No mRNA data.
    function(chr = NULL, start = NULL, end = NULL,
             local = TRUE,
             qtl = FALSE,
             mrnadir = mrnadir_val) 
      NULL
  } else {
    function(chr = NULL, start = NULL, end = NULL,
             local = TRUE,
             qtl = FALSE,
             mrnadir = mrnadir_val)
      read_mrna(chr, start, end, dbfile, local, qtl, mrnadir)
  }
}
