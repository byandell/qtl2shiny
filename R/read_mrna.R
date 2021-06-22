#' Read mrna expression object from file
#'
#' Uses fst to read mrna expression object and associated annotations.
#'
#' @param chr_id vector of chromosome identifiers
#' @param start_val,end_val start and end values in Mbp
#' @param datapath name of folder with Derived Data
#' @param local read only mRNA values local to region if \code{TRUE};
#' otherwise include distal mRNA values that map to region
#' @param qtl read only mRNA values with QTL peak in region if \code{TRUE}
#' @param mrnadir name of directory with mRNA data
#'
#' @details Reads `expr`, `peaks` and `annot` information on mRNA and combines into a list.
#'     The `expr` and `peaks` elements are stored in `fst` database,
#'     while `peaks` is an `RDS` database. 
#'     
#' @return list with \code{expr} = matrix of expression mRNA values in region and \code{annot} = data frame of annotations for mRNA.
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @importFrom dplyr filter group_by inner_join mutate n rename summarize ungroup
#' @importFrom fst read_fst
#' @importFrom rlang .data
#'
read_mrna <- function(chr_id=NULL, start_val=NULL, end_val=NULL, datapath,
                      local = TRUE, qtl = FALSE, 
                      mrnadir = "RNAseq") {

  if(is.null(chr_id))
    stop("must supply chr_id")
  if(is.null(start_val))
    start_val <- 0
  if(is.null(end_val))
    end_val <- Inf
  
  # Identify mRNA located in region or with QTL peak in region.
  peaks.mrna <- fst::read_fst(file.path(datapath, mrnadir, "peaks.mrna.fst"))
  if(local) {
    peaks.mrna <- dplyr::filter(peaks.mrna,
                                .data$gene_chr == chr_id,
                                pmax(.data$gene_start, .data$gene_end) >= start_val,
                                pmin(.data$gene_start, .data$gene_end) <= end_val)
  } else {
    peaks.mrna <- dplyr::filter(peaks.mrna,
                                ((.data$gene_chr == chr_id &
                                    pmax(.data$gene_start, .data$gene_end) >= start_val &
                                    pmin(.data$gene_start, .data$gene_end) <= end_val) |
                                   (.data$qtl_chr == chr_id &
                                      .data$qtl_pos >= start_val &
                                      .data$qtl_pos <= end_val)))
  }
  mrna_ids <- unique(peaks.mrna$gene_id)

  # Get annotations for unique mRNA IDs
  annot.mrna <-
    dplyr::rename(
      dplyr::filter(
        readRDS(file.path(datapath, mrnadir, "annot.mrna.rds")),
        .data$id %in% mrna_ids),
      pos = .data$middle_point)

  annot.mrna <-
    dplyr::mutate(
      dplyr::inner_join(
        annot.mrna,
        dplyr::rename(
          dplyr::ungroup(
            dplyr::summarize(
              dplyr::group_by(peaks.mrna, .data$gene_id),
              qtl_ct = dplyr::n(),
              info = paste0(.data$qtl_chr, "@",
                            round(.data$qtl_pos), ":",
                            round(.data$lod), collapse = ","),
              qtl_pos = ifelse(any(.data$qtl_chr == chr_id &
                                     .data$qtl_pos >= start_val &
                                     .data$qtl_pos <= end_val),
                               .data$qtl_pos[.data$qtl_chr == chr_id], NA),
              lod = ifelse(is.na(.data$qtl_pos),
                           NA, .data$lod))),
          id = .data$gene_id),
        by = "id"),
      local = !is.na(.data$qtl_pos) &
        .data$chr == chr_id &
        pmax(.data$start, .data$end) >= start_val &
        pmin(.data$start, .data$end) <= end_val)

  if(qtl) {
    annot.mrna <- dplyr::filter(annot.mrna, !is.na(.data$qtl_pos))
    # Reduce to mRNA having peak in region.
    expr_id <- annot.mrna$id
    peaks.mrna <- dplyr::filter(peaks.mrna, .data$gene_id %in% expr_id)
  } else {
    expr_id <- annot.mrna$id
  }
  if(nrow(annot.mrna) == 0 || nrow(peaks.mrna) == 0)
    return(NULL)

  # Get expression data.
  expr.mrna <- read_fast(file.path(datapath, mrnadir, "expr.mrna.fst"),
                         expr_id, rownames = TRUE)

  list(expr = expr.mrna, annot = annot.mrna, peaks = peaks.mrna)
}
