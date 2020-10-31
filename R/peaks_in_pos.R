#' @importFrom dplyr arrange desc distinct filter mutate select
#' @importFrom rlang .data
#' 
peaks_in_pos <- function(analyses, peaks, use_pos = TRUE,
                         chr_id=NULL, pos_Mbp=NULL, win=NULL) {
  phenames <- analyses$pheno
  peaks <- dplyr::filter(peaks,
                         .data$pheno %in% phenames)
  if(use_pos) {
    if(!(is.null(chr_id) | is.null(pos_Mbp) | is.null(win))) {
      if(win > 0) {
        ## Filter peaks
        peaks <- dplyr::filter(peaks, 
                               .data$chr == chr_id,
                               .data$pos >= pos_Mbp - win,
                               .data$pos <= pos_Mbp + win)
      }
    }
  }
  dplyr::mutate(
    dplyr::select(
      dplyr::distinct(
        dplyr::arrange(
          peaks,
          dplyr::desc(.data$lod)),
        .data$pheno, .data$lod, .data$pheno_group, .data$pheno_type, .data$chr, .data$pos),
      .data$pheno, .data$lod, .data$pheno_type, .data$pheno_group, .data$chr, .data$pos),
    lod = round(.data$lod, 1),
    pos = round(.data$pos, 2))
}
