peaks_in_pos <- function(analyses, peaks, use_pos = TRUE,
                         chr_id=NULL, pos_Mbp=NULL, win=NULL) {
  phenames <- analyses$pheno
  peaks <- dplyr::filter(peaks,
                         pheno %in% phenames)
  if(use_pos) {
    if(!(is.null(chr_id) | is.null(pos_Mbp) | is.null(win))) {
      if(win > 0) {
        ## Filter peaks
        peaks <- dplyr::filter(peaks, 
                               chr == chr_id,
                               pos >= pos_Mbp - win,
                               pos <= pos_Mbp + win)
      }
    }
  }
  dplyr::mutate(
    dplyr::select(
      dplyr::distinct(
        dplyr::arrange(
          peaks,
          dplyr::desc(lod)),
        pheno, lod, pheno_group, pheno_type),
      pheno, lod, pheno_type, pheno_group),
    lod = round(lod, 1))
}
