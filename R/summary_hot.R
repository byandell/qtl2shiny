#' @importFrom dplyr rename select
#' @importFrom rlang .data
#' 
summary_hot <- function(peak_set, scan_obj) {
  # Used chosen datasets, or all if not chosen.

  map <- scan_obj$map
  scan <- scan_obj$scan
  
  # Match lod columns to those present.
  lodcol <- match(peak_set, colnames(scan))
  lodcol <- lodcol[!is.na(lodcol)]

  if(length(lodcol) & (nrow(scan) == length(unlist(map)))) {
    chr_id <- names(map)
    dplyr::select(
      dplyr::rename(
        summary(
          subset(scan, lodcolumn = lodcol),
          map, chr = chr_id),
        count = .data$lod),
      -.data$marker)
  } else {
    NULL
  }
}