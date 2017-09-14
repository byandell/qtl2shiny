hotspot <- function(map, peaks, peak_window = 1, minLOD = 5.5) {
  # Set up list by chr of postions and peaks.
  round_pos <- function(x) {
    rng <- round(range(x))
    out <- seq(rng[1],rng[2])
    names(out) <- out
    out
  }
  map_pos <- purrr::map(map, round_pos)
  # Kludge
  for(chr in names(map_pos)) {
    names(map_pos[[chr]]) <- paste(chr, names(map_pos[[chr]]), sep = ":")
  }
  
  peaks <- dplyr::filter(peaks,
                         lod >= minLOD)
  if(!nrow(peaks))
    return(NULL)
  out_chr <- purrr::transpose(list(pos = map_pos, 
                                   peaks = split(peaks, peaks$chr)))

  peaks_type <- function(mapi, peaks, peak_window=1) {
    # count peaks at position by type
    peaks_type <- split(peaks, peaks$pheno_type)
    out <- as.data.frame(purrr::map(peaks_type, 
                                    outer_window, 
                                    mapi, peak_window))
    all <- as.matrix(out)
    out$all <- apply(all, 1, sum)
    rownames(out) <- mapi
    dplyr::select(out, all, dplyr::everything())
  }
  outer_window <- function(posi, mapi, peak_window = 1) {
    posi <- dplyr::filter(posi)$pos
    apply(outer(posi, mapi,
                function(x,y,z) abs(x-y) <= z,
                peak_window),
          2, sum)
  }
  
  out_peaks <- 
    dplyr::bind_rows(
      purrr::map(out_chr,
                 function(x, peak_window) peaks_type(x$pos, x$peaks, peak_window),
                 peak_window))
  out_peaks <- as.matrix(out_peaks)
  out_peaks[is.na(out_peaks)] <- 0
  rownames(out_peaks) <- unlist(sapply(map_pos, names))
  class(out_peaks) <- c("scan1", "matrix")

  out <- list(scan = out_peaks, map = map_pos)
  class(out) <- c("hotspot", "list")
  out
}
subset.hotspot <- function(x, chr, ...) {
  out <- list(scan = subset(x$scan, x$map, chr),
              map = x$map[chr])
  class(out) <- c("hotspot", "list")
  out
}
hotspot_old <-  function(map_chr, peak_window, pheno_types, map, peaks) {
  out_peaks <- NULL
  if(!is.null(peak_window)) {
    out_peaks <- matrix(0, length(unlist(map)),
                        length(pheno_types),
                        dimnames = list(qtl2scan:::map2markernames(map), pheno_types))
    
    index1 <- 1
    for(chri in map_chr) {
      mapi <- map[[chri]]
      index2 <- index1 + length(mapi) - 1
      for(phenoj in pheno_types) {
        if(phenoj=="all")
          posi <- peaks
        else
          posi <- dplyr::filter(peaks, pheno_type == phenoj)
        posi <- dplyr::filter(posi, chr==chri)$pos
        out_peaks[seq(index1, index2),phenoj] <-
          apply(outer(posi, mapi,
                      function(x,y,z) abs(x-y) <= z,
                      peak_window),
                2, sum)
      }
      index1 <- index2 + 1
    }
    class(out_peaks) <- c("scan1", "matrix")
  }
  out_peaks
}
