hotspot <- function(map, peaks, peak_window = 1, minLOD = 5.5) {
  # Set up list by chr of postions and peaks.
  round_pos <- function(x) {
    rng <- round(range(x))
    out <- seq(rng[1],rng[2])
    names(out) <- out
    out
  }
  chr_pos <- purrr::map(map, round_pos)
  # Kludge
  for(chr in names(chr_pos)) {
    names(chr_pos[[chr]]) <- paste(chr, names(chr_pos[[chr]]), sep = ":")
  }
  
  peaks <- dplyr::filter(peaks,
                         lod >= minLOD)
  if(!nrow(peaks))
    return(NULL)
  out_chr <- purrr::transpose(list(pos = chr_pos, 
                                   peaks = split(peaks, peaks$chr)))

  peaks_type <- function(posi, peaks, peak_window=1) {
    # count peaks at position by type
    peaks_by_type <- split(peaks, peaks$pheno_type)
    out <- data.frame(purrr::map(peaks_by_type, 
                                 outer_window, 
                                 posi, peak_window),
                      check.names = FALSE)
    if(!nrow(out))
      return(NULL)
    groups <- dplyr::distinct(peaks, pheno_group, pheno_type)
    groups <- split(groups$pheno_type, groups$pheno_group)
    grps <- data.frame(
      purrr::map(groups,
                 function(x, out) {
                   apply(out[, x, drop = FALSE], 1, sum)
                   },
                 out),
      check.names = FALSE)
    out$all <- apply(out, 1, sum)
    out <- dplyr::bind_cols(out, grps)
    if(max(out) == 0)
      return(NULL)
    rownames(out) <- posi
    dplyr::select(out, all, 
                  dplyr::one_of(names(groups)),
                  dplyr::everything())
  }
  outer_window <- function(peaksi, posi, peak_window = 1) {
    peaksi <- dplyr::filter(peaksi)$pos
    apply(outer(peaksi, posi,
                function(x,y,z) abs(x-y) <= z,
                peak_window),
          2, sum)
  }
  
  # Want to identify what purrr::map are NULL and adjust map
  out_peaks <- purrr::map(out_chr,
                          function(x, peak_window) peaks_type(x$pos, x$peaks, peak_window),
                          peak_window)
  chr_pos <- chr_pos[!sapply(out_peaks, is.null)]
  if(!length(chr_pos))
    return(NULL)
  out_peaks <- dplyr::bind_rows(out_peaks)
  if(!nrow(out_peaks))
    return(NULL)
  
  out_peaks <- as.matrix(out_peaks)
  out_peaks[is.na(out_peaks)] <- 0
  # Breaks here when threshold is too large.
  rownames(out_peaks) <- unlist(sapply(chr_pos, names))
  class(out_peaks) <- c("scan1", "matrix")

  out <- list(scan = out_peaks, map = chr_pos)
  class(out) <- c("hotspot", "list")
  out
}
subset.hotspot <- function(x, chr = NULL, nonzero = NULL, ...) {
  if(!is.null(chr)) {
    x <- list(scan = subset(x$scan, x$map, chr),
                map = x$map[chr])
    class(x) <- c("hotspot", "list")
  }
  # drop chr with all zeroes
  if(!is.null(nonzero)) {
    cts <- apply(x$scan[, nonzero, drop = FALSE], 1, sum)
    chrs <- unlist(tapply(cts, 
                          ordered(rep(names(x$map), sapply(x$map, length)),
                                  names(x$map)),
            function(x) !all(x == 0)))
    if(!all(chrs))
      x <- subset(x, chr = names(chrs)[chrs])
  }
  x
}
##########################################################
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
