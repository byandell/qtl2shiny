# Want alternative table of top phenotypes if only one chr.
# Use set_par$dataset and input$chr_ct
# length(chr_ct) == 1 & !("all" %in% chr_ct)
# dplyr::filter(peaks_tbl(), chr == chr_ct, pheno_type %in% dataset)
# dplyr::select(pheno, pheno_type, chr, pos, lod)
peakDataTable <- function(scan_tbl, peaks_tbl) {
  if(length(chrs <- unique(scan_tbl$chr)) > 1) {
    dplyr::arrange(scan_tbl, desc(count))
  } else {
    pheno_types <- make.names(unique(scan_tbl$pheno))
    peaks_tbl <- dplyr::filter(peaks_tbl, chr %in% chrs)
    if("all" %in% pheno_types)
      pheno_types <- make.names(unique(peaks_tbl$pheno_type))
    dplyr::arrange(
      dplyr::select(
        dplyr::filter(peaks_tbl,
                      make.names(pheno_type) %in% pheno_types),
        pheno, pheno_type, chr, pos, lod),
      desc(lod))
  }
}