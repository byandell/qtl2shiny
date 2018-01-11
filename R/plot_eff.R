# Plot of effects.
plot_eff <- function(pheno, eff_obj, map, scan_obj, xlim = NULL,
                     addlod = FALSE, allele_info) {
  if(is.null(eff_obj) | is.null(pheno) | is.null(scan_obj))
    return(NULL)
  
  main <- pheno
  effcol <- match(pheno, names(eff_obj))
  if(any(is.na(effcol)))
    return(plot_null("effect name mismatch"))
  
  lodcol <- match(pheno, colnames(scan_obj))
  if(any(is.na(lodcol))) # Probably all sex
    lodcol <- match("AddSex", colnames(scan_obj))
  if(any(is.na(lodcol)))
    return(plot_null("scan name mismatch"))
  colnames(scan_obj)[lodcol] <- pheno
  
  # Pull colors and shortnames from allele_info
  colors <- allele_info$color
  names(colors) <- allele_info$shortname
  
  chr_id <- names(map)[1]
  if(!addlod) {
    max_pos <- max(scan_obj, map, lodcolumn=lodcol)$pos[1]
    ggplot2::autoplot(
      eff_obj[[effcol]], map,
      xlim=xlim,
      main = main,
      legend.position = "right",
      maxpos = max_pos, maxcol = lodcol,
      colors = colors)
  } else { # coef_and_lod
    ggplot2::autoplot(
      eff_obj[[effcol]], map,
      scan1_output = subset(scan_obj, map,
                            chr = chr_id,
                            lodcolumn = lodcol),
      xlim = xlim,
      maxcol = lodcol,
      legend.position = "none",
      legend.position_lod = "none",
      main = main,
      colors = colors)
  }
}
