select_phenames <- function(phenames, peaks_df, local,
                            chr_id, peak_Mbp, window_Mbp) {
  selected <- phenames
  if(shiny::isTruthy(peaks_df) && nrow(peaks_df)) {
    peaks_df <- dplyr::filter(peaks_df, 
                              chr == chr_id)
    if(shiny::isTruthy(local)) {
      peaks_df <- dplyr::filter(peaks_df, 
                                pos >= peak_Mbp - window_Mbp,
                                pos <= peak_Mbp + window_Mbp)
    }
    phenames <- dplyr::distinct(
      dplyr::arrange(
        peaks_df,
        dplyr::desc(lod)),
      pheno)$pheno
  }

  if("all" %in% selected)
    selected <- c(selected[!(selected %in% c("all","none"))],
                  phenames)
  if("none" %in% selected)
    selected <- ""
  if(!is.null(selected)) {
    selected <- sort(unique(selected))
    selected <- selected[selected %in% phenames]
  }
  
  ## Update phenames to include selected (but not "")
  phenames <- unique(c(selected, phenames))
  phenames <- phenames[phenames != ""]
  # Limit to first 1000
  nphe <- length(phenames)
  phenames <- phenames[seq_len(min(1000, nphe))]
  
  choices <- c("all","none", phenames)
  label = ifelse(nphe <= 1000,
                 "Choose phenotypes",
                 paste("Top 1000 of", nphe))
  list(label = label, choices = choices, selected = selected)
}