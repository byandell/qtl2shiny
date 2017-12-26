select_phenames <- function(phenames, peaks_df) {
  selected <- phenames
  if(shiny::isTruthy(peaks_df)) {
    phenames <- unique(peaks_df$pheno)
    # Limit to first 1000
    nphe <- length(phenames)
    phenames <- phenames[seq_len(min(1000, nphe))]
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
  
  choices <- c("all","none", phenames)
  label = ifelse(nphe <= 1000,
                 "Choose phenotypes",
                 paste("Top 1000 of", nphe))
  list(label = label, choices = choices, selected = selected)
}