# Make chr_pos from chr, peak and window
#
# Create character string of chr:left-right with left=peak-window, right=peak+window
make_chr_pos <- function(chr_id=NULL, peak_Mbp=NULL, window_Mbp=NULL,
                         range = c(peak_Mbp - window_Mbp,
                                   peak_Mbp + window_Mbp)) {
  if(is.null(chr_id))
    chr_id <- "?"
  if(length(range) < 2 )
    range <- rep("?", 2)
  paste(chr_id, range[1], range[2], sep = "_")
}
