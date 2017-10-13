hotspot_wrap <- function(map, peaks, peak_window = 1, minLOD = 5.5) {
  # This uses global datapath
  if(peak_window == 1 & minLOD == 5.5) 
    readRDS(file.path(datapath, "filtered", "hotspot.rds"))
  else
    DOread::hotspot(map, peaks, peak_window, minLOD)
}