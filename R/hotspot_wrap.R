hotspot_wrap <- function(map, peaks, peak_window = 1, minLOD = 5.5) {
  # This uses global hotspots
  if(peak_window == 1 & minLOD == 5.5) 
    hotspots
  else
    qtl2pattern::hotspot(map, peaks, peak_window, minLOD)
}