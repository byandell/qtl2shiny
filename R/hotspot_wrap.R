hotspot_wrap <- function(map, peaks, peak_window = 1, minLOD = 5.5,
                         project_info) {
  # This uses global hotspots
  if(peak_window == 1 & minLOD == 5.5) 
    read_project_rds(project_info, "hotspot")
  else
    qtl2pattern::hotspot(map, peaks, peak_window, minLOD)
}