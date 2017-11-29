hotspot_wrap <- function(map, peaks, peak_window = 1, minLOD = 5.5,
                         project_info) {
  # This uses global hotspots
  if(peak_window == 1 & minLOD == 5.5) 
    qtl2shiny::qtl2shiny_read(project_info$project, "hotspot")
  else
    qtl2pattern::hotspot(map, peaks, peak_window, minLOD)
}