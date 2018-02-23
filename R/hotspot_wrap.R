hotspot_wrap <- function(map, peaks, peak_window = 1, minLOD = 5.5,
                         project_info) {
  # This uses global hotspots
  if(peak_window == 1 & minLOD == 5.5 & nrow(project_info)) 
    read_project(project_info, "hotspot")
  else {
    if(shiny::isTruthy(map) && shiny::isTruthy(peaks)) {
      qtl2pattern::hotspot(map, peaks, peak_window, minLOD)
    } else
      NULL
  }
}