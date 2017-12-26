changing project_info should change multiple things

- chr and pos (based on hotspots)
- put project name in upper left

need to make sure pheno names are reset when project changes
  does not seem to do updateSelectInput properly
  
hotspot_wrap has minLOD = 5.5 hardwired but Recla min is 6?

changing project without visiting phenotypes can lead to errors later on with phenotype not in set

may need to move phenos into setup: rethink phenos, peaks, setup