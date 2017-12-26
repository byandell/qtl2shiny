set_analyses <- function(dataset, data_group, dat) {
  ## Filter by dataset.

  # Start with phenotypes in data groups.
  if(shiny::isTruthy(data_group)) {
    dat <- dplyr::filter(dat, pheno_group %in% data_group)
  }
  # For identified datasets, drop other datasets from those data groups.
  if(shiny::isTruthy(dataset)) {
    if(!("all" %in% dataset)) {
      dat_sets <- dplyr::distinct(dat, pheno_type, pheno_group)
      dat_groups <- unique(dplyr::filter(dat_sets,
                                         pheno_type %in% dataset)$pheno_group)
      dat <- dplyr::filter(dat, 
                           (pheno_type %in% dataset) |
                             !(pheno_group %in% dat_groups))
    }
  }
  dat
}