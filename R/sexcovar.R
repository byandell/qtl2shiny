sexcovar <- function(addcovar, sex_type) {
  if(!all(sort(unique(addcovar$sex)) %in% c("M","F")) & sex_type %in% c("M","F")) {
    cat("cannot handle levels of sex not M and F", file = stderr())
    return(NULL)          
  }
  switch(sex_type,
         "F" = addcovar <- addcovar[addcovar$sex == "F",, drop = FALSE],
         "M" = addcovar <- addcovar[addcovar$sex == "M",, drop = FALSE])
  if(sex_type %in% c("F","M")) {
    addcovar <- dplyr::select(addcovar, -sex)
    if(ncol(addcovar) == 0)
      addcovar <- NULL
  }
  addcovar
}