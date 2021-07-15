#' Wrapper for 'qtl2' function fit1
#' 
#' @param driver genotype probabilities in matrix
#' @param target numeric vector of target values
#' @param ... additional parameters for \code{\link[qtl2]{fit1}}
#' 
#' @export
#' 
#' @importFrom qtl2 fit1
#' @importFrom stringr str_replace
#' 
fitQtl2 <- function(driver,
                    target,
                    ...) {

  if(is.null(rownames(driver)))
    rownames(driver) <- seq_len(nrow(driver))
  
  out <- qtl2::fit1(driver, target, ...)
  
  # Replace lod names with LR
  names(out) <- stringr::str_replace(names(out), "lod", "LR")
  names(out) <- stringr::str_replace(names(out), "_LR", "LR")

  # Rescael to make them likelihoods (or likelihood ratios)
  out$LR <- out$LR * log(10)
  out$indLR <- out$indLR * log(10)
  
  # Add df for later use
  out$df <- ncol(driver) - 1
  
  # Residuals
  fitted <- rep(NA, length(target))
  names(fitted) <- if(is.matrix(target)) {
    rownames(target)
  } else {
    names(target)
  }
  fitted[names(out$fitted)] <- out$fitted
  out$resid <- target - fitted
  
  out
}