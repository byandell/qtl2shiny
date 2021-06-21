# Copied from github.com/kbroman/broman package
winsorize <- function (x, q = 0.006) 
{
  assertthat::assert_that(is.numeric(x))
  assertthat::assert_that(assertthat::is.number(q), q >= 0, q <= 1)
  lohi <- stats::quantile(x, c(q, 1 - q), na.rm = TRUE)
  if (diff(lohi) < 0) 
    lohi <- rev(lohi)
  x[!is.na(x) & x < lohi[1]] <- lohi[1]
  x[!is.na(x) & x > lohi[2]] <- lohi[2]
  x
}