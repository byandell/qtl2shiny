#' Normal scores nqrank
#'
#' Copied from github.com/kbroman/broman package
#' @importFrom stats sd runif qnorm
#' @export
nqrank <- function (x, jitter = FALSE)
{
  ## qtl::nqrank(x, jitter)
  y <- x[!is.na(x)]
  themean <- mean(y, na.rm = TRUE)
  thesd <- stats::sd(y, na.rm = TRUE)
  y[y == Inf] <- max(y[y < Inf]) + 10
  y[y == -Inf] <- min(y[y > -Inf]) - 10
  if (jitter)
    y <- rank(y + stats::runif(length(y)) / (stats::sd(y) * 10^8))
  else y <- rank(y)
  x[!is.na(x)] <- stats::qnorm((y - 0.5)/length(y))
  x * thesd / stats::sd(x, na.rm = TRUE) - mean(x, na.rm = TRUE) + themean
}
