# Summary across traits
summary_na <- function(phe) {
  tmpfn <- function(x) {
    out <- summary(x)
    out["NA's"] <- sum(is.na(x))
    out
  }
  t(sapply(phe, tmpfn))
}
