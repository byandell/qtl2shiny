# Summary across traits
summary_na <- function(phe) {
  tmpfn <- function(x) {
    out <- summary(x)
    out["NA's"] <- sum(is.na(x))
    out
  }
  out <- t(apply(phe, 2, tmpfn))
  data.frame(pheno = rownames(out), out, check.names = FALSE)
}
