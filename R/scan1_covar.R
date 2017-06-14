#' scan1 for multiple traits with possibly different covariates
#' 
#' @export
#' 
#' @importFrom qtl2scan scan1
#' @importFrom dplyr select
#'
scan1_covar <- function(phe_df, cov_mx, probs_obj, K_chr, analyses_df, ...) {
  # This is set up for different types of models (e.g. "binary"),
  # but qtl2scan::scan1 ignores "binary" this right now if kinship is provided.
  # So below, force kinship to NULL if any "binary", 
  # and force attr "hsq" to null as well.
  models <- analyses_df$model
  analyses_df <- which_covar(analyses_df)
  ## Collapse to unique identifier for each row = each phenotype.
  covarset <- apply(analyses_df, 1, function(x) paste(1 * x, collapse = ""))
  ucov <- unique(covarset)

  wh <- which(covarset == ucov[1])
  scans <- scanfn(probs_obj, phe_df, K_chr, cov_mx, analyses_df, wh, models, ...)
  attr(scans, "hsq") <- NULL
  if(length(ucov) > 1) for(i in ucov[-1]) {
    wh <- which(covarset == i)
    tmp <- scanfn(probs_obj, phe_df, K_chr, cov_mx, analyses_df, wh, models, ...)
    attr(tmp, "hsq") <- NULL
    scans <- cbind(scans, tmp)
  }
  # reorder by decreasing max lod
  modify_object(scans, scans[,order(-apply(scans,2,max)), drop=FALSE])
}

scanfn <- function(probs_obj, phe_df, K_chr, cov_mx, analyses_df, wh, models,
                   sex_type = c("A","I","F","M","all"), ...) {
  
  sex_type <- match.arg(sex_type)
  if(sex_type == "all" & ncol(phe_df) > 1)
    sex_type <- "A"
  
  # scan1 for wh phenotypes using their covariates.
  phe_df <- phe_df[, wh, drop=FALSE]
  covars <- unlist(analyses_df[wh[1],])
  cov_mx <- cov_mx[, covars, drop=FALSE]
  models <- models[wh]
  if(all(models == models[1])) {
    kinship <- if(models[1] == "binary") NULL else K_chr
    if(length(covars)) {
      scansex(probs_obj, phe_df, kinship, cov_mx,
              models[1], sex_type)
    } else { # no covariates (unlikely)
      qtl2scan::scan1(probs_obj, phe_df, kinship,
                      model = models[1])
    }
  } else { # multiple models
    umod <- unique(models)
    whm <- which(models == umod[1])
    kinship <- if(umod[1] == "binary") NULL else K_chr
    if(length(covars)) {
      out <- scansex(probs_obj, phe_df[, whm, drop=FALSE], 
                     kinship, cov_mx,
                     umod[1], sex_type)
      attr(out, "hsq") <- NULL
      for(mod in umod[-1]) {
        whm <- which(models == mod)
        kinship <- if(mod == "binary") NULL else K_chr
        tmp <- scansex(probs_obj, phe_df[, whm, drop=FALSE], 
                       kinship, cov_mx, mod, sex_type)
        attr(tmp, "hsq") <- NULL
        out <- cbind(out, tmp)
      }
    } else { # no covariates (unlikely)
      kinship <- if(umod[1] == "binary") NULL else K_chr
      out <- qtl2scan::scan1(probs_obj, 
                             phe_df[, whm, drop=FALSE], 
                             kinship, 
                             model = umod[1])
      attr(out, "hsq") <- NULL
      for(mod in umod[-1]) {
        whm <- which(models == mod)
        kinship <- if(mod == "binary") NULL else K_chr
        tmp <- qtl2scan::scan1(probs_obj, 
                               phe_df[, whm, drop=FALSE], 
                               kinship, 
                               model = mod)
        attr(tmp, "hsq") <- NULL
        out <- cbind(out, tmp)
      }
    }
    out
  }
}

scansex <- function(genoprobs, pheno, kinship, addcovar = NULL,
                    model, sex_type) {
  if(sex_type == "all") {
    out <- qtl2scan::scan1(genoprobs, pheno, kinship, 
                           addcovar,
                           model = model)
    out <- cbind(out,
                 qtl2scan::scan1(genoprobs, pheno, kinship, 
                                 addcovar, 
                                 intcovar = addcovar[,"sex", drop = FALSE],
                                 model = model))
    out <- cbind(out,
                 qtl2scan::scan1(genoprobs, pheno, kinship, 
                                 addcovar[addcovar[, "sex"] == 0,, drop = FALSE],
                                 model = model))
    out <- cbind(out,
                 qtl2scan::scan1(genoprobs, pheno, kinship, 
                                 addcovar[addcovar[, "sex"] == 1,, drop = FALSE],
                                 model = model))
    colnames(out) <- c("AddSex","IntSex","Female","Male")
    out
  } else {
    intcovar <- NULL
    switch(sex_type,
           "I" = {
             intcovar <- addcovar[,"sex", drop = FALSE]
           },
           "F" = {
             addcovar <- addcovar[addcovar[, "sex"] == 0,, drop = FALSE]
           },
           "M" = {
             addcovar <- addcovar[addcovar[, "sex"] == 1,, drop = FALSE]
           })
    qtl2scan::scan1(genoprobs, pheno, kinship, addcovar, 
                    intcovar = intcovar,
                    model = model)
  }
}

