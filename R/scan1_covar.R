#' scan1 for multiple traits with possibly different covariates
#' 
#' @param phe_df data frame of phenotypes
#' @param cov_df data frame of covariates
#' @param kinship kinship matrix or list of kinship matrices
#' @param analyses_df data frame of analyses information
#' 
#' @export
#' 
#' @importFrom qtl2scan scan1
#' @importFrom dplyr select
#'
scan1_covar <- function(phe_df, cov_df, probs_obj, kinship, analyses_df, ...) {
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
  scans <- scanfn(probs_obj, phe_df, kinship, cov_df, analyses_df, wh, models, ...)
  attr(scans, "hsq") <- NULL
  if(length(ucov) > 1) for(i in ucov[-1]) {
    wh <- which(covarset == i)
    tmp <- scanfn(probs_obj, phe_df, kinship, cov_df, analyses_df, wh, models, ...)
    attr(tmp, "hsq") <- NULL
    scans <- cbind(scans, tmp)
  }
  # reorder by decreasing max lod
  modify_object(scans, scans[,order(-apply(scans,2,max)), drop=FALSE])
}

scanfn <- function(probs_obj, phe_df, kinship, cov_df, analyses_df, wh, models,
                   sex_type = c("A","I","F","M","all"), ...) {
  
  sex_type <- match.arg(sex_type)
  if(sex_type == "all" & ncol(phe_df) > 1)
    sex_type <- "A"
  
  # scan1 for wh phenotypes using their covariates.
  phe_df <- phe_df[, wh, drop=FALSE]
  cov_df <- wh_covar(analyses_df, wh, cov_df)

  models <- models[wh]
  if(all(models == models[1])) {
    kinship <- if(models[1] == "binary") NULL else kinship
    if(ncol(cov_df)) {
      scansex(probs_obj, phe_df, kinship, cov_df,
              models[1], sex_type)
    } else { # no covariates (unlikely)
      qtl2scan::scan1(probs_obj, phe_df, kinship,
                      model = models[1])
    }
  } else { # multiple models
    umod <- unique(models)
    whm <- which(models == umod[1])
    kinship <- if(umod[1] == "binary") NULL else kinship
    if(ncol(cov_df)) {
      out <- scansex(probs_obj, phe_df[, whm, drop=FALSE], 
                     kinship, cov_df,
                     umod[1], sex_type)
      attr(out, "hsq") <- NULL
      for(mod in umod[-1]) {
        whm <- which(models == mod)
        kinship <- if(mod == "binary") NULL else kinship
        tmp <- scansex(probs_obj, phe_df[, whm, drop=FALSE], 
                       kinship, cov_df, mod, sex_type)
        attr(tmp, "hsq") <- NULL
        out <- cbind(out, tmp)
      }
    } else { # no covariates (unlikely)
      kinship <- if(umod[1] == "binary") NULL else kinship
      out <- qtl2scan::scan1(probs_obj, 
                             phe_df[, whm, drop=FALSE], 
                             kinship, 
                             model = umod[1])
      attr(out, "hsq") <- NULL
      for(mod in umod[-1]) {
        whm <- which(models == mod)
        kinship <- if(mod == "binary") NULL else kinship
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
  intcovar <- NULL
  if(!is.null(addcovar)) {
    if(sex_type == "all") {
      out <- scansex(genoprobs, pheno, kinship, addcovar, model, "A")
      out <- cbind(out,
                   scansex(genoprobs, pheno, kinship, addcovar, model, "I"))
      out <- cbind(out,
                   scansex(genoprobs, pheno, kinship, addcovar, model, "F"))
      out <- cbind(out,
                   scansex(genoprobs, pheno, kinship, addcovar, model, "M"))
      colnames(out) <- c("AddSex","IntSex","Female","Male")
      return(out)
    } else {
      intcovar <- NULL
      switch(sex_type,
             "I" = {
               intcovar <- model.matrix(~ sex, addcovar)[, -1, drop = FALSE]
             },
             "F" = {
               addcovar <- addcovar[addcovar$sex == "F",, drop = FALSE]
             },
             "M" = {
               addcovar <- addcovar[addcovar$sex == "M",, drop = FALSE]
             })
    }
    addcovar <- covar_df_mx(addcovar)
  }
  qtl2scan::scan1(genoprobs, pheno, kinship, addcovar, 
                  intcovar = intcovar,
                  model = model)
}
covar_df_mx <- function(addcovar) {
  if(is.data.frame(addcovar)) {
    f <- formula(paste("~", paste(names(addcovar), collapse = "+")))
    addcovar <- model.matrix(f, addcovar)[,-1, drop = FALSE]
    m <- match("sexM", colnames(addcovar))
    if(!is.na(m))
      colnames(addcovar)[m] <- "sex"
  }
  addcovar
}

