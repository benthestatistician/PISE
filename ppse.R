##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title SE of propensity-paired differences on a fitted propensity score
##' @param object fitted propensity score model, of or inheriting from class \code{glm}
##' @param covariance.extractor function to extract covarance of fitted model coefficients
##' @param data.matrix covariate matrix, ordinarily as retrieved by \code{model.matrix} (the default) 
##' @return scalar, interpretable as standard error
##' @author Mark M. Fredrickson, Ben B. Hansen
ppse <- function(object, covariance.extractor, data.matrix,...)
    UseMethod("ppse")

ppse.glm <- function(object, covariance.extractor=vcov, data.matrix=NULL, ...) 
{
    if (is.null(data.matrix))
        {
            form <- formula(object)
            form <- terms(form, specials="strata")
            mf <- model.frame(object)
            data.matrix <- model.matrix(form, mf)
        }  
    ppse.default(object, covariance.extractor=covariance.extractor, data.matrix=data.matrix,...)
}

ppse.default <- function(object, covariance.extractor=vcov, data.matrix=NULL, ...)
    {
        if (is.null(names(coef(object)))) stop("propensity coefficients have to have names")
        if (is.null(data.matrix)) {
            mf <- model.frame(object)
            tt <- attr(mf, "terms")
            data.matrix <- model.matrix(tt, mf)
        }
        stopifnot(!is.null(colnames(data.matrix)),
                  setequal(colnames(data.matrix), names(coef(object)))
                  )
  covb <- covariance.extractor(object)

### `sandwich` gives me covs for `brglm`'s that lack dimnames.
### so, dropped below in favor of what comes after.
#  covb <- covb[,!colnames(covb) == "(Intercept)", drop = FALSE]
#  covb <- covb[!rownames(covb) == "(Intercept)",, drop = FALSE]
  coeffs <- coef(object)
  coeff.not.NA <- is.finite(coeffs)
  coeffs <- coeffs[coeff.not.NA]

### with glms, NA coeffs won't have counterparts in cov matrix
### with e.g. coxph's that's not the case.  Make uniform.
  if (!all(coeff.not.NA) && length(coeff.not.NA)==nrow(covb))
    covb <- covb[coeff.not.NA, coeff.not.NA, drop=FALSE]
  
  not.an.intercept <- names(coeffs) != "(Intercept)"
  covb <- covb[not.an.intercept, not.an.intercept, drop=FALSE]
    coeffs <- coeffs[not.an.intercept]
    
  data <- data.matrix[,coeff.not.NA & colnames(data.matrix) != "(Intercept)",
                      drop = FALSE]

    covx <- cov(data)

    Sbetahat <- covx %*% coeffs
    varPShat <- crossprod(coeffs, Sbetahat)[1,1] # indexing to turn a 1x1 matrix into a scalar
    Sperp <- covx - (1/varPShat)*Sbetahat %*% t(Sbetahat)
  # calculate the correction for the expected difference in propensity scores
  
  sqrt(2 * sum(Sperp * covb) )
}
