##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title SE of propensity-paired differences on a fitted propensity score
##' @param object fitted propensity score model, of or inheriting from class \code{glm}
##' @param covariance.extractor function to extract covarance of fitted model coefficients
##' @param data a data frame
##' @param tt a terms object
##' @return scalar, interpretable as standard error
##' @author Mark M. Fredrickson, Ben B. Hansen
ppse <- function(object, covariance.extractor, data,...)
    UseMethod("ppse")

ppse.glm <- function(object, covariance.extractor=vcov, data=NULL, ...) 
{
    if (is.null(data))
        {
            form <- formula(object)
            form <- terms(form, specials="strata")
            data <- model.frame(object)
        }  
    ppse.default(object, covariance.extractor=covariance.extractor, data=data, tt=form,...)
}

ppse.default <- function(object, covariance.extractor=vcov, data=model.frame(object), tt=terms(object), ...)
    {
        if (is.null(names(coef(object)))) stop("propensity coefficients have to have names")

        data.matrix <- model.matrix(tt, data)
        
        stopifnot(!is.null(colnames(data.matrix)),
                  setequal(colnames(data.matrix), names(coef(object))),
                  !is.null(attr(data.matrix, "assign"))
                  )
        
        covb <- covariance.extractor(object)
        is.intercept <- attr(data.matrix, "assign")==0
        covb <- covb[!is.intercept, !is.intercept]
        covx <- cov(data.matrix)
        covx <- covx[!is.intercept, !is.intercept]
        coeffs <- coef(object)
        coeffs <- coeffs[!is.intercept]
        
        terms.to.sweep.out <- survival:::untangle.specials(tt, "strata")$terms ## strata, if present
        
        cols.to.keep <- !(attr(data.matrix, "assign") %in% terms.to.sweep.out)
        cols.to.keep <- cols.to.keep[!is.intercept]
        
        covb <- covb[cols.to.keep, cols.to.keep, drop=FALSE]
        coeffs <- coeffs[cols.to.keep]
        
        
        S <- if (!all(cols.to.keep))
        {
          n <- nrow(data.matrix)
          K <- sum(!cols.to.keep)
          
            (covx[cols.to.keep, cols.to.keep] -
                covx[cols.to.keep,!cols.to.keep, drop=FALSE] %*%
                    solve(covx[!cols.to.keep, !cols.to.keep, drop=FALSE],
                          covx[!cols.to.keep, cols.to.keep, drop=FALSE])
             ) *((n-1)/(n-K))
        } else covx
        
### `sandwich` gives me covs for `brglm`'s that lack dimnames.
### so, dropped below in favor of what comes after.
        ##  covb <- covb[,!colnames(covb) == "(Intercept)", drop = FALSE]
        ##  covb <- covb[!rownames(covb) == "(Intercept)",, drop = FALSE]
        coeff.not.NA <- is.finite(coeffs)
        coeffs <- coeffs[coeff.not.NA]

### with glms, NA coeffs won't have counterparts in cov matrix
### with e.g. coxph's that's not the case.  Make uniform.
if (!all(coeff.not.NA) && length(coeff.not.NA)==nrow(covb))
    covb <- covb[coeff.not.NA, coeff.not.NA, drop=FALSE]
S  <- S[coeff.not.NA, coeff.not.NA, drop=FALSE]

Sbetahat <- S %*% coeffs
varPShat <- crossprod(coeffs, Sbetahat)[1,1] # indexing to turn a 1x1 matrix into a scalar
Sperp <- S - (1/varPShat)*Sbetahat %*% t(Sbetahat)

## calculate the correction for the expected difference in propensity scores
sqrt(2 * sum(Sperp * covb) )
}
