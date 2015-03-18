##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title SE of propensity-paired differences on a fitted propensity score
##' @param propensity.glm fitted propensity score model, of or inheriting from class \code{glm}
##' @param covariance.extractor function to extract covarance of fitted model coefficients
##' @param data.matrix covariate matrix, ordinarily as retrieved by \code{model.matrix} (the default) 
##' @return scalar, interpretable as standard error
##' @author Mark M. Fredrickson, Ben B. Hansen
ppse <- function(propensity.glm, data.matrix=model.matrix(propensity.glm), covariance.extractor=vcov) {

#  stopifnot(inherits(propensity.glm, "glm")) #use at own risk 
  stopifnot(!is.null(colnames(data.matrix)),
            !is.null(names(coef(propensity.glm))),
            setequal(colnames(data.matrix), names(coef(propensity.glm)))
            )
  data <- data.matrix[,!colnames(data.matrix) == "(Intercept)", drop = FALSE]

  covx <- cov(data)

  covb <- covariance.extractor(propensity.glm)

### `sandwich` gives me covs for `brglm`'s that lack dimnames.
### so, dropped below in favor of what comes after.
#  covb <- covb[,!colnames(covb) == "(Intercept)", drop = FALSE]
#  covb <- covb[!rownames(covb) == "(Intercept)",, drop = FALSE]
  coeffs <- coef(propensity.glm)
  coeff.not.NA <- !is.na(coeffs)
  coeffs <- coeffs[coeff.not.NA]

### expecting that NA coeffs won't have counterparts in cov matrix  
  not.an.intercept <- names(coeffs) != "(Intercept)"
  covb <- covb[not.an.intercept, not.an.intercept]
  
### could we ever get a covb with different dim than covx?
### is this worth checking? --sure, when there are NA coefs.
### addressing that possibility as follows.  
  coeff.not.NA <- coeff.not.NA[not.an.intercept]
  covx <- covx[coeff.not.NA, coeff.not.NA]
  
  # calculate the correction for the expected difference in propensity scores
  ps <- predict(propensity.glm)
  rho.beta <- apply(data, 2, function(xi) { cor(xi, ps) })
  sd.x <- sqrt(diag(covx))

  srb <- sd.x * rho.beta
  
  sqrt(2 * sum(covx * covb) - 2 * (t(srb) %*% covb %*% srb)[1,1]) # indexing to turn a 1x1 matrix into a scalar
}
