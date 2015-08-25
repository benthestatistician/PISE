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
            covx[cols.to.keep, cols.to.keep] -
                covx[cols.to.keep,!cols.to.keep, drop=FALSE] %*%
                    solve(covx[!cols.to.keep, !cols.to.keep, drop=FALSE],
                          covx[!cols.to.keep, cols.to.keep, drop=FALSE])
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

Sperp <- makeSperp(S, betas=coeffs)

## calculate the correction for the expected difference in propensity scores
sqrt(2 * sum(Sperp * covb) )
}
##' Covariance of projections of X onto X'beta, calculated from Cov(X) and beta
##'
##' A key quantity figuring in paired standard error calculations
##' @title Calculate S-perp from S and betas
##' @param S covariance of X
##' @param betas coefficients
##' @return positive definite matrix
##' @author Ben B Hansen
makeSperp <- function(S, betas) {
    stopifnot(is.matrix(S), nrow(S)==ncol(S), is.numeric(betas), nrow(S)==length(betas))
    Sbetas <- S %*% betas
    varPShat <- crossprod(betas, Sbetas)[1,1] # indexing to turn a 1x1 matrix into a scalar
    S - (1/varPShat)*Sbetas %*% t(Sbetas)
}

getglmQweights <- function(eta, prior.weights=NULL, family=binomial())
    {
        nobs <- length(eta)
        if (is.null(prior.weights)) 
            prior.weights <- rep.int(1, nobs)
        good <- prior.weights >0
        variance <- family$variance
        linkinv <- family$linkinv
        if (!is.function(variance) || !is.function(linkinv)) 
            stop("'family' argument seems not to be a valid family object", 
                 call. = FALSE)
        mu.eta <- family$mu.eta
        mu <- linkinv(eta)
        mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
        good <- prior.weights >0 & mu.eta.val !=0
        ifelse(good,prior.weights*mu.eta.val^2/variance(mu),0)
    }
ppse.qr <- function(object, covariance.extractor=vcov, data=NULL, fitted.model, ...) 
{
    stopifnot(inherits(fitted.model, "glm"))
    if (!identical(covariance.extractor,vcov)) stop('ppse.qr only supports vcov as covariance extractor')
    tt <- terms(fitted.model)
##    stopifnot(!is.null(fitted.model$weights))
    ## b/c I'm not sure how this interacts w/ factor expansion in next few lines
    if (is.null(data)) data <- model.frame(fitted.model)
    stopifnot(is.null(model.offset(data))) # assume away offsets (for now!)    
    offset <- rep(0, nrow(data))
    eta <- offset +
        fitted.model$linear.predictors # lazy, for now. 
    weights <- getglmQweights(eta, # refigure weights, since glm.fit doesn't update them
                              prior.weights=fitted.model$prior.weights, # after final iteration
                              family=fitted.model$family)
    stopifnot(length(weights)==nrow(data))
    w <- sqrt(weights) 
    good <- !is.na(w) & w>0
    
    data.matrix <- model.matrix(tt, data)
    resids <- fitted.model$residuals
##    data.matrix <- data.matrix[good,]
    resids <- resids[good]
    eta <- eta[good]
    w <- w[good]
    z <- eta + resids # "z" as in `glm.fit`
    

## Next 2 lines seem to be idle...   
###    if (ncol(data.matrix) !=length(object$pivot)) stop("length of QR's pivot doesn't match ncol(data.matrix).")
###    data.matrix <- data.matrix[,object$pivot]
## Next calc reflected erroneous linear algebra on my part.    
###    qtilde.prime <- backsolve(qr.R(object), # qtilde ~= XR^(-1), where X is actual data matrix, XW = QR
###                        x=t(data.matrix), # specifically, qtilde solves
###                              transpose=TRUE) # X= qtilde %*% R
    qmat <- qr.Q(object)[good,]
    ##    qcoeffs <- qr.qty(object,z*w) # doesn't work, returns a object of length length(z).  Likewise for qr.qy(...)
    qcoeffs <- crossprod(qmat, z*w)

##  qtilde is the reweighting of the Q-matrix that corresponds to a rotated X. (Q ~ rotated W*X)    
    qtilde <- w^(-1) * qmat
    covqtilde <- cov(qtilde)


    Sqperp <-  makeSperp(covqtilde, qcoeffs)

    ## On to estimating (co)variance of the qcoeffs...
    nobs <- sum(good)
    rdf <- nobs - object$rank
###    stopifnot(all.equal(fitted.model$linear.predictors[good], drop(qtilde %*% qcoeffs))) # fails...
    ## but discrepancy is plausibly due fact that fitted.model$weights isn't updated after final iteration.
    ## it would be nice to have a cross-check against glm results, but I don't see how to do that.
        browser()

    ## if we later tinker w/ eta based on the QR, we'll have to adjust z, w in the below
    qfitted <- qmat %*% qcoeffs
    rss <- sum((z*w - qfitted)^2)
    varqcoeffs <- rss/rdf # q cols have sums of squares = to 1, so division by that is implicit
    sqrt(2 * varqcoeffs*sum(diag(Sqperp)))
}
