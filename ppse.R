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
ppse.qr <- function(object, covariance.extractor=vcov, data=NULL, fitted.model,
                    coeffs.from.fitted.model=FALSE, tol.coeff.alignment=Inf,
                    return.extras=FALSE,...) 
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

    ## To better reproduce internals of glm, I tried:
    ## (1) using the weights based on penultimate glm fit,
    ## not the updated weights as figured by my `getglmQweights`.
###    weights <- fitted.model$weights
    ## (2) insisting on convergence,
###    stopifnot(fitted.model$converged)
    ## Rationale for requiring convergence: Since we don't also have access to
    ## penultimate eta's, we're relying on the convergence being far enough along that their differences
    ## from the final etas are numerically negligible.  I get the impression that the `glm.fit` convergence
    ## criteria are such as to ensure this.
    ## None of this made much difference, with the aglm example used in optmatch, I got differences
    ## in coefficients of order 2e-6 either way. Chalking the discrepancy up (uneasily) to differences
    ## between numerical calcs used within `C_Cdqrls` and the QR-based method I'm using. 
    
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
    ## NB: successful `getglmQweights` confirms that `linkinv` is there and is a function
    linkinv <- fitted.model$family$linkinv
    mu <- linkinv(eta)
    y <- model.response(data, type="double")
    y <- y[good]
    mu.eta <- fitted.model$family$mu.eta
    mu.etaval <- mu.eta(eta)
    resids <- (y-mu)/mu.etaval
    
    z <- (eta -offset) + resids # "z" as in `glm.fit`
    

## Next 2 lines seem to be idle...   
###    if (ncol(data.matrix) !=length(object$pivot)) stop("length of QR's pivot doesn't match ncol(data.matrix).")
###    data.matrix <- data.matrix[,object$pivot]
    qmat <- qr.Q(object)[good,]
    ## NB: `qcoeffs <- qr.qty(object,z*w)` doesn't work, returns a object of length length(z),
    ## not object$rank.  Likewise for qr.qy(...)
    qcoeffs.from.QR <- crossprod(qmat, z*w)
    qcoeffs.from.fitted.model <- drop(qr.R(object)%*%coef(fitted.model)[object$pivot])
    if (is.finite(tol.coeff.alignment))
        {
            coeff.diffs <- qcoeffs.from.fitted.model - qcoeffs.from.QR
            if (max(abs(coeff.diffs)) >= tol.coeff.alignment)
            stop(paste("QR/reported coefficients differ by up to", prettyNum(abs(coeff.diffs))))
        }
    qcoeffs <- if (coeffs.from.fitted.model) qcoeffs.from.fitted.model else qcoeffs.from.QR

    ##  qtilde is the reweighting of the Q-matrix that corresponds to a rotated X. (Q ~ rotated W*X)    
    qtilde <- w^(-1) * qmat
    covqtilde <- cov(qtilde)

    Sqperp <-  makeSperp(covqtilde, qcoeffs)

    ## On to estimating (co)variance of the qcoeffs...
    nobs <- sum(good)
    stopifnot((rdf <- nobs - object$rank)>0) # i.e., fitted.model$df.residual
    
    est.disp <- !(fitted.model$family$family %in% c("poisson", "binomial"))
    dispersion <-
        if (est.disp)
        {              
            qfitted <- qmat %*% qcoeffs
            rss <- sum((z*w - qfitted)^2)
            rss/rdf # q cols have sums of squares = to 1, so division by that is implicit
        } else 1
    ## because the Q matrix is orthogonal, corresponding nominal Cov-hat is dispersion * Identity
    ans <- sqrt(2 * dispersion *sum(diag(Sqperp)))
    if (return.extras)
        {
    attr(ans, "dispersion") <- dispersion
    attr(ans, "scaled.variances") <- (nobs-1)*(mean(w)^-2)*
        cbind(Winv.dot.Q=diag(covqtilde), Winv.dot.Q.perp= diag(Sqperp))
}
    ans
}

ppse.array <- function(object, covariance.extractor=NULL, data=NULL,...)
  {
    stopifnot(inherits(object, "array"), length(dim(object))==3,
              !is.null(dimnames(object)), "cov.beta" %in% dimnames(object)[[3]],
              "Sperp" %in% dimnames(object)[[3]])
covb <- object[,,"cov.beta", drop=TRUE]
Sperp <- object[,,"Sperp", drop=TRUE]
sqrt(2 * sum(covb * Sperp))
  }
