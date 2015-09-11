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

ppse.default <- function(object, covariance.extractor=vcov, data=model.frame(object), tt=terms(object), simplify=TRUE, ...)
    {
        if (is.null(names(coef(object)))) stop("propensity coefficients have to have names")

        data.matrix <- model.matrix(tt, data)
        
        stopifnot(!is.null(colnames(data.matrix)),
                  setequal(colnames(data.matrix), names(coef(object))),
                  isTRUE(all.equal(colnames(data.matrix), names(coef(object)))),
                  !is.null(attr(data.matrix, "assign"))
                  )

        coeffs <- coef(object)
        coeffnames <- names(coeffs)
        coeff.NA <- !is.finite(coef(object))
        
        vnames <- coeffnames[!(attr(data.matrix, "assign")==0 | #exclude intercept
                           coeff.NA)]
        coeffs <- coeffs[vnames]

        
        stopifnot(!is.null(dimnames(covb <- covariance.extractor(object))),
                  setequal(dimnames(covb)[[1]], coeffnames[!coeff.NA]))
        
        if (is.null(dimnames(covb)))
          {
            dimnames(covb) <- list(coeffnames[!coeff.NA], coeffnames[!coeff.NA])
          }
        covb <- covb[vnames, vnames]

        covx <- cov(data.matrix)
        covx <- covx[vnames, vnames]

        
        terms.to.sweep.out <- survival:::untangle.specials(tt, "strata")$terms ## strata, if present
        cols.to.sweep.out <- attr(data.matrix, "assign") %in% c(0, terms.to.sweep.out)
    cols.to.sweep.out <- colnames(data.matrix)[cols.to.sweep.out]

        cols.to.keep <- !(vnames %in% cols.to.sweep.out)
        
        covb <- covb[cols.to.keep, cols.to.keep, drop=FALSE]
        coeffs <- coeffs[cols.to.keep]
        
        
        S <- if (sum(cols.to.keep)!=length(vnames))
        {
          n <- nrow(data.matrix)
          K <- length(vnames) - sum(cols.to.keep)
          
            (covx[cols.to.keep, cols.to.keep] -
                covx[cols.to.keep,!cols.to.keep, drop=FALSE] %*%
                    solve(covx[!cols.to.keep, !cols.to.keep, drop=FALSE],
                          covx[!cols.to.keep, cols.to.keep, drop=FALSE])
             ) * ((n-1)/(n-K))
        } else covx

        Sperp <- makeSperp(S, betas=coeffs)

        ans <- list("cov.beta"= covb, "Sperp"=Sperp)
        if (simplify) ppse(ans,...) else ans
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
                    simplify=TRUE,...) 
{
    stopifnot(inherits(fitted.model, "glm"),
              length(object$pivot)<=length(coef(fitted.model)),
              as.logical(object$rank))
    if (!identical(covariance.extractor,vcov)) stop('ppse.qr only supports vcov as covariance extractor')

    glm.family.uses.estimated.dispersion <-
        !(substr(fitted.model$family$family, 1, 17) %in%  # borrowed from sandwich:::bread.glm
          c("poisson", "binomial", "Negative Binomial"))
    tt <- terms(fitted.model)
    fitted.model.coeffs <- coef(fitted.model)
    fitted.model.coeffs.NA <- is.na(fitted.model.coeffs)
    fitted.model.coeffs[fitted.model.coeffs.NA] <- 0
##    stopifnot(!is.null(fitted.model$weights))
    ## b/c I'm not sure how this interacts w/ factor expansion in next few lines
    if (is.null(data)) data <- model.frame(fitted.model)
    stopifnot(is.null(model.offset(data))) # assume away offsets (for now!)    
    offset <- rep(0, nrow(data))
    eta <- fitted.model$linear.predictors 

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
##    resids <- fitted.model$residuals
    data.matrix <- data.matrix[good,]
##    resids <- resids[good]
    eta <- eta[good]
    w <- w[good]
    ## NB: successful `getglmQweights` confirms that `linkinv` is there and is a function
    linkinv <- fitted.model$family$linkinv
    mu <- linkinv(eta)
    y <- model.response(data, type="double")
    twocol.response <- !is.null(dim(y))
    y <- y[good]
    mu.eta <- fitted.model$family$mu.eta
    mu.etaval <- mu.eta(eta)
    resids <- if (twocol.response)
      {
        fitted.model$residuals[good]
    }
    else (y-mu)/mu.etaval # Roll our own if it's easy enough
    
    z <- (eta -offset) + resids # "z" as in `glm.fit`
    

## Next 2 lines seem to be idle...   
###    if (ncol(data.matrix) !=length(object$pivot)) stop("length of QR's pivot doesn't match ncol(data.matrix).")
###    data.matrix <- data.matrix[,object$pivot]
    qmat <- qr.Q(object)[good,]
    colnames(qmat) <- qnames <- paste0("Q.",colnames(qr.R(object)))


    ## NB: `qcoeffs <- qr.qty(object,z*w)` doesn't work, returns a object of length length(z),
    ## not object$rank.  Likewise for qr.qy(...)
    qcoeffs.from.QR <- crossprod(qmat, z*w)

    if (coeffs.from.fitted.model || is.finite(tol.coeff.alignment))
      {
        qcoeffs.from.fitted.model <- drop(qr.R(object)%*%fitted.model.coeffs[object$pivot])
        names(qcoeffs.from.fitted.model) <- qnames
      }
    qcoeffs <- if (coeffs.from.fitted.model) qcoeffs.from.fitted.model else qcoeffs.from.QR

    keep.these.Qcolumns <- seq_len(object$rank)
    qcoeffs <- qcoeffs[keep.these.Qcolumns]
    qmat <- qmat[,keep.these.Qcolumns]

    ## this is here for testing purposes 
    if (is.finite(tol.coeff.alignment))
        {
            qcoeffs.from.QR[!keep.these.Qcolumns] <- 0
            coeff.diffs <- qcoeffs.from.fitted.model - qcoeffs.from.QR
            if (max(abs(coeff.diffs)) >= tol.coeff.alignment)
            stop(paste("QR/reported coefficients differ by up to", prettyNum(max(abs(coeff.diffs)))))
        }

    
    ##  qtilde is the reweighting of the Q-matrix that corresponds to a rotated X. (Q ~ rotated W*X)    
    qtilde <- w^(-1) * qmat
    covqtilde <- cov(qtilde)

    Sqperp <-  makeSperp(covqtilde, qcoeffs)

    ## On to estimating (co)variance of the qcoeffs...
    nobs <- sum(good)
    stopifnot((rdf <- nobs - object$rank)>0) # i.e., fitted.model$df.residual
    
    dispersion <-
        if (glm.family.uses.estimated.dispersion)
        {              
            qfitted <- qmat %*% qcoeffs
            rss <- sum((z*w - qfitted)^2)
            rss/rdf # q cols have sums of squares = to 1, so division by that is implicit
            ## Instead follow: sandwich:::bread.glm?  They use
            ## wres <- as.vector(residuals(x, "working")) * weights(x, "working")
            ## sum(wres^2)/sum(weights(x, "working"))
        } else 1
    ## because the Q matrix is orthogonal, corresponding nominal Cov-hat is dispersion * Identity
    ans <- list("cov.beta"=dispersion, "Sperp.diagonal"=diag(Sqperp))
    if (simplify) ppse(ans,...) else ans
}

ppse.list <- function(object, covariance.extractor=NULL, data=NULL,...)
  {
    stopifnot(all(!is.na(pmatch(c("cov.beta","Sperp"),names(object)))),
              is.numeric(object$cov.beta),
              is.numeric(object$Sperp),
              length(object$cov.beta)==1 ||
              length(object$cov.beta)==length(object$Sperp))

    ## calculate the correction for the expected difference in propensity scores
    sqrt(2 * sum(object$cov.beta * object$Sperp))
  }

## this is tricky: for `glm`, R gets its QR decomposition from LINPACK,
## presumably b/c it likes the LINPACK approach to singular matrices.
## LINPACK does the minimal pivoting to remove the singularity, hewing closely
## to the ordering in which the RHS variables were supplied; but it does no
## additional pivoting.  (You can try to force LINPACK
## to remove more of them by mucking with the tolerance, but because it doesn't select
## the variable most deserving of elimination I find that the result can be terrible
## as an approximation to the original problem solution.) 
## LAPACK pivots so as to arrange an R matrix the
## diagonal entries of which are strictly decreasing in magnitude, which is nice, but
## it appears not to null out enough to avoid singularities.  
## By starting with
## ordinary R's LINPACK QR decomp we identify the columns that have to be remove, then
## we update to LAPACK for a follow-up solve in order to identify the next best
##  candidates for elimination.  

redo_qr  <- function(object, LAPACK=TRUE, tol=1e-07) #, precentering=FALSE
{
    stopifnot(inherits(object, "glm"),"qr" %in% names(object), is.qr(object$qr),
              as.logical(object$qr$rank), length(object$qr$pivot)>=object$qr$rank)
    force(LAPACK)
    force(tol)

    oldrank <- object$qr$rank
    
    tt <- terms(object, specials="strata")
    terms.to.sweep.out <- survival:::untangle.specials(tt, "strata")$terms ## strata, if present

    weights <- getglmQweights(object$linear.predictors, 
                              prior.weights=object$prior.weights, # after final iteration
                              family=object$family)
    data <- model.frame(object)
    stopifnot(length(weights)==nrow(data))
    w <- sqrt(weights) 
    good <- !is.na(w) & w>0
    data.matrix <- model.matrix(tt, data)

    cols.to.sweep.out <- attr(data.matrix, "assign") %in% c(0, terms.to.sweep.out)
    cols.to.sweep.out <- colnames(data.matrix)[cols.to.sweep.out]
    cols.to.keep <- setdiff(colnames(data.matrix), cols.to.sweep.out)
    
##    resids <- fitted.model$residuals    
    stopifnot(all(cols.to.keep %in% ( cols.R <- colnames(qr.R(object$qr)))),
              all(!duplicated(cols.R)))
    data.matrix <- data.matrix[good,]

    piv <- seq_len(oldrank)
    whichcol <- object$qr$pivot[piv]
    cols.to.zero.out <- cols.R[-whichcol]

    Xw <- w*data.matrix
    ## Now we need to solve for/sweep out the intercept & strata terms
    Xw.red <- lm.fit(Xw[,cols.to.sweep.out, drop=FALSE], Xw[,cols.to.keep])$residuals
    colnames(Xw.red) <- cols.to.keep
    ## ...or do we find the qr first, then use qr.solve?
    Xw.red[, cols.to.zero.out] <- 0
    ans <- qr(Xw.red, LAPACK=LAPACK, tol=tol)

     
    oldrank.adjusted <- oldrank - sum(match(cols.to.sweep.out, cols.R, nomatch=Inf) <= oldrank)

    ## If LAPACK=T, then rank is always returned as ncol(Xw). Override this.
    ## Note that the LAPACK answer is the one w/ the convenient property that cols
    ## we zeroed out beforehand will be pivoted to the very back, even in the presence of
    ## additional singularities -- see test.qr_trim.R.
    ans$rank <- min(ans$rank, oldrank.adjusted)
    ans
}
    
drop1_ppse_stats <- function(theglm, data=NULL)
  {
      stopifnot(theglm$qr$rank>=2)
      kappatri <- if (getRversion() <='2.15.1') kappa.tri else .kappa_tri
    newqr <- redo_qr(theglm)
    ans <- ppse(newqr, fitted.model=theglm, simplify=FALSE)
    ans$ppse <- ppse(ans)

    qmat <- qr.Q(newqr)
    qmat <- qmat[,seq_len(newqr$rank)]
    colnames(qmat) <- qnames <- paste0("Q.",colnames(qr.R(newqr)))

    ## condition number calc adapted from kappa.qr
    R.dropped <- R <- newqr$qr[1L:min(dim(newqr$qr)), seq_len(newqr$rank),drop=FALSE]
    R[lower.tri(R)] <- 0
    ans$kappa <- kappatri(R)

## Which column to drop? The rightmost one in R -- except don't drop intercept 
## column!  (If there are stratifying/exact matching vars then they should
## be excluded also -- not dealing with that just now)
    if (grepl("(Intercept)", names(ans[["Sperp.diagonal"]])[newqr$rank], fixed=F))
      {
        ans$skipped.intercept <- TRUE
        colshuffle <- newqr$rank - 0:1
        if (newqr$rank>2) colshuffle <- c(1L:(newqr$rank-2), colshuffle)
        qmat <- qmat[,colshuffle]
        R.dropped <- R.dropped[,colshuffle]
        R.dropped <- if (nrow(R.dropped)>=newqr$rank) R.dropped[colshuffle,]
        else R.dropped[colshuffle[1L:nrow(R.dropped)],]
      } else  ans$skipped.intercept <- FALSE

    R.dropped <- R.dropped[-newqr$rank,-newqr$rank]
    R.dropped[lower.tri(R.dropped)] <- 0
    ans$kappa.drop1 <- kappatri(R.dropped)



    if (is.null(data)) data <- model.frame(theglm)
    stopifnot(is.null(model.offset(data))) # assume away offsets (for now!)    
    offset <- rep(0, nrow(data))
    eta <- offset +
        theglm$linear.predictors # lazy, for now. 

## Now let's reconstruct the linear predictor, with and without that last Q-column

    w <- sqrt(theglm$weights)
    good <- !is.na(w) & w>0
    eta <- eta[good]
    linkinv <- theglm$family$linkinv
    mu <- linkinv(eta)
    y <- model.response(data, type="double")
    y <- y[good] 
    qmat <- qmat[good,]
    mu.eta <- theglm$family$mu.eta
    mu.etaval <- mu.eta(eta)
    resids <- (y-mu)/mu.etaval
    
    z <- (eta -offset) + resids # "z" as in `glm.fit`
    qcoeffs <- crossprod(qmat, z*w)

    ##  qtilde is the reweighting of the Q-matrix that corresponds to a rotated X. (Q ~ rotated W*X)    
    qtilde <- w^(-1) * qmat


## qcols.to.keep <- seq_len(newqr$rank)
##eta.from.q <- qtilde[,qcols.to.keep] %*% qcoeffs[qcols.to.keep]

ans$delta.ps.sq.over.ppse.sq <- sum((qtilde[,newqr$rank] * qcoeffs[newqr$rank])^2)/
  ((sum(good)-1) *sum(ans[["Sperp.diagonal"]]))
ans
  }
