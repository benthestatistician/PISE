##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title SE of propensity-paired differences on a fitted propensity score
##' @param object fitted propensity score model, of or inheriting from class \code{glm}
##' @param covariance.estimator which of \code{vcov}, \code{sandwich} to use to estimate model coefficient covariances?
##' @param simplify return a scalar (\code{simplify==TRUE}, the default) or a list of quantities combine to make that scalar?
##' @param data a data frame
##' @return scalar, in units of \code{object}'s linear predictor, if simplify=T; otherwise a list
##' @author Mark M. Fredrickson, Ben B. Hansen
ppse <- function(object, covariance.estimator, data, simplify, ...)
    UseMethod("ppse")

ppse.glm <- function(object, covariance.estimator=c("vcov", "sandwich")[1],
                     data=NULL, simplify=TRUE, coeffs.from.fitted.model=FALSE, ...) 
{
    ppse_via_qr(object, covariance.estimator=covariance.estimator,
                data=data, simplify=simplify,
                coeffs.from.fitted.model=coeffs.from.fitted.model,...)
}

ppse.bayesglm <- function(object, covariance.estimator=c("vcov", "sandwich")[1],
                          data=NULL, simplify=TRUE, coeffs.from.fitted.model=TRUE, ...) 
{
    ppse_via_qr(object, covariance.estimator=covariance.estimator,
                data=data, simplify=simplify,
                coeffs.from.fitted.model=coeffs.from.fitted.model,...)
}



##'
##' .. content for \details{} ..
##' @title SE of propensity-paired differences on a fitted propensity score: version -1
##' @param object as in \code{ppse}
##' @param covariance.estimator as in \code{ppse}
##' @param data as in \code{ppse}
##' @param tt as in \code{ppse}
##' @param simplify as in \code{ppse}
##' @param terms.to.sweep.out terms in calling formula that should be swept out
##' @param ... 
##' @return 
##' @author Mark M. Fredrickson, Ben B Hansen
ppse_notstabilized <- function(object, covariance.estimator="vcov",
                         data=model.frame(object), simplify=TRUE,
                         tt=terms(formula(object), specials="strata"),
                         terms.to.sweep.out=survival:::untangle.specials(tt, "strata")$terms,
                         cluster=NULL,...)
    {
        if (is.null(names(coef(object)))) stop("propensity coefficients have to have names")

        data.matrix <- model.matrix(tt, data)
        
        stopifnot(!is.null(colnames(data.matrix)),
                  setequal(colnames(data.matrix), names(coef(object))),
                  isTRUE(all.equal(colnames(data.matrix), names(coef(object)))),
                  !is.null(attr(data.matrix, "assign")),
                  covariance.estimator %in% c("vcov", "sandwich"),
                  covariance.estimator=="sandwich" | is.null(cluster),
                  !is.null(cluster) || length(cluster)!=nrow(data)
                  )
        if (!is.null(cluster)) stop("cluster arg not currently supported")
        if (covariance.estimator=="sandwich") stopifnot(require("sandwich"))
        covariance.extractor <- eval(parse(text=covariance.estimator))
        
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

        ans <- list("cov.betahat"= covb, "betahat"=coeffs, "cov.X"=S)
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
ppse_via_qr <-
    function(object, covariance.estimator=c("vcov", "sandwich")[1],
             data=NULL, simplify=TRUE,
             coeffs.from.fitted.model=FALSE, QR=object$qr, 
             tt=terms(formula(object),specials="strata"),
             terms.to.sweep.out=survival:::untangle.specials(tt, "strata")$terms,
             tol.coeff.alignment=Inf,
             cluster=NULL, ...) 
{
    
    stopifnot(inherits(object, "glm"),
              !is.null(QR), 
              length(QR$pivot)<=length(coef(object)),
              as.logical(QR$rank),
              covariance.estimator %in% c("vcov", "sandwich"),
              covariance.estimator=="sandwich" | is.null(cluster), 
              is.null(terms.to.sweep.out) | 
              is.numeric(terms.to.sweep.out) &
              all(as.integer(terms.to.sweep.out)==terms.to.sweep.out))
    terms.to.sweep.out <- unique(terms.to.sweep.out)
    nterms <-  length(attr(tt, "term.labels"))
    if (any(terms.to.sweep.out > nterms)) stop("terms.to.sweep.out values too big")
    if (!is.null(terms.to.sweep.out) && length(terms.to.sweep.out))
    {
        first_non_sweep <- min(setdiff(1L:nterms, terms.to.sweep.out))
        dontsweep <- (terms.to.sweep.out >= first_non_sweep)
        if (first_non_sweep < max(terms.to.sweep.out))
        warning(paste0(sum(dontsweep), 
	" strata() directive(s) ignored. To fix, put at beginning of model formula."))

        terms.to.sweep.out <- terms.to.sweep.out[!dontsweep]
        }


    glm.family.uses.estimated.dispersion <-
        !(substr(object$family$family, 1, 17) %in%  # borrowed from sandwich:::bread.glm
          c("poisson", "binomial", "Negative Binomial"))
    fitted.model.coeffs <- coef(object)
    fitted.model.coeffs.NA <- is.na(fitted.model.coeffs)
    fitted.model.coeffs[fitted.model.coeffs.NA] <- 0
    if (is.null(data)) data <- model.frame(object)
    data.matrix <- model.matrix(tt, data)
    Xcols_to_terms <- attr(data.matrix, "assign")
    Xcols_to_sweep_out <- Xcols_to_terms  %in% c(0,terms.to.sweep.out)
    Qcols_to_keep <- seq_len(QR$rank)
    if (any(Xcols_to_sweep_out))
        {
    Xcols_to_sweep_out <- names(fitted.model.coeffs)[Xcols_to_sweep_out]
    Qcols_to_sweep_out <- match(Xcols_to_sweep_out, colnames(qr.R(QR)))
    Qcols_to_sweep_out <-
        Qcols_to_sweep_out[Qcols_to_sweep_out==seq_along(Qcols_to_sweep_out)]
    Qcols_to_keep <- setdiff(Qcols_to_keep, Qcols_to_sweep_out)
        }

    ## Now we start a bunch of calculations that re-figure 
    ## internal components of the glm fit: z, eta, weights. 
    stopifnot(is.null(model.offset(data))) # assume away offsets (for now!)    
    offset <- rep(0, nrow(data))
    eta <- object$linear.predictors 

    
    ## refigure weights, since glm.fit doesn't update them after last iteration
    weights <- getglmQweights(eta, 
                              prior.weights=object$prior.weights, 
                              family=object$family)
    stopifnot(length(weights)==nrow(data))
    ## NB: To address discrepancies between coeffs reported w/ model fit
    ## and coeffs constructed from the returned QR decomposition, I tried:
    ## (1) using the weights based on penultimate glm fit,
    ## not the updated weights as figured by my `getglmQweights`.
###    weights <- object$weights
    ## (2) insisting on convergence,
###    stopifnot(object$converged)
    ## Rationale for requiring convergence: Since we don't also have access to
    ## penultimate eta's, we're relying on the convergence being far enough along
    ## that their differences  from the final etas are numerically negligible.
    ## I get the impression that the `glm.fit` convergence criteria are such as to
    ## ensure this.
    ## Neither did what I was hoping for: with the aglm example used in optmatch,
    ## I got differences  in coefficients of order 2e-6 either way. Chalking the
    ## discrepancy up (uneasily) to differences between numerical calcs used within
    ## `C_Cdqrls` and the QR-based method I'm using. 


    w <- sqrt(weights) 
    good <- !is.na(w) & w>0

    data.matrix <- data.matrix[good,]
    eta <- eta[good]
    w <- w[good]
    if (!is.null(cluster)) cluster <- factor(cluster[good])
    ## NB: successful `getglmQweights` confirms that `linkinv` is there and is a function
    linkinv <- object$family$linkinv
    mu <- linkinv(eta)
    y <- model.response(data, type="double")
    twocol.response <- !is.null(dim(y))
    y <- y[good]
    mu.eta <- object$family$mu.eta
    mu.etaval <- mu.eta(eta)
    resids <- if (twocol.response)
      {
        object$residuals[good]
    }
    else (y-mu)/mu.etaval # Roll our own if it's easy enough
    
    z <- (eta -offset) + resids # "z" as in `glm.fit`
    ## End reconstruction of internals of glm-fitting

    qmat <- qr.Q(QR)[good,]
    colnames(qmat) <- qnames <- paste0("Q.",colnames(qr.R(QR)))


    ## NB: `qcoeffs <- qr.qty(QR,z*w)` doesn't work, returns object of length length(z),
    ## not QR$rank.  Likewise for qr.qy(...)

    if (coeffs.from.fitted.model || is.finite(tol.coeff.alignment))
      {
        qcoeffs.from.fitted.model <- drop(qr.R(QR)%*%fitted.model.coeffs[QR$pivot])
        names(qcoeffs.from.fitted.model) <- qnames
        qcoeffs.from.QR <- NULL
        qcoeffs <- qcoeffs.from.fitted.model
      } else qcoeffs.from.QR <- crossprod(qmat, z*w)
    if (!coeffs.from.fitted.model) {
        if (is.null(qcoeffs.from.QR))
                qcoeffs.from.QR <- crossprod(qmat, z*w)
        qcoeffs <- qcoeffs.from.QR
                   }

    ## this is here for testing purposes 
    if (is.finite(tol.coeff.alignment))
    {
        if (is.null(qcoeffs.from.QR)) qcoeffs.from.QR <- crossprod(qmat, z*w)
            qcoeffs.from.QR[!Qcols_to_keep] <- 0
            coeff.diffs <- qcoeffs.from.fitted.model - qcoeffs.from.QR
            if (max(abs(coeff.diffs)) >= tol.coeff.alignment)
            stop(paste("QR/reported coefficients differ by up to", prettyNum(max(abs(coeff.diffs)))))
        }

    qcoeffs <- qcoeffs[Qcols_to_keep]
    qmat <- qmat[,Qcols_to_keep]
    
    ##  xtilde is the rotation of X that expresses it in the same coordinates 
    ##  as the Q matrix.  I.e., if diag(w)X = QR, Xtilde = X %*% R^(-1)
    ## First I was calculating this as diag(w^(-1))%*%Q , i.e.
    ## xtilde <- w^(-1) * qmat
    ## but I decided against this when I saw that in easy examples
    ## (cf "aglm" in tests) it was creating an (Intercept) column w/ nonzero variance.
    ## The below instead calculates  X %*% R^(-1) directly.
    ## This requires a little more computing but seemed to avoid the problem w/ the
    ## (Intercept) column.
    xtilde <- t(backsolve(qr.R(QR),t(data.matrix),
                          k=QR$rank, transpose=T)
                )
    xtilde <- xtilde[,Qcols_to_keep]
    
    covxtilde <- cov(xtilde)

    ## On to estimating (co)variance of the qcoeffs...
    nobs <- sum(good)
    stopifnot((rdf <- nobs - QR$rank)>0) # i.e., fitted.model$df.residual
    
    dispersion <-
        if (glm.family.uses.estimated.dispersion &&
            covariance.estimator=="vcov" #Won't use dispersion in sandwich calc
            )
        {              
            ## Following summary.glm
            df.r <- object$df.residual
            if (df.r>0)
                sum((resids^2* weights)[weights>0])/df.r else NaN
        } else 1
    
    ans <- if (covariance.estimator=="vcov")
               { ## because the Q matrix is orthogonal, corresponding
                 ## nominal Cov-hat is dispersion * Identity
                   list("cov.betahat"=dispersion, "betahat"=qcoeffs, "cov.X"=covxtilde)
               } else { # in this case covariance.estimator=="sandwich"
                   esteqns <- # this calc should be the same as 
                      qmat * (resids * w) #  xtilde * weights * resids
                   if (!is.null(cluster)) {
                       esteqns <- aggregate(esteqns, by = list(cluster), FUN = sum)[,-1]
                       esteqns <- as.matrix(esteqns)
                   }
                   meatmatrix.unscaled <- # (unscaled bread being the identity)
                        crossprod(esteqns)
                   ## Per KISS principle, no d.f. adjustments. For now. 
                   list("cov.betahat"=meatmatrix.unscaled,
                        "betahat"=qcoeffs, "cov.X"=covxtilde) 
           }
    if (simplify) ppse(ans,...) else ans
}
##' Convert separated ppse calculations into a scalar ppse
##'
##' Calculating the ppse involves estimation of a matrix
##' matrix of regression coefficients and calculation of
##' corresponding data covariance, both done with respect to
##' the same coordinates. In the process, estimates of those
##' regression coefficients are also obtained, in the same
##' coordinate system. This function converts a list with these
##' components into a ppse. 
##' and of a 
##' @title 
##' @param object 
##' @param covariance.estimator 
##' @param data 
##' @param ... 
##' @return 
##' @author Ben B Hansen
ppse.list <- function(object, covariance.estimator=NULL, data=NULL,...)
  {
    stopifnot(all(!is.na(pmatch(c("cov.betahat","betahat", "cov.X"),names(object)))),
              is.numeric(object$cov.betahat),
              is.numeric(object$betahat),
              is.numeric(object$cov.X),
              is.matrix(object$cov.X),
              length(object$betahat)==nrow(object$cov.X),
              length(object$cov.betahat)==1 ||
              length(object$cov.betahat)==length(object$cov.X))

    cov_beta <- if (length(object$cov.betahat)==1) { # this really means dispersion * Identity
                    object$cov.betahat * diag(nrow(object$cov.X))
                } else object$cov.betahat
    
    Sperp <-  makeSperp(object$cov.X, object$betahat)

    ## calculate the correction for the expected difference in propensity scores
    sqrt(2 * sum(cov_beta * Sperp))
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
    
##    resids <- object$residuals    
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

     
    oldrank.adjusted <- oldrank - sum(match(cols.to.sweep.out, cols.R, nomatch=(oldrank+1)) <= oldrank)

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
    ans <- ppse_via_QR(theglm, QR=newqr, simplify=FALSE)
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

    ##  xtilde is the rotation of X that expresses it in the same coordinates as the Q matrix.
    ##  I.e., if diag(w)X = QR, write Xtilde for  X %*% R^(-1)
    xtilde <- t(backsolve(qr.R(QR),t(data.matrix),
                          k=QR$rank, transpose=T)
                )


## qcols.to.keep <- seq_len(newqr$rank)
##eta.from.q <- xtilde[,qcols.to.keep] %*% qcoeffs[qcols.to.keep]

ans$delta.ps.sq.over.ppse.sq <- sum((xtilde[,newqr$rank] * qcoeffs[newqr$rank])^2)/
  ((sum(good)-1) *sum(ans[["Sperp.diagonal"]]))
ans
  }
