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
