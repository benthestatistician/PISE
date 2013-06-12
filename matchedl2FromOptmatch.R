##' <description>
##'
##' <details>
##' @title Matched covariances of vectors or matrices
##' @param object Numeric vector or matrix, or fitted from which a design matrix is to be extracted
##' @param strata an optmatch object
##' @param ...  Additional arguments used by methods (currently ignored)
##' @return Scalar or matrix, representing a within-matched-set variance or covariance, resp, pooled across matched sets
##' @author Ben Hansen
matchedl2sq <- function(object, strata, ...) UseMethod("matchedl2sq")
matchedl2sq.default <- function(object, strata,...)
  {
stopifnot(inherits(strata, "optmatch"),
          is.null(dim(object)))
wts <- getl2weights(strata)
res <- numeric(length(wts))
names(res) <- names(strata)
res[!is.na(strata)] <- residuals(slm(object~strata, na.action=na.omit))
sum(res^2 *wts)
  }
matchedl2sq.matrix <- function(object, strata,...)
  {
stopifnot(inherits(strata, "optmatch"))
wts <- getl2weights(strata)
res <-   matrix(0,
                nrow(object), ncol(object),
                dimnames=dimnames(object))
res[!is.na(strata),] <- residuals(slm(object~strata, na.action=na.omit))
crossprod(res*wts, res)
  }
matchedl2sq.lm <- function(object, strata,...)
  {
    mat <- model.matrix(object)
    mat <- mat[,colnames(mat)!="(Intercept)"]
    matchedl2sq(mat,strata=strata,...)
  }
getl2weights <- function(strata)
  {
stopifnot(inherits(strata, "optmatch"))
tx <- as.numeric(attr(strata, "contrast.group"))
txfrac <- numeric(length(tx))
txfrac[!is.na(strata)] <- fitted(slm(tx ~ strata, na.action=na.omit))
ms <- table(strata)
n <- 2*sum((ms-1)/ms) # shortcut to sum of harmonic means of cardinalities of intersections of matched sets with tx and with control groups 
pmin(txfrac, 1-txfrac)/n
  }

##' <description>
##'
##' <details>
##' To compute standard errors be computed by Huber-White and with conditioning on treatment margins of matched sets, set \code{covariance.extractor} to \code{SandwichVcov}, as defined in file sandwich_glm.R.  \code{covariance.extractor} could also be set to \code{vcov} or \code{sandwich}, from package \sQuote{sandwich}.
##' For a brief justification, see hwn2012-04-06
##' @title  Expected value of mean square difference between true and estimated propensity scores
##' @param amatch  An optmatch object, ie the value of a call to fullmatch() or pairmatch()
##' @param PSmodel A fitted propensity score model, e.g. object of class glm
##' @param covariance.extractor Function used to extract covariance estimate from \code{PSmodel}.
##' @return A scalar, the estimated expected value of \code{crossprod(x%*%(betahat - beta))}.
##' @author Ben Hansen
Exp_matched_l2sq_EPSminusPS <- function(amatch, PSmodel, covariance.extractor=vcov)
  {
    stopifnot(inherits(amatch, "optmatch"),
              all(names(amatch) %in% row.names(mf <- model.frame(PSmodel)))
              )
    thecoef <- coef(PSmodel)
    zz <- as.logical(model.response(mf))
    mm <- model.matrix(PSmodel)
    stopifnot(length(notanintercept <- attr(mm, "assign")!=0L)==length(thecoef),
              isTRUE(all.equal(names(thecoef), colnames(mm)))
              )
    if (any(isNA <- is.na(thecoef)))
      {
        thecoef <- thecoef[!isNA]
        mm <- mm[,!isNA]
        notanintercept <- notanintercept[!isNA]
      }
    thecoef <- thecoef[notanintercept]
    mm <- mm[,notanintercept]
    if (!isTRUE(all.equal(names(thecoef), colnames(mm))))
        stop("Something is wrong -- please debug!")
    coefnames <- names(thecoef)
    thecov <- covariance.extractor(PSmodel, strata=amatch)
    thecov <- thecov[coefnames, coefnames]
    theml2 <- matchedl2sq(mm, amatch)
    theml2 <- theml2[colnames(thecov),colnames(thecov)]
    sum(theml2*thecov)
  }
