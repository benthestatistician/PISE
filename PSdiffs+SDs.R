##' .. content for \description{} (no empty lines) ..
##' Function to extract from fitted propensity score models treatment-by-control
##' matrices of: (a) separations on the linear propensity score; and
##' (b) corresponding standard errors.
##' .. content for \details{} ..
##' @title Propensity score difference maker
##' @param x A fitted propensity score model
##' @param subset subset of the observations to use
##' @param covariance.extractor function to extract covariance of fitted model.
##' @param usedimnames logical: return matrices with dimnames?
##' @return list: \sQuote{differences}, a matrix of PS differences; \sQuote{SEs}, matrix of corresp. standard errors.
##' @author Ben Hansen
getPSdiffs <- function(x, covariance.extractor=vcov, subset=NULL, usedimnames=TRUE)
  {
    stopifnot(is.null(subset) || is.logical(subset),
              length(subset)==nrow(mf <- model.frame(x)) || is.null(subset))
    if (is.null(subset)) subset <- rep(TRUE, nrow(mf))
    thecoef <- coef(x)
    thecov <- covariance.extractor(x)
    zz <- as.logical(model.response(mf))
    mm <- model.matrix(x)
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
    thecov <- thecov[notanintercept,notanintercept]
    mm.t <- mm[zz&subset,notanintercept]
    mm.c <- mm[!zz&subset,notanintercept]
    n.t <- nrow(mm.t)
    n.c <- nrow(mm.c)

    if (!isTRUE(all.equal(names(thecoef), colnames(mm.t)))||
        !isTRUE(all.equal(names(thecoef), colnames(mm.c))) )
        stop("Something is wrong in getPSdiffs().  Please debug!")

    diffandSE <- function(vec)
      {
      thediffs <- sweep(mm.c, 2, vec)
      thePSdiffs <- -thediffs%*%thecoef
      thePSdiffSEs <- rowSums((thediffs%*%thecov)*thediffs)
      thePSdiffSEs <- sqrt(thePSdiffSEs)
      c(thePSdiffs,thePSdiffSEs)
    }

    y <- apply(mm.t, 1, diffandSE)
    diffMat <- t(y[1L:n.c,])
    SEMat <- t(y[n.c+(1L:n.c),])
    if (usedimnames)
      {
        thern <- row.names(mf)
        dimnames(diffMat) <-
          dimnames(SEMat) <-
            list(Tx=thern[zz&subset],
                 Ctl=thern[!zz&subset])
      }
    list(differences=diffMat, SEs=SEMat)
  }

##' .. content for \description{} (no empty lines) ..
##'  Function to assemble SEs of propensity score differences as an
##' \code{optmatch.dlist} object, ie in a form useful for matching on them
##' or creating calipers in terms of them.
##' .. content for \details{} ..
##' @title SEs of propensity score differences, as a matching \sQuote{distance}
##' @param x A fitted propensity score model
##' @param structure.fmla Optional formula indicating subclasses within which to match
##' @param covariance.extractor Function to extract covariance from propensity model
##' @return SEs of (linear) PS differences, as an \code{optmatch.dlist}.
##' @author Ben Hansen
PSSEdist <- function(x, structure.fmla = NULL, covariance.extractor=vcov)
  {
    stopifnot(inherits(x, "glm"),
              all(substr(names(mf <- model.frame(x)),1,3)!="._.")
              )

    if (is.null(structure.fmla)) structure.fmla <- update.formula(x, ~1)
    thecoef <- coef(x)
    thecov <- covariance.extractor(x)
    zz <- as.logical(model.response(mf))
    mm <- model.matrix(x)
    stopifnot(length(notanintercept <- attr(mm, "assign")!=0L)==length(thecoef),
              isTRUE(all.equal(names(thecoef), colnames(mm)))
              )
    if (any(isNA <- is.na(thecoef)))
      {
        thecoef <- thecoef[!isNA]
        mm <- mm[,!isNA]
        notanintercept <- notanintercept[!isNA]
      }
    mm <- mm[,notanintercept]
    thecoef <- thecoef[notanintercept]
    thecov <- thecov[notanintercept,notanintercept]
    mdnms <- paste("._.",1L:ncol(mm),sep="")
    colnames(mm) <- mdnms
    ndat <- data.frame(mf, `.ZzZ`=zz, mm)
    strfmla <- update.formula(structure.fmla, .ZzZ~.)
    dimnames(thecov) <- list(mdnms,mdnms)

PSSEmat <- function(treatments,controls, thecov=thecov,mdnms=mdnms)
  {
    PSSErow <- function(vec, mat, thecov)
      {
        thediffs <- as.matrix(sweep(mat, 2, vec))
        theSEs <- rowSums((thediffs%*%thecov)*thediffs)
        theSEs <- sqrt(theSEs)
        theSEs
      }
    dmat <- apply(treatments[mdnms],1,PSSErow,mat=controls[mdnms],thecov=thecov)
    t(dmat)
  }

mdist(PSSEmat, structure.fmla=strfmla, data=ndat, thecov=thecov, mdnms=mdnms)
  }
