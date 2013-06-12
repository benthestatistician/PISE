##' <description>
##' Extract and invert observed information from a logistic regression fit.
##' <details>
##' @title Retrieve inverse of observed information
##' @param PSmodel A logistic regression fit produced with \code{glm(...,family=binomial)}
##' @return  Inverse of observed information, 1/n scaling factor omitted 
##' @author Ben Hansen
GetInvObservedInfo <- function(PSmodel)
  {
stopifnot(PSmodel$rank>0, PSmodel$family$family=="binomial")
p <- PSmodel$rank
p1 <- 1L:p
Qr <- PSmodel$qr
coef.p <- PSmodel$coefficients[Qr$pivot[p1]]
covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
covmat.unscaled
  }

##' <description>
##' For a fitted logistic regression model, calculate the B matrix of the Huber-White variance estimator for the beta hats, with optional conditioning on stratum totals.
##' <details>
##' The function returns the unscaled version of the B matrix, i.e. division by n is omitted.  Somewhat oddly, passing a stratum argument leads the function to estimate this matrix in a way appropriate for estimating conditional variance of a the unconditional logistic regression estimates in \code{PSmodel}.
##' @title For logistic glm's, sum of outer products of score statistics with themselves
##' @param PSmodel A logistic regression fit produced with \code{glm(...,family=binomial)}
##' @param strata Factor variable, ordinarily the result of a call to \code{optmatch::fullmatch} or \code{optmatch::pairmatch}
##' @return Empirical B matrix minus 1/n scaling, i.e. sum of outer products of score statistics with themselves. 
##' @author Ben Hansen
EmpiricalScoreCov <- function(PSmodel, strata=NULL)
  {
    stopifnot(inherits(PSmodel, "glm"), PSmodel$family$family=="binomial",
              is.null(strata) || inherits(strata, "factor")
              )
    z.minus.pi <- residuals(PSmodel, type="response")
    
    thecoefnames <- names(coef(PSmodel))
    mm <- model.matrix(PSmodel)
    stopifnot(nrow(mm)==length(z.minus.pi), ncol(mm)==length(thecoefnames),
              isTRUE(all.equal(thecoefnames, colnames(mm)))
              )
    sumscores <- z.minus.pi*mm

    if (!is.null(strata)) {
      noNA <- !is.na(strata)
      sumscores <- sumscores[noNA,]
      stratMM <- SparseMMFromFactor(strata[noNA])
      sumscores <- t(stratMM) %*% sumscores
      sumscores <- as.matrix(sumscores)
      colnames(sumscores) <- thecoefnames
    }
      
      crossprod(sumscores)
  }

##' <description>
##'  Sandwich estimate of covariance matrix of the coefficients of a logistic regression model, optionally conditioning on the margins of a stratification or post-stratification.
##' <details>
##' Somewhat oddly, passing a stratum argument leads the function to estimate this matrix in a way appropriate for estimating conditional variance of a the unconditional logistic regression estimates in \code{PSmodel}.
##' @title Sandwich estimate of covariance matrix of the coefficients of a logistic regression
##' @param PSmodel 
##' @param strata 
##' @param inv.obs.info 
##' @return Sandwich estimate of covariance
##' @author Ben Hansen
SandwichVcov <- function(PSmodel,
                         strata=NULL,
                         inv.obs.info=GetInvObservedInfo(PSmodel))
  {
    meat <- EmpiricalScoreCov(PSmodel, strata=strata)
    stopifnot(all(colnames(inv.obs.info) %in% rownames(meat)))
    meat <- meat[colnames(inv.obs.info), colnames(inv.obs.info)]
    inv.obs.info %*% meat %*% inv.obs.info
  }

##' <description>
##'
##' <details>
##'  The expectation and variance are calculated with conditioning on matched set margins (\mathcal{Z}).  Justifications are in HW notes labeled "mid-June 2011", the framing of which differs somewhat from my current (April 2012) preferences -- an update might help.  Not sure the extent to which these results depend on the PS parameters being asymptotically Normal; particularly the variance of summed variances calculation seems most likely to have that dependence.
##' @title Matched variation in estimated propensity score, with permutation-based benchmarks
##' @param strata 
##' @param PSmodel 
##' @param inv.obs.info 
##' @return List containing: \code{sigmahatsq}, sum of matched set variances of estimated PS; \code{E_0} and \code{V_0}, estimated expectation and variance of sum of matched set variances in (PS-EPS); and \code{z}, a corresponding t-type statistic
##' @author Ben Hansen
SumSqPShatresids <- function(strata, PSmodel,
                             inv.obs.info=GetInvObservedInfo(PSmodel))
{
require("SparseM")
covbeta <- SandwichVcov(PSmodel=PSmodel, strata=strata,
                        inv.obs.info=inv.obs.info)

subset <- if (inherits(strata, "optmatch")) {
  optmatch::matched(strata)
} else !is.na(strata)
mm <- model.matrix(PSmodel)[subset, colnames(covbeta)]
thetahat <- predict(PSmodel)[subset]

strata.csr <- SparseMMFromFactor(strata[subset])
pstilde <- resid(SparseM::slm.fit(strata.csr, thetahat))
xtilde <- resid(SparseM::slm.fit(strata.csr, mm))
calS <- crossprod(xtilde)
calS.sqrt <- chol(calS, pivot=TRUE)
calS.sqrt <- calS.sqrt[, order(attr(calS.sqrt, "pivot"))]
covA <- calS.sqrt %*% covbeta %*% t(calS.sqrt)
eigvals <- eigen(covA, symmetric=TRUE, only.values = TRUE)$values
ans <- list(sigmahatsq=sum(pstilde^2), `E_0`=sum(eigvals), `V_0`=2*sum(eigvals^2))
ans$z <- with(ans, (sigmahatsq - E_0)/sqrt(V_0))
unlist(ans)
}
