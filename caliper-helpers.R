##' <description>
##' Standard error of the difference in propensity scores for a treatment-control pair with covariate differences set at their expected values
##' <details>
##' Given a fitted propensity model and two values of the corresponding covariate vector, a standard error attaches to the difference  of the two observations' propensity scores (on the linear predictor scale) in a straightforward way. This function sets the difference between those two vectors to be equal to the square root of the expected value of the square of said difference, under random sampling of one treatment and one control.  Well, almost: with \code{center.} and \code{scale.} arguments at their default values, the means and sd's of treatment and control groups are estimated in a resistant fashion, in anticipation of the likelihood of the matching process removing at least a small portion of each of those groups.
##' Another good choice for \code{covariance.extractor} is \code{sandwich}, from package \sQuote{sandwich}.
##' The default \code{scale.} function is modeled on \code{mad}, but arranged to permit \code{center.} to decide the centers of the covariate distributions.
##'  By default, output is scaled by reciprocal of pooled sd in the estimated propensity, for consistency with established ways of specifying propensity calipers.
##' @title std_psdiff_se
##' @param model A fitted propensity score model
##' @param covariance.extractor Function to extract covariance of fitted propensity coefficients.
##' @param scale. Function with which to estimate standard deviations of covariates.
##' @param center. Function with which to estimate centers of covariate distributions.
##' @param inverse.sd.scaling Divide result by the pooled sd of estimated propensity scores?
##' @param ... Additional arguments needed by \code{model.frame} or \code{model.matrix} (e.g. a data argument)
##' @return A number representing the covariance of the difference in estimated propensity scores between two units selected at random from the sample contributing to estimation of same propensity score.
##' @author Ben Hansen
std_psdiff_se <- function(model,covariance.extractor=vcov,inverse.sd.scaling=TRUE,scale.=function(x) {1.4826*median(abs(x))},center.=function(x){mean.default(x, trim=.025)},...) UseMethod("std_psdiff_se")
std_psdiff_se.default <- function(model,covariance.extractor=vcov,inverse.sd.scaling=TRUE,scale.=function(x) {1.4826*median(abs(x))},center.=function(x){mean.default(x, trim=.025)},...)
  {
    thez <- model.response(model.frame(model,...))
    thez <- as.logical(thez)
    if (any(is.na(thez)) || all(thez) || all(!thez)) stop("Based on its response vector, model\n doesn't appear to be a fitted propensity score.")
    
    thecov <- covariance.extractor(model)
    vnms <- colnames(thecov)[colnames(thecov)!="(Intercept)"]
    thecov <- thecov[vnms,vnms]
    
    thexes <- model.matrix(model)
    thexes <- thexes[,vnms]

    GetCenterAndScale <- function(matrix.)
      {
        centers <- apply(matrix.,2,center.)
        matrix. <- sweep(matrix.,2,centers)
        scales <- apply(matrix.,2, scale.)
        data.frame(centers=centers, scales=scales)
      }

    csc.tx <- GetCenterAndScale(thexes[thez,])
    csc.ctl <- GetCenterAndScale(thexes[!thez,])
    typical.seps <- (csc.tx$centers - csc.ctl$centers)^2 +
      csc.tx$scales^2 + csc.ctl$scales^2
    typical.seps <- sqrt(typical.seps)
    
    caliper.sq <- crossprod(typical.seps, thecov) %*% typical.seps
    ans <- sqrt(drop(caliper.sq))

    scl <- if (inverse.sd.scaling)
      {
        cf <- coef(model)[vnms]
        theps <- thexes %*% cf
        optmatch:::szn.scale(theps, thez)
      } else 1

    ans/scl
    }
