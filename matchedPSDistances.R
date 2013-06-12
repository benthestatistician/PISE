makepooledsd <- function(x, Tx, standardization.scale) 
{
  stopifnot(is.numeric(x)||is.logical(x), is.logical(Tx))
  if (is.null(standardization.scale)) return(1)
  if (is.numeric(standardization.scale)) return(standardization.scale)
  if (is.function(standardization.scale)) 
    return(optmatch:::szn.scale(x=x, Tx=Tx, standardizer=standardization.scale))
}
groupCenter <- function(object, ...) UseMethod("groupCenter")
groupCenter.glm <- function(object, strata, samples=0,...) {
stopifnot(is.null(samples) || is.numeric(samples) || is.logical(samples) || is.data.frame(samples),
          !is.data.frame(samples) || nrow(samples)==length(object$linear.predictor)) 
  if (is.null(samples) || samples==0) 
  return(data.frame(PShat=groupCenter(object$linear.predictor, strata=strata)))
  
  if (is.numeric(samples) || is.logical(samples))
    {
       thePSes <- data.frame(PShat=object$linear.predictors,
                              samplePropensities(object, samples=as.numeric(samples)))
       nsamples <- as.numeric(samples)
  } else {
    thePSes <- data.frame(PShat=object$linear.predictors, samples)
    nsamples <- length(samples)
  }
       names(thePSes) <- c("PShat", paste("PSdraw", 1:nsamples, sep=""))
       return(groupCenter(thePSes, strata))
  }
groupCenter.clogit <- function(object,...)
  {
    linPSes <- as.vector(object$linear.predictors)
    names(linPSes) <- names(object$residuals)
    groupCenter(linPSes,...)
  }
groupCenter.data.frame <- function(object, strata, ...) {
mfd <- if (inherits(strata, "optmatch")) optmatch::matchfailed(strata) else logical(nrow(ans))
mtd.nms <- names(strata)[!is.na(strata) & !mfd]
theYs <- as.matrix(object[mtd.nms,])
ans <- resid(lm(theYs~strata[mtd.nms]))
as.data.frame(ans)
} #sapply(object, groupCenter, strata=strata)
groupCenter.numeric <- function(object, strata, ...) {
  stopifnot(inherits(strata, "factor"), length(object)==length(strata) || length(object)==sum(optmatch::matched(strata)), length(object)!=length(strata) ||!any(is.na(object) & !is.na(strata)), length(object)==length(strata) || !is.null(names(object)), length(object)==length(strata) || !any(is.na(object)))
mfd <- if (inherits(strata, "optmatch")) optmatch::matchfailed(strata) else logical(nrow(ans))
  mtd <- !is.na(strata) & !mfd
  if (length(object)!=length(strata)) mtd <- names(strata)[mtd]
resid(lm(object[mtd]~strata[mtd]))
}
PSl2sqtimesn <- function(object, omatch, alreadycentered=FALSE)
  {
    stopifnot(inherits(omatch, "optmatch"))
    mtd <- optmatch::matched(omatch)
    thepses <- if (is.data.frame(object) && alreadycentered)
      { object
      } else as.data.frame(groupCenter(object, strata=omatch))
    omatch <- omatch[mtd, drop=TRUE]
    
    thelens <- table(omatch)
    theweights <- unsplit(lapply(thelens, function(x) rep(1/x, x)),omatch)
     sapply(thepses, function(x) sum(theweights*x^2))
  }
PSoutlierquotient <- function(object, omatch, alreadycentered=FALSE)
  {
    stopifnot(inherits(omatch, "optmatch"))
    thepses <- if (is.data.frame(object) && alreadycentered)
      { object
      } else as.data.frame(groupCenter(object, strata=omatch))

    l2sqtimesn <- PSl2sqtimesn(thepses, omatch, alreadycentered=TRUE)
###    nn <- nlevels(omatch[optmatch::matched(omatch), drop=TRUE])

    sapply(thepses, function(x) max(abs(x)))/sqrt(l2sqtimesn)
  }
getMatchedPSDists <- function(object, ...) UseMethod("getMatchedPSDists")
getMatchedPSDists.data.frame <- function(object, omatch, standardization.scale=mad,...) 
  {
     stopifnot(class(omatch)[1]=="optmatch", length(omatch)==nrow(object), !any(is.na(object)))
     if (all(as.vector(object)>0 & as.vector(object)<1)) 
         warning("getMatchedPSDists() works on linear PSes, ie logits of est'd probabilities")
         theps <- object
         psnms <- names(object)
           thetx <- attr(omatch, "contrast.group")
     theses <- sapply(theps, makepooledsd, Tx=thetx, standardization.scale=standardization.scale)
     
    theps.tx <- split(theps[thetx,], omatch[thetx,drop=F], drop=F)
    theps.ctl <- split(theps[!thetx,], omatch[!thetx,drop=F], drop=F)
    legits <- sapply(theps.tx, nrow) & sapply(theps.ctl, nrow)
    thediffs <- levels(omatch)[legits]; names(thediffs) <- thediffs

    thediffs <- lapply(thediffs, function(nm) {sapply(psnms, function(nms) {abs(theps.tx[[nm]][[nms]] - theps.ctl[[nm]][[nms]])/theses[[nms]]})})         
###     mssizes <- sapply(thediffs, nrow)
###     weights <- rep(mssizes^-1, mssizes)/length(thediffs)
     ans <- as.data.frame(do.call(rbind, thediffs))
     ans
  }
getMatchedPSDists.glm <- function(object, omatch, standardization.scale=mad, ...)
  {
    stopifnot(class(omatch)[1]=="optmatch", length(omatch)==length(fitted(object)))
    theps <- object$linear.predictors
    thetx <- attr(omatch, "contrast.group")
    thesd <- makepooledsd(theps, Tx=thetx, standardization.scale=standardization.scale)
    theps.tx <- split(theps[thetx], omatch[thetx,drop=F], drop=F)
    theps.ctl <- split(theps[!thetx], omatch[!thetx,drop=F], drop=F)
    thediffs <- levels(omatch); names(thediffs) <- levels(omatch)
    thediffs <- lapply(thediffs, function(nm) {theps.tx[[nm]] - theps.ctl[[nm]]})
    abs(unlist(thediffs))/thesd

  }

