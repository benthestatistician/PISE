# C-c C-o
##' .. content for \description{} (no empty lines) ..
##' Function to extract estimates of center from randomization distributions or propensity-tilted randomization distributions.
##' @title Dispersion estimates from (tilted) matched randomization ECDFs
##' @param rdObj Object of class \code{rdObj} (containing simulated matched random assignments), as produced by \code{rdist}.
##' @param FUN  Function to be applied to \code{rdObj$statistic[[i]]} before taking (weighted or unweighted) means.
##' @param propensity A fitted propensity score model, for weighting reflecting tilting of assignments, or \code{NULL}, for no weighting.
##' @return List with components \sQuote{muhat}, the (weighted) mean of values of \code{FUN}, and \sQuote{se}, the corresponding standard error(s).
##' @author Ben Hansen
estimateMean <- function(rdObj, FUN=function(x) x, propensity=NULL)
  {
stopifnot(inherits(rdObj, "rdObj"),
          all.equal(names(rdObj$statistic), names(rdObj$dP.over.dShuffle.selection.probabilities))
          )
nrep <- length(rdObj$statistic)
dP.over.dS <- rdObj$dP.over.dShuffle.selection.probabilities
dQ.over.dS <- if (is.null(propensity)) rep(1, nrep) else dP.over.dS * tilt.from.permdist(rdObj, propensity=propensity)

statsM1 <- FUN(rdObj$statistic[[1]])*dQ.over.dS[1]
statsM2 <- statsM1^2

if (nrep>1)
  {
    for (ii in 2:nrep)
      {
        incr <- FUN(rdObj$statistic[[ii]])*dQ.over.dS[ii]
        statsM1 <- statsM1 + incr
        statsM2 <- statsM2 + incr^2
        NULL
      }
  }

statsM1 <- statsM1/nrep
statsM2 <- statsM2/nrep
SE <- sqrt((statsM2 - statsM1^2)/nrep)
return(list(muhat=statsM1, se=SE))
}

##' .. content for \description{} (no empty lines) ..
##' Function to extract dispersion estimates from randomization distributions or propensity-tilted randomization distributions.
##' .. content for \details{} ..
##' The standard error formula applies delta method along with result that V(s\^2) =  sigma\^4(2/(n-1) + kurtosis/n).  As of spring 2011 the function was giving funny answers -- not sure why.
##' @title Dispersion estimates from (tilted) matched randomization ECDFs
##' @param rdObj Object of class \code{rdObj} (containing simulated matched random assignments), as produced by \code{rdist}.
##' @param FUN  Function to be applied to \code{rdObj$statistic[[i]]} before taking (weighted or unweighted) means.
##' @param propensity A fitted propensity score model, for weighting reflecting tilting of assignments, or \code{NULL}, for no weighting.
##' @return List with components \sQuote{sigmahat}, estimated sigma for values of \code{FUN}, and \sQuote{se}, the corresponding standard error(s).
##' @author Ben Hansen
estimateSD <- function(rdObj, FUN=function(x) x, propensity=NULL)
  {
stopifnot(inherits(rdObj, "rdObj"),
          all.equal(names(rdObj$statistic), names(rdObj$dP.over.dShuffle.selection.probabilities)),
          (nrep <- length(rdObj$statistic))>1 )

dP.over.dS <- rdObj$dP.over.dShuffle.selection.probabilities
dQ.over.dS <- if (is.null(propensity)) rep(1, nrep) else dP.over.dS * tilt.from.permdist(rdObj, propensity=propensity)

statsM1 <- FUN(rdObj$statistic[[1]])*dQ.over.dS[1]
for (ii in 2:nrep)
  {
    statsM1 <- statsM1 + FUN(rdObj$statistic[[ii]])*dQ.over.dS[ii]
  }
muhats <- statsM1/nrep

M2 <- (FUN(rdObj$statistic[[1]]) - muhats)^2*dQ.over.dS[1]
M4 <- (FUN(rdObj$statistic[[1]]) - muhats)^4*dQ.over.dS[1]
for (ii in 2:nrep)
  {
    incr <- FUN(rdObj$statistic[[ii]]) - muhats
    M2 <- M2 + incr^2*dQ.over.dS[ii]
    M4 <- M4 + incr^4*dQ.over.dS[ii]
  }
sigmahat <- sqrt(M2/(nrep-1))
M2 <- M2/nrep
M4 <- M4/nrep
kurt <- M4/(M2^2) - 3
return(list(sigmahat=sigmahat,
            se=.5*sigmahat*sqrt(2/(nrep-1) + kurt/nrep))
### formula applies delta method along with result that V(s^2) =  \sigma^4(2/(n-1) + \kappa/n)
       )
}

##' .. content for \description{} (no empty lines) ..
##' Turn matched set shuffles into EDCFs corresponding to matched randomization.
##' .. content for \details{} ..
##' There should ordinarily be an argument \code{statistic.}, a function accepting a treatment vector as its first argument.  This function will also get each of the arguments in ... that are passed to \code{rdist}.
##' @title Create EDCFs of matched randomization distributions.
##' @param object an \code{msShuffles} object or the ingredients of one, i.e. an \code{assignmentRecipe} object or an \code{optmatch} object. 
##' @param ... Additional arguments to dispatched methods, ordinarily including argument \code{statistic.} of mode \code{function}.  See details.
##' @return Object of class \code{rdObj}, to which \code{estimateMean}, \code{estimateSD},... can be applied.
##' @author Ben Hansen
rdist <- function(object,...)  UseMethod("rdist")
rdist.msShuffles <- function(object,statistic.,...)
{
nreps <- length(object)
Tx.templ <- attr(object, "contrast.group.template")
manyone <- attr(object, "manyone")
stopifnot(all.equal(row.names(object), names(manyone)))

ans <- list(shuffles=object,
            dP.over.dShuffle.selection.probabilities=1/tilt.from.permdist(object)
            )

cstat <- function(compactTx)
    statistic.(reconstituteTx(compactTx, manyone=manyone, unpackedTx=Tx.templ),...)
###thestats <-
ans$statistic <- lapply(object, cstat)

### ans$statistic <-
###    if (is.array(thestats[[1]]))
###    {
###        ds <- dim(thestats[[1]])
###        dn <- dimnames(thestats[[1]])
###        array(unlist(thestats),
###              dim=c(ds, nreps),
###              dimnames=c(dn, list(rep=1L:nreps))
###              )
###    } else if (length(thestats[[1]])>1)
###    {
###        dn <- names(thestats[[1]])
###        matrix(unlist(thestats), length(thestats[[1]]), length(thestats),
###               dimnames=c(list(statistics=dn), list(reps=1:nreps)))
###    } else {unlist(thestats) }

class(ans) <- c("rdObj", class(ans))
ans
}
rdist.default <- function(object,...)
{
newobject <- msShuffler(object,...)
rdist(newobject,...)
}

msShuffler <- function(object,...) UseMethod("msShuffler")
msShuffler.optmatch <- function(object,propensity=NULL,...)
{
    newobject <- makeAssignmentRecipe(propensity,object)
    msShuffler(newobject,...)
}
msShuffler.assignmentRecipe <- function(object,shuffles=2,...)
{
stopifnot(all.equal(names(object), names(manyone <- attr(object, "manyone"))))

samples <- sapply(object, function(lpses)
       sample(x=as.integer(names(lpses)), # arrange so that this is integer
              size=shuffles,replace=TRUE,prob=exp(lpses)))

ans <- as.data.frame(t(samples))
names(ans) <- 1L:length(ans)

attributes(ans) <-
    c(attributes(ans),
      attributes(object)[c("names.all.units", "contrast.group.template",
                           "manyone", "mssize", "contrast.group",
                           "positions.of.matched.units",
                           "matched.set.mean.linear.propensities",
                           "linear.propensities")])

attributes(object) <- attributes(object)['names']
attr(ans, "matched.set.linear.propensities") <- object

class(ans) <- c("msShuffles", class(ans))
ans
}
msShuffler.default <- function(object,...) stop("msShuffler requires an assignmentRecipe or optmatch argument.")

tilt.from.permdist <- function(x,...) UseMethod("tilt.from.permdist")
tilt.from.permdist.msShuffles <- function(x, propensity=NULL,...)
  {
    stopifnot(is.null(propensity)|| is.numeric(propensity) || inherits(propensity, "glm"))
    if (is.numeric(propensity) && length(propensity)==1) return(rep(1, length(x)))
    lpses <- if (is.null(propensity))
      {
        attr(x, "matched.set.linear.propensities")
      } else {
        lps <- if (inherits(propensity, "glm")) propensity$linear.predictors else propensity
        lapply(attr(x,"matched.set.linear.propensities"),
               function(x)
               {
                 aa <- lps[as.integer(names(x))]
                 names(aa) <- names(x)
                 aa
               } )
      }
    diff.of.cumfcts <- sum( sapply(lpses, function(lps) log( sum(exp(lps))/length(lps) ) ) )
    lpses.mses.by.shuffles <-
      sapply(row.names(x), function(msnm) lpses[[msnm]][as.character(x[msnm,])])
    rownames(lpses.mses.by.shuffles) <- colnames(x)
    shuffle.sum.lpses <- apply(lpses.mses.by.shuffles, 1, sum)
    exp(shuffle.sum.lpses - diff.of.cumfcts)
  }
tilt.from.permdist.rdObj <- function(x,...) tilt.from.permdist(x$shuffles,...)
tilt.from.permdist.default <- function(x,...) stop("tilt.from.permdist() requires input of class 'msShuffles'.")

reconstituteTx <- function(compactTx, manyone, unpackedTx)
  {
    unpackedTx[compactTx[!manyone]] <- TRUE
    unpackedTx[compactTx[manyone]] <- FALSE
    as.numeric(unpackedTx)
  }
makeAssignmentRecipe <- function(object,strata,...) UseMethod("makeAssignmentRecipe")
makeAssignmentRecipe.glm <- function(object,strata,...)
{
stopifnot(object$family$family=="binomial")
makeAssignmentRecipe(object$linear.predictors,strata,...)
}
makeAssignmentRecipe.numeric <- function(object, strata, truncate.linear.ps.at=Inf,...)
{

     stopifnot(inherits(strata, "optmatch"),
               is.numeric(truncate.linear.ps.at),
               truncate.linear.ps.at>0,
               all.equal(names(strata), names(object)))
     strata <- presplitOptmatch(strata)
     msmeanPSes <- numeric(length(strata))
     names(msmeanPSes) <- names(strata)
     manyone <- attr(strata, "manyone")

     for (lev in names(strata))
         {
             mo <- manyone[lev]
             snms <- names(strata[[lev]])
             lps <- object[as.integer(snms)]
             msmeanPSes[lev] <- msmnPS <- mean(lps)
             strata[[lev]] <- (1- 2*mo)*(lps - msmnPS)
             if (is.finite(truncate.linear.ps.at))
               strata[[lev]] <-
                 pmin(truncate.linear.ps.at, pmax(-1*truncate.linear.ps.at, strata[[lev]]))

             names(strata[[lev]]) <- snms
             NULL
         }

     attr(strata, "matched.set.mean.linear.propensities") <-
       if (is.finite(truncate.linear.ps.at)) NULL else msmeanPSes

     attr(strata, "linear.propensities") <- object
     class(strata) <- c("assignmentRecipe", class(strata))
     strata
}
makeAssignmentRecipe.default <- function(object, strata,...)
{ # To be invoked when propensity slot is NULL

if (!is.null(object)) stop("I don't know how to handle 1st args (propensities) of this type.")

     strata <- presplitOptmatch(strata)
     msmeanPSes <- numeric(length(strata))
     names(msmeanPSes) <- names(strata)
     attr(strata, "linear.propensities") <- 0
     attr(strata, "matched.set.mean.linear.propensities") <- msmeanPSes
     class(strata) <- c("assignmentRecipe", class(strata))
     strata

}

presplitOptmatch <- function(omatch)
{
stopifnot(inherits(omatch, "optmatch"))

tx <- attr(omatch, "contrast.group")
mtd.int <- which(matched(omatch))
tx.sm <- tx[mtd.int]
om.sm <- omatch[mtd.int,drop=T]
lpses <- numeric(length(mtd.int))
names(lpses) <- mtd.int #names(om.sm)
ans <- split(lpses, om.sm)

names(mtd.int) <- names(om.sm)
attr(ans, "positions.of.matched.units") <- mtd.int
attr(ans, "contrast.group") <- attr(omatch, "contrast.group")
attr(ans, "mssize") <- sapply(ans, length)
attr(ans, "manyone") <- manyone <-
    tapply(tx.sm,om.sm,function(tvec) as.logical(sum(tvec)-1))
attr(ans, "contrast.group.template") <- tx
for (lev in names(ans))
    attr(ans, "contrast.group.template")[as.integer(names(ans[[lev]]))] <- rep(manyone[lev],length(ans[lev]))
attr(ans, "names.all.units") <- names(omatch)
ans
}
