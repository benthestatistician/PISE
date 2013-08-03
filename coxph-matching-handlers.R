setOldClass("Surv")
setOldClass("coxph")
setOldClass(c("coxph.penal", "coxph"), S4Class="coxph")

convertSurv2to3 <- function(x)
  {
    stopifnot(inherits(x, "Surv"), packageVersion("survival") >= '2.37')
    if (ncol(x)==2)
      {
    colnames(x) <- gsub("time", "stop", colnames(x), fixed=TRUE)
    data.frame(start=rep(-Inf, nrow(x)), as.matrix(x))
  } else as.data.frame(as.matrix(x))
  }
# seems to require survival version 2.37 or better -- version 2.36-14 didn't have an as.matrix method, but 2.37-4 did.
setMethod("exactMatch", "Surv", function(x, row.names=NULL,...) #returnEventTimes=FALSE,
  {
stopifnot(#
          inherits(x, "Surv"),
          ncol(x) %in% 2L:3L,
          identical(colnames(x), c("time",   "status")) ||
          identical(colnames(x), c("start",  "stop",   "status")),
          is.null(row.names) || is.character(row.names),
          is.null(row.names) || length(row.names)==nrow(x)
          )

n <- nrow(x)
xT <- which(as.logical(x[,'status']))
xC <- setdiff(1L:n, xT)
x <- convertSurv2to3(x)

at_risk_when <- function(eventtime, x.controls)
  {
    which.of.controlrows <- with(x.controls,
                                 (start < eventtime) & (eventtime <= stop))
    which(which.of.controlrows)
  }

csForTs <- lapply(drop(x[xT,"stop"]), at_risk_when, x.controls=x[xC,])
cols <- unlist(csForTs)
tmp <- sapply(csForTs, length)
rows <- rep(1L:(length(csForTs)), tmp)
matrixcontents <- #if (returnEventTimes) rep(drop(x[xT,"stop"]), tmp)  else
  rep(0, length(rows))
nameConverter <- if (is.null(row.names)) as.character else function(x) row.names[x] 
tmp <- optmatch:::makeInfinitySparseMatrix(matrixcontents, cols = cols, rows =
                                rows, rownames = nameConverter(xT), colnames = nameConverter(xC))
return(tmp)
}
)


# NB: If standardization.scale has a different default here than with the other match_on methods,
# then if this is brought into optmatch the documentation of the other methods will have to be adjusted,
# in that the default value for standardization.scale will have to be noted in method-specific Details
# rather than in the param description.  Also, note that here the scale estimate is pooled across the
# interaction of any strata with the treatment variable, to bring it a bit closer to the spirit of the Cox
# model -- not that any of the rescaling stuff really makes sense here (IMHO).

#' method to assemble risk score distances in format suitable for matching
#'
#' Given a fitted \code{coxph} object modeling time-specific risks or hazards of falling into the treatment group, this function assembles into a distance suitable for matching units that at some point fell into the treatment group to units that never did so.  The risk differences extracted for matching are those based on the treatment group member's and the corresponding potential control group members' covariate profiles at the time that treatment group member fell into the treatment group.
#'
#' @title match_on method for coxph's
#' @param x a fitted model resulting from a call to \code{\link{survival::coxph}}
#' @param within specification of potential matches to be excluded a priori, typically produced by \code{\link{exactMatch}}, \code{\link{caliper}} or some combination thereof.
#' @param ... additional arguments for the \code{\link{match_on}}'s numeric method
#' @examples
#' a <- 1
#' class(a) <- "lm"
#' MyHappyFunction(a)
#'
#' @rdname match_on-methods
#' @docType methods
setMethod("match_on", "coxph", function(x, within = NULL, caliper = NULL, standardization.scale=NULL, ...)
          {
### First apply exactMatch() to the Surv within the coxph, creating 
### a within object encoding Tx vs Cntl, which observations at risk when each event happened, time of event.
### b/c coxph's don't store dimnames, assign names based on the row numbers of its data frame _before na handling happened_
            mf <- model.frame(x)
            tms.x <- terms(x)
            rnms <- row.names(mf)
            if (!is.null(attr(tms.x, "specials")$cluster))
              warning("Cox model involves clustering. Note that match_on() doesn't presently accommodate this in any way.")
            if (any(substr(attr(tms.x,"term.labels"),1,8)=="frailty("))
              warning("Cox model has a frailty term. Note that match_on() doesn't address this in any way.")
            stratfac <- if (is.null(attr(tms.x, "specials")$strata)) NULL else
            interaction(mf[attr(tms.x, "specials")$strata]) # if model has strata, risks are comparable only within them.
            within0 <- exactMatch(x$y, row.names=rnms) #, returnEventTimes=TRUE
            experiencedEvent <- as.logical(x$y[,"status"])
            names(experiencedEvent) <- rnms
            if (!is.null(stratfac))
              {
                within1 <- exactMatch(stratfac, experiencedEvent)
                within0 <- within0 + within1
              }            
            if (!is.null(within)) within0 <- within0 + within
            
            fitteds <- predict(x, type="lp")
            names(fitteds) <- rnms

            if (!is.null(standardization.scale))
              {
                stratfac <- if (is.null(stratfac)) experiencedEvent else interaction(stratfac,experiencedEvent)
                stratfac <- factor(stratfac) # makes sure unused levels get ditched
                spreads <- tapply(fitteds, stratfac,
                                  function(x, scaler) {(length(x)-1) * scaler(x)^2},
                                  scaler=standardization.scale,simplify=TRUE)
                scale <- sqrt(sum(spreads)/(length(stratfac)-nlevels(stratfac)))
                fitteds <- fitteds/scale
              }
            match_on(x=fitteds, within=within0, caliper=caliper, z=experiencedEvent, ...)
          }
)


##' Method to create smaller \code{InfinitySparseMatrix}-es from larger ones, by 
##' collapsing blocks of a given \code{InfinitySparseMatrix} using a user-supplied function.
##' 
##' The \code{InfinitySparseMatrix}, \code{x}, must have non-NULL row and column names. \code{by} is required to have matching names: If a factor is supplied as the \code{by} argument, it must have names; a data frame is required to have row names.  In either case, the names must include each character appearing in the row and column names of x.
##'
##' For each level of \code{by}, unit names corresponding to that level either must all appear within the row names of \code{x} or must all appear within the column names of \code{x}.  This behaviour may change once I figure out appropriate conventions for handling situations when this is not true.
##'
##' (In the future it might be nice to add an optional \code{data} argument, permitting \code{by} to be specified as a formula or expression that then gets evaluated within \code{data}.)
##' @title 
##' @param x An InfinitySparseMatrix object, as created by \code{\link{match_on}}, \code{\link{caliper}} or \code{\link{exactMatch}}.
##' @param by Either a factor or a data frame containing a factor (and nothing else); see details.
##' @param FUN A function with which to collapse the blocks of x.  Must return a scalar.
##' @param ... Additional arguments to be passed to \code{FUN}.
##' @return An \code{InfinitySparseMatrix}, with row and column names taken from levels of \code{by}
##' @author Ben Hansen
setMethod("aggregate", "InfinitySparseMatrix",  function(x, by, FUN, ...)
  {
    stopifnot(!is.null(x@colnames), !is.null(x@rownames),
              is.factor(by) || is.character(by) || is.data.frame(by),
              !is.data.frame(by) || (length(by)==1 && (is.factor(by[[1]]) || is.character(by[[1]]))),
              !(is.factor(by) || is.character(by)) || !is.null(names(by)),
              !is.data.frame(by) || !is.null(row.names(by))
              )
    stop("aggregate method for InfinitySparseMatrix-es not implemented yet")
  }
    )

boxplot.coxph <- function(x, data=NULL, xlab="Group", ylab=expression(paste(X, symbol("\242"), hat(beta))), main="Overlap on fitted scores",varwidth=TRUE,...)
  {
    if (!is.null(data)) warning("data argument currently ignored")
theSurv <- convertSurv2to3(x$y)
groups <- theSurv[[3]]
linear.score <- x$linear.predictors
boxplot(linear.score ~ groups, xlab=xlab, ylab=ylab,main=main,...)
  }
