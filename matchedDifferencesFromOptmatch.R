##' Matched differences in an arbitrary real variable from the match
##'
##' This function builds a matrix to left-multiply a vector or design matrix by,
##' yielding a matrix of matched differences.  It also passes on some information
##' tabulated in the process of making this matrix that is useful for aligning and
##' checking the results, and for averaging matched differences together in ways that
##' respect differences in sizes of the matched sets.
##'
##' @param amatch Factor variable encoding a match, as created by optmatch functions
##' @param type Class of matrix to be created, defaulting to csr class of package SparseM
##' @return List with components: \code{diffMaker}, the matrix to left-multiply vectors or
##' matrices by, either in ordinary dense matrix format or in the 'csr' format, as implemented
##' by SparseM package; \code{isOneMany}, logical vector of length \# of matched sets
##' indicating whether matched set is one-many, ie on treatment and one or more controls, or a
##' many-one; \code{theOne}, character vector of length \# of matched sets recording name of
##' "The One," ie the treatment if it's a one-many matched set and otherwise the single control;
##' \code{theOthers.HowMany}, integer vector indicating how many there are in the matched set
##' other than the one, i.e. size of matched set minus one; \code{rownames}, \code{colnames},
##' character vectors giving row and column names of the matrix.
##' @author Ben Hansen \email{ben.hansen@@umich.edu}
##' @examples
##' require(optmatch); data(nuclearplants)
##' afm <- fullmatch(caliper(.2, glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants) ) )
##' stratumStructure(afm)
##' (afm.md.d <- MatchedDiffMaker(afm, type='dense'))
##'
MatchedDiffMaker <- function(amatch, type=c("SparseM-csr","dense"))
  {
stopifnot( inherits(amatch, "optmatch"), !is.na(pmatch(type[1], c("dense", "SparseM-csr", "sparsem-csr"))) )

TYPES <- c("SparseM-csr","dense")
type <- pmatch(type, c(TYPES, tolower(TYPES)))[1]
type <- rep(TYPES,2)[type]

N <- length(amatch)
theColNames <- names(amatch)
thematch <- amatch[matched(amatch), drop=TRUE]
K <- nlevels(thematch)
tx <- attr(thematch, "contrast.group")
isManyOne <- tapply(tx, thematch, function(x) sum(x)>1)

names(tx) <- names(thematch)
theOne <- tapply(tx, thematch, function(x)
                 { if (sum(x)==1) names(x)[which.max(x)] else names(x)[which.min(x)] })
theRest <- tapply(tx, thematch, function(x)
                 { if (sum(x)==1) names(x)[-which.max(x)] else names(x)[-which.min(x)] },
                  simplify=FALSE)
nc <- sapply(theRest, length) # the no. of controls or, for many-one matches, the no. of tx
NROW <- sum(nc)
theRowNames <- unlist(theRest)

if (type[1]=="dense")
  {
mat <- matrix(integer(NROW*N), nrow=NROW, ncol=N) #, dimnames=list(theRowNames, theColNames)
}
if (type[1]=="SparseM-csr")
  {
    if (!require("SparseM")) stop("Can't create sparse matrix because can't load package SparseM.")

mat <- as.matrix.csr(rep(0,NROW*N), nrow=NROW, ncol=N)
  }

for (i in names(theOne))
  {
    therows <- match(theRest[[i]], theRowNames)
    theOneCol <- match(theOne[i], theColNames)
    theRestCols <- match(theRest[[i]], theColNames)
### NB:In SparseM case should be possible to speed up the below substantially by filling the matrices directly.
### The csr format is by row, and this loop is essentially filling by row.
    if (!isManyOne[i])
      {
        mat[therows,theOneCol]<- 1L
        mat[therows,theRestCols] <- diag(-1L, nrow=nc[i] )
      } else
    {
      mat[therows,theOneCol] <- -1L
      mat[therows,theRestCols] <- diag(1L, nrow=nc[i] )
    }
  }

list(diffMaker=mat, isOneMany=!isManyOne, theOne=theOne, theOthers.HowMany=nc, rownames=theRowNames, colnames=theColNames)
}

##' Helper function to \code{\link{MatchedDiffMaker}} to create readable tables of matched differences
##'
##' @param x A design matrix, aligned by row with the optmatch object used to create \code{matchedDiffList}
##' @param matchedDiffList The value of a call to \code{\link{MatchedDiffMaker}}
##' @return Data frame displaying name of treatment, name of control, size of matched set and matched differences
##' @author Ben Hansen
##' @examples
##' require(optmatch); data(nuclearplants)
##' afm <- fullmatch(caliper(.2, glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants) ) )
##' stratumStructure(afm)
##' afm.md.d <- MatchedDiffMaker(afm, type='SparseM-csr')
##' MatchedDiffTable(nuclearplants[c("t1", "t2", "pt")], afm.md.d)
##' MatchedDiffTable(nuclearplants$cost, afm.md.d)
##'
MatchedDiffTable <- function(x, aMDM)
  {
stopifnot(is.numeric(as.matrix(x)),
          is.null(dim(x)) || isTRUE(all.equal(aMDM$colnames, row.names(x))))

xnm <-  if (is.null(dim(x)))
  {
    paste(deparse(substitute(x), width.cutoff = 500L), collapse = " ")
  } else if (is.null(colnames(x))) paste("V",1L:ncol(x)) else colnames(x)

ans <- with(aMDM,
            data.frame(theTx=ifelse(rep(isOneMany, theOthers.HowMany), rep(theOne, theOthers.HowMany), rownames),
                       theCtl=ifelse(rep(isOneMany, theOthers.HowMany), rownames, rep(theOne, theOthers.HowMany)),
                       msSize=1+rep(theOthers.HowMany, theOthers.HowMany)
                       )
            )
diffs <- aMDM$diffMaker%*%as.matrix(x)
diffs <- as.matrix(diffs)

dim(diffs) <- c(nrow(ans), length(xnm))
colnames(diffs) <- xnm
data.frame(ans, diffs)
  }

##' <description>
##' Function to calculate matched differences on linear propensity score, with associated SEs.
##' <details>
##' To compute standard errors be computed by Huber-White and with conditioning on treatment margins of matched sets, set \code{covariance.extractor} to \code{SandwichVcov}, as defined in file sandwich_glm.R.  \code{covariance.extractor} could also be set to \code{vcov} or \code{sandwich}, from package \sQuote{sandwich}.
##' @title Matched differences on given fitted propensity score, with SEs.
##' @param PSmodel A fitted propensity score model, e.g. object of class glm
##' @param amatch  An optmatch object, ie the value of a call to fullmatch() or pairmatch()
##' @param covariance.extractor Function used to extract covariance estimate from \code{PSmodel}.
##' @return data frame with components: \code{theTx} and \code{theCtl},
##' character vectors giving names of the treatment and control units
##' being compared; \code{diff}, treatment-minus-control difference on
##' the fitted linear propensity score; \code{SE}, corresponding standard error.
##' @author Ben Hansen  
PSdiffs <- function(PSmodel, amatch, covariance.extractor=vcov)
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
        stop("Something is wrong in PSdiffs().  Please debug!")
    coefnames <- names(thecoef)
    thecov <- covariance.extractor(PSmodel, strata=amatch)
    thecov <- thecov[coefnames, coefnames]

    aMDM <- if (require("SparseM"))
      {MatchedDiffMaker(amatch, type="SparseM-csr")
     } else MatchedDiffMaker(amatch, type="dense")
    xdiffs <- as.matrix(aMDM$diffMaker %*% mm)
    PSdiffs <- xdiffs %*%thecoef
    MakePSSE <- function(thecov) {
    PSSEs <- rowSums( (xdiffs%*%thecov)*xdiffs )
    sqrt(PSSEs)
  }
#    ps.ses <- as.data.frame(lapply(thecovlist, MakePSSE))
    ps.se <- MakePSSE(thecov)
    theOne <- with(aMDM, rep(theOne,theOthers.HowMany))
    theOther <- aMDM$"rownames"
    isOneMany <- with(aMDM, rep(isOneMany,theOthers.HowMany))
    ans <- data.frame(theTx=as.character(ifelse(isOneMany, theOne, theOther)),
                      theCtl=as.character(ifelse(isOneMany, theOther, theOne)),
                      diff=PSdiffs, #*ifelse(isOneMany,-1,1),
                      SE=ps.se,#ps.ses,
                      stringsAsFactors=FALSE)
    ans
  }


model.matrix.bayesglm <- function(x,...)
  {
model.matrixBayes(terms(x),data=model.frame(x),
                        contrasts.arg=x$contrasts,
        		keep.order=x$keep.order,
                        drop.baseline=x$drop.baseline,...)
  }


  
