##' <description>
##' Turn a factor variable into a sparse matrix of 0's and 1's, such that if observation i 
##' has the jth level then there is a 1 at position (i,j) (but nowhere else in row i).
##' <details>
##' NA's give rise to rows with no 1s.  
##' As the result is only meaningful in the context of the SparseM package,
##' function requires that SparseM be loaded.
##' @title Represent factor levels as columns of a sparse matrix
##' @param thefactor Factor variable, or object inheriting from class factor
##' @return Sparse csr matrix the columns of which are dummy variables for levels of thefactor 
##' @author Ben Hansen
SparseMMFromFactor <- function(thefactor)
  {
stopifnot(require("SparseM"), inherits(thefactor, "factor")
          )
theNA <- ##if (inherits(thefactor, "optmatch")) !matched(thefactor) else
is.na(thefactor)

  if (all(theNA)) stop("No non-NA's in thefactor") else {
  if (any(theNA) && !inherits(thefactor, "optmatch")) warning("NA's found in thefactor.")
}

nlev <- nlevels(thefactor)
nobs <- length(thefactor)
theint <- as.integer(thefactor)
if (any(theNA)) theint[theNA] <- 1L#nlev + 1L:sum(theNA)
new("matrix.csr",
    ja=theint,
    ia=1L:(nobs+1L),
    ra=(1L-theNA),
    dimension = c(nobs, nlev) #+sum(theNA)
    )
  }

##' .. content for \description{} (no empty lines) ..
##' Given a factor variable, generate a sparse matrix 
##' suitable for left-multiplying data matrices in order
##' to assemble within-category sums.
##' .. content for \details{} ..
##' @title Sparse matrix for stratum sums
##' @param thefactor 
##' @return Sparse matrix of class csr (SparseM package)
##' @author Ben B Hansen
stratum_summing_SparseM <- function(thefactor)
  {
    stopifnot(require("SparseM"), inherits(thefactor, "factor")
              )
    if (all(is.na(thefactor))) warning("No non-NA's in thefactor") 

    nn <- length(thefactor)
    ind.splitted <- split(1L:nn, thefactor)
    stratsizes <- sapply(ind.splitted, length)
    rowstarts <- 1L+c(0L,as.integer(cumsum(stratsizes)))
    nobs <- sum(stratsizes)
    
    new("matrix.csr",
        ja=unlist(ind.splitted),
        ia=rowstarts,
        ra=rep(1L,nn),
        dimension=c(nlevels(thefactor), nn)
        )
  }

