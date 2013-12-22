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
##' .. content for \description{} (no empty lines) ..
##' Stratum-wise means for a numeric variable or matrix of numeric variables.
##' Based on sparse matrix multiplication, to limit memory consumption.
##' .. content for \details{} ..
##' The within-stratum calculation is effectively \code{mean(x,na.rm=T)}.
##' Stratum-variable combinations for which there are no non-NA 
##' observations receive the value NA; this is the only circumstance
##' under which NAs are returned.
##' @title Stratum wise means
##' @param mat numeric vector or matrix (\code{matrix.csr}s are allowed)
##' @param fac factor defining the strata
##' @return A matrix (dense) of dimension \code{nlevels(fac)} by \code{ncol(mat)}
##' @author Ben B Hansen
stratMeans <- function(mat, fac)
                       {
                           stopifnot(is.numeric(mat),
                                     inherits(mat, "matrix") || inherits(mat, "matrix.csr") || is.null(dim(mat)),
                                     is.null(dim(mat)) || nrow(mat)==length(fac),
                                     !is.null(dim(mat)) || length(mat)==length(fac)
                                     )
                                     
                           # calculate sums w/in each stratum/matched set
                           ssmat <- stratum_summing_SparseM(fac)
                           mat_notisNA <- !is.na(mat)
                           mat_NAas0 <- ifelse(mat_notisNA, mat, 0)
                           ssums <- ssmat%*%mat_NAas0
                           ssums <- as.matrix(ssums)

                           # stratum counts of non-missing obses
                           dm <- dim(mat_notisNA)
                           mat_notisNA <- as.numeric(mat_notisNA)
                           dim(mat_notisNA) <- dm
                           snotNA <- ssmat%*%mat_notisNA
                           snotNA <- as.matrix(snotNA)
                           snotNA[snotNA==0] <- NA
                           ssums/snotNA
                       }
