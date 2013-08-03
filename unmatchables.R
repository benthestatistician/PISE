###Some functions for inferring directly from a distance who can be matched and who can't.
###I think this should be rewritten a little, to support two related functions sharing a help page.
###the workhorse function would be unmatchables(), which would return a list with elements `0` and `1`,
### character vectors giving names of tx and control units that are unmatchable.  This would be suitable
### for passing as exclude= argument to caliper().  Layered over that would be a function matchable(), which
### taking an (optional, I think) second data= argument and returning a logical vector identifying matchable units.  The
### main usage case for this one is when re-fitting a propensity model to a data frame excluding the unmatchables.
### Re usage w/o an explicit data frame being passed: the main case I have in mind for this is in refitting a PS
### giving argument "subset=matchable(myISM)".  In this case it would ideally use the data= arg to the propensity
### scoring function for row names, throwing an error if no such arg were given.  Perhaps there's some way to
### to leverage the lazy evaluation model to achieve this behavior.
###need tests and definitions!
setGeneric("matchable", def= function(x, data=NULL) standardGeneric("matchable"))
unmatchable <- function(...) !matchable(...)
setMethod("matchable", "InfinitySparseMatrix", function(x,data=NULL) {
if (is.null(x@rownames) || is.null(x@colnames)) stop("Match distances must have row and column names.")

if (is.null(data))
  {
    result <- logical(sum(dim(x)))
    names(result) <- unlist(dimnames(x))
  } else
{
  data <- as.data.frame(data)
  result <- logical(nrow(data))
  names(result) <- row.names(data)
  if (!all(x@rownames %in% row.names(data))) stop("Match distance row names don't appear in rownames(data)")
  if (!all(x@colnames %in% row.names(data))) stop("Match distance column names don't appear in colnames(data)")
}

matchable.txes <- x@rows[is.finite(x@.Data)]
matchable.ctls <- x@cols[is.finite(x@.Data)]
result[x@rownames[matchable.txes]] <- TRUE
result[x@colnames[matchable.ctls]] <- TRUE
result
})
setMethod("matchable", "DenseMatrix", function(x,data=NULL) {
if (is.null(rownames(x)) || is.null(colnames(x))) stop("Match distances must have row and column names.")

if (is.null(data))
  {
    result <- logical(sum(dim(x)))
    names(result) <- unlist(dimnames(x))
  } else
{
  data <- as.data.frame(data)
  result <- logical(nrow(data))
  names(result) <- row.names(data)
  if (!all(rownames(x) %in% row.names(data))) stop("Match distance row names don't appear in rownames(data)")
  if (!all(colnames(x) %in% row.names(data))) stop("Match distance column names don't appear in colnames(data)")
}
matchable.txes <- apply(x, 1, function(x) any(is.finite(x)))
matchable.ctls <- apply(x, 2, function(x) any(is.finite(x)))
result[rownames(x)[matchable.txes]] <- TRUE
result[colnames(x)[matchable.ctls]] <- TRUE
result
})
setMethod("matchable", "optmatch.dlist", function(x,data=NULL) {
if (is.null(data))
  {
    result <- logical(length(attr(x, "row.names")))
    names(result) <- attr(x, "row.names")
  } else
{
  data <- as.data.frame(data)
  result <- logical(nrow(data))
  names(result) <- row.names(data)
  if (!all(attr(x, "rownames") %in% row.names(data))) stop("Matching unit names don't all appear in rownames(data)")
}
matchables <- tapply(x, function(dmat) {
matchable.txes <- apply(x, 1, function(x) any(is.finite(x)))
matchable.ctls <- apply(x, 2, function(x) any(is.finite(x)))
c(rownames(x)[matchable.txes], colnames(x)[matchable.ctls])
}  )
matchables <- unlist(matchables)
result[matchables] <- TRUE
result
})
