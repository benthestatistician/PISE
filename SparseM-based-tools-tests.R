require(SparseM)
##require("optmatch")
data(nuclearplants, package="optmatch")

test <- function(t, m = "Error!") {
  if (!t) {
    stop(m)
  } 
}
test.eq <- function(x,y, digits=12, m = "Error!") {
  if (any(abs(x-y)>= 10^(-digits))) { stop(m) } 
}


shouldError <- function(expr, msg = "Exception should be thrown") {
  r <- try(expr, silent = T)
  if (!inherits(r, "try-error")) {
    stop(msg)  
  }
}


themat <- as.matrix(subset(nuclearplants, select=t1:cap))
thefac <- as.factor(nuclearplants$pt)

### Using sparse matrices to strip out fixed effects
stratmeans1 <- fitted(lm(themat~thefac))
dimnames(stratmeans1) <- NULL
all.equal(stratmeans1,
          fitted(SparseM:::slm.fit(SparseMMFromFactor(thefac),
                                 themat)
                 )
          )

### NAs in the factor variable just don't get associated w/ any stratum
thefac1 <- thefac
thefac1[1] <- NA

all.equal(rep(0, nlevels(thefac)),
          as.matrix(SparseMMFromFactor(thefac1))[1,]
          )

factab <- as.vector(table(thefac))
stratum_summing_SparseM(thefac) %*%themat

ftds <- (function(mat, fac)
         {
           ssmat <- stratum_summing_SparseM(fac)
           dmat <- SparseMMFromFactor(fac)
           ssums <- ssmat%*%mat
           smns <- as.matrix(ssums)/factab
           ans <- as.matrix(dmat %*% smns)
           subset(ans, !is.na(fac))
         })(themat, thefac)

all.equal(ftds,
          fitted(SparseM:::slm.fit(SparseMMFromFactor(thefac),
                                 themat)
                 )
          )
