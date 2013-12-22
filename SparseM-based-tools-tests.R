source("SparseM-based-tools.R")
library(testthat)
library(SparseM)
data(nuclearplants, package="optmatch")

context("SparseM fixed effect design matrices, stratum sums.")

test_that("Return value of SparseMMFromFactor is the f.e. design matrix\n corresponding to provided factor", {
themat <- as.matrix(subset(nuclearplants, select=t1:cap))
thefac <- as.factor(nuclearplants$pt)

### Using sparse matrices to strip out fixed effects
stratmeans1 <- fitted(lm(themat~thefac))
dimnames(stratmeans1) <- NULL
expect_equal(stratmeans1,
             fitted(SparseM:::slm.fit(SparseMMFromFactor(thefac),
                                      themat)
                    )
             )
})

test_that("NAs in factor variable provided to SparseMMFromFactor simply ignored\n (rows not dropped, but are all 0s)",
          {
              
              themat <- as.matrix(subset(nuclearplants, select=t1:cap))
              thefac1 <- as.factor(nuclearplants$pt)
              thefac1[1] <- NA
              expect_warning(thefac1_SM <- SparseMMFromFactor(thefac1),
                             "NA")

              expect_equal(length(thefac1), nrow(as.matrix(thefac1_SM)))
              expect_equal(rep(0, nlevels(thefac)),
                           as.matrix(thefac1_SM)[1,]
                           )
          })

test_that("SparseM based tabulation of stratum sums (whence stratum means also)",
          {
              
              themat <- as.matrix(subset(nuclearplants, select=t1:cap))
              thefac <- as.factor(nuclearplants$pt)
              stratsums <- sapply(split(as.data.frame(themat), thefac), colSums)
              stratsums <- t(stratsums)
              dimnames(stratsums) <- NULL
              expect_equal( stratsums,
             as.matrix(stratum_summing_SparseM(thefac) %*%
                       themat) # dim=strata x vars
              )
          })

test_that("stratMeans calculates by-stratum means of one or more numeric vars",
{
    expect_equal(matrix(tapply(themat[,1,drop=T], thefac, mean), nlevels(thefac),1),
                 stratMeans(themat[,1,drop=T], thefac))

    smns <- stratMeans(themat, thefac)
              
              expect_equal(as.matrix(SparseMMFromFactor(thefac) %*% smns),
                           fitted(SparseM:::slm.fit(SparseMMFromFactor(thefac),
                                                    themat)
                                  )
                           )
          })
                           
### Should add tests of NA behavior too.
