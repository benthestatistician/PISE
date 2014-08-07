## test coercion of optmatch to matrix.csr:
## - gives proper mean diffs w/ matched pairs
## - gives proper mean diffs w/ matches w/ mult controls
## - gives proper mean diffs w/ matches w/ varying num controls
## - gives mean diffs of 0 for matched sets w/ no tx or no ctl
## - associates levels(from)[i] with row i of result

library(optmatch)
library(testthat)

source("mlm.R")

test_that("parseMatchingProblem", {

  data(nuclearplants)
  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)

  n2 <- cbind(nuclearplants, mmm)
  res <- parseMatchingProblem(cost ~ pr + pt + mmm, n2)
  
  expect_is(res$mf, "data.frame")
  expect_is(res$match, "integer")
  expect_equal(length(res$match), 1)

  expect_equivalent(res$mf[, res$match], mmm)
  expect_equal(dim(res$mf)[2], 4) # cost, pr, pt, and mmm

  mmm2 <- fullmatch(pr ~ t1, data = nuclearplants)
  n3 <- cbind(n2, mmm2)
  res2 <- parseMatchingProblem(cost ~ pr + mmm, n3)
  expect_equivalent(res2$mf[, res2$match], mmm)
  expect_error(parseMatchingProblem(cost ~ pr + mmm + mmm2, n3), "one")
               
  expect_error(parseMatchingProblem(cost ~ pr, n3), "include")
})

test_that("optmatch -> matrix.csr", {
  data(nuclearplants)

  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  csr <- as(mmm, "matrix.csr")

  expect_is(csr, "matrix.csr")
  expect_equal(dim(csr), c(nlevels(mmm), length(mmm)))

  csrm <- as.matrix(csr) # cast it back to dense matrix for testing
  grps <- apply(csrm, 2, function(col) { which(col != 0) })

  tmp <- table(grps, mmm)
  expect_true(all(diag(tmp) != 0))
  diag(tmp) <- 0
  expect_true(all(tmp == 0))

  expect_true(all(rowSums(csrm) == 0))

  expect_equal(as.vector(csr %*% rep(1, length(mmm))), rep(0, nlevels(mmm)))
  
})
