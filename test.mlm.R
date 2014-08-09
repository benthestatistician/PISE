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
  res <- parseMatchingProblem(cost ~ pr*pt + mmm, n2)
  
  expect_is(res$fmla, "formula")
  expect_is(res$mf, "data.frame")
  expect_is(res$match, "optmatch")

  expect_equivalent(res$match, mmm)
  expect_equal(dim(res$mf)[2], 3) # cost, pr, pt, but no mmm

  mmm2 <- fullmatch(pr ~ t1, data = nuclearplants)
  n3 <- cbind(n2, mmm2)
  res2 <- parseMatchingProblem(cost ~ pr + mmm, n3)
  expect_equivalent(res2$match, mmm)
  expect_error(parseMatchingProblem(cost ~ pr + mmm + mmm2, n3), "one")
               
  expect_error(parseMatchingProblem(cost ~ pr, n3), "include")
})

test_that("optmatch -> matrix.csr", {
  data(nuclearplants)

  ### start with the most straightforward case: everyone is matched
  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)

  csr <- as(mmm, "matrix.csr")

  expect_is(csr, "matrix.csr")
  expect_equal(dim(csr), c(nlevels(mmm), length(mmm)))

  csrm <- as.matrix(csr) # cast it back to dense matrix for testing
  grps <- apply(csrm, 2, function(col) {
    tmp <- which(col != 0)
    if (length(tmp) == 0) {
      return(NA)
    }
    return(tmp)
  })

  tmp <- table(grps, mmm)
  expect_true(all(diag(tmp) != 0))
  diag(tmp) <- 0
  expect_true(all(tmp == 0))

  expect_true(all(rowSums(csrm) == 0))

  expect_equal(as.vector(csr %*% rep(1, length(mmm))), rep(0, nlevels(mmm)))

  ### next case: not everyone is match
  mmm <- pairmatch(pr ~ t1 + t2 + cap, data = nuclearplants)

  csr <- as(mmm, "matrix.csr")

  expect_is(csr, "matrix.csr")
  expect_equal(dim(csr), c(nlevels(mmm), length(mmm)))

  csrm <- as.matrix(csr) # cast it back to dense matrix for testing
  grps <- apply(csrm, 2, function(col) {
    tmp <- which(col != 0)
    if (length(tmp) == 0) {
      return(NA)
    }
    return(tmp)
  })

  tmp <- table(grps, mmm)
  expect_true(all(diag(tmp) != 0))
  diag(tmp) <- 0
  expect_true(all(tmp == 0))

  expect_true(all(rowSums(csrm) == 0))

  expect_equal(as.vector(csr %*% rep(1, length(mmm))), rep(0, nlevels(mmm)))

  ### knock a specific group
  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  mmm <- mmm[mmm != levels(mmm)[1]]

  csr <- as(mmm, "matrix.csr")

  expect_is(csr, "matrix.csr")
  expect_equal(dim(csr), c(nlevels(mmm), length(mmm)))

  csrm <- as.matrix(csr) # cast it back to dense matrix for testing
  grps <- apply(csrm, 2, function(col) {
    tmp <- which(col != 0)
    if (length(tmp) == 0) {
      return(NA)
    }
    return(tmp)
  })

  tmp <- table(grps, mmm)
  expect_true(all(tmp[, 1] == 0))
  expect_true(all(diag(tmp[,-1]) != 0))

  expect_true(all(rowSums(csrm) == 0))

  expect_equal(as.vector(csr %*% rep(1, length(mmm))), rep(0, nlevels(mmm)))
  
  ### remove only a treated member
  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  mmm <- mmm[!(mmm == levels(mmm)[1] & nuclearplants$pr == 1)]

  csr <- as(mmm, "matrix.csr")

  expect_is(csr, "matrix.csr")
  expect_equal(dim(csr), c(nlevels(mmm), length(mmm)))

  csrm <- as.matrix(csr) # cast it back to dense matrix for testing
  grps <- apply(csrm, 2, function(col) {
    tmp <- which(col != 0)
    if (length(tmp) == 0) {
      return(NA)
    }
    return(tmp)
  })

  tmp <- table(grps, mmm)
  expect_true(all(tmp[, 1] == 0))
  expect_true(all(diag(tmp[,-1]) != 0))

  expect_true(all(rowSums(csrm) == 0))

  expect_equal(as.vector(csr %*% rep(1, length(mmm))), rep(0, nlevels(mmm)))
})

test_that("mlm", {
  
  data(nuclearplants)
  mmm <- fullmatch(pr ~ t1 + t2 + cap, data = nuclearplants)

  n2 <- cbind(nuclearplants, mmm)

  # these are non-failure tests. just checking for errors
  # fails: mlm(cost ~ mmm, data = n2)
  mlm(cost ~ t1 + t2 + mmm, data = n2)
  mlm(cost ~ t1 + t2 + mmm, data = n2, fit.type = "robust")
  mlm(cost ~ t1 + t2 + mmm, data = n2, ms.weights = harmonic)
 
  ppp <- pairmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  n3 <- cbind(nuclearplants, ppp)

  mlm(cost ~ t1 + ppp, data = n3)
  mlm(cost ~ t1 + ppp, data = n3, fit.type = "robust")
  mlm(cost ~ t1 + ppp, data = n3, ms.weights = harmonic)
})

test_that("summary.lm extracts proper SEs", {
  ppp <- pairmatch(pr ~ t1 + t2 + cap, data = nuclearplants)
  n3 <- cbind(nuclearplants, ppp)

  expect_equal(coef(summary(mlm(cost ~ pr + ppp, data = n3)))[,],
               coef(summary(lm(cost ~ pr + ppp, data = n3)))["pr",])
})

