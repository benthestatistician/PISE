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
  expect_is(res$match, "optmatch")

  expect_equal(res$match, mmm)
  expect_equal(dim(res$mf)[2], 3) # cost, pr, pt, but not mmm

  mmm2 <- fullmatch(pr ~ t1, data = nuclearplants)
  n3 <- cbind(n2, mmm2)
  expect_equal(parseMatchingProblem(cost ~ pr + mmm, n3)$match, mmm)
  expect_error(parseMatchingProblem(cost ~ pr + mmm + mmm2, n3), "one")
               
  expect_error(parseMatchingProblem(cost ~ pr, n3), "include")
})
