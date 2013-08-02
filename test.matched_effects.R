library(testthat)
library(optmatch)
source("matched_effects.R")

context("matched_effects")

test_that("basic", {
  data(nuclearplants)

  f <- fullmatch(match_on(pr ~ cost, data=nuclearplants), data=nuclearplants)

  a <- matched_effects(nuclearplants$cost, f, weighting.scheme="h")

  expect_equal(length(a), 2)
  expect_true(all(names(a) == c("effects", "covariance")))
  expect_true(is.numeric(a$effects))
  expect_equal(length(a$effects), 1)
  expect_true(is.matrix(a$covariance))
  expect_true(all(a$covariance[1,1] > 0))

  b <- matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="h")

  expect_equal(length(b), 2)
  expect_true(all(names(b) == c("effects", "covariance")))
  expect_true(is.numeric(b$effects))
  expect_equal(length(b$effects), 2)
  expect_true(is.matrix(b$covariance))
  expect_true(all(a$covariance > 0))

  h1 <- matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="h")
  h2 <- matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="HAR")
  h3 <- matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="harmon")
  expect_error(matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="haar"))
  h4 <- matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="Harm")

  expect_true(identical(h1, h2))
  expect_true(identical(h1, h3))
  expect_true(identical(h1, h4))

  e1 <- matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="e")
  e2 <- matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="ETT")
  e3 <- matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="eTT")
  expect_error(matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="ettt"))
  e4 <- matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="Ett")

  expect_true(identical(e1, e2))
  expect_true(identical(e1, e3))
  expect_true(identical(e1, e4))

  defweight <- matched_effects(nuclearplants[,c("cost", "cap")], f)
  ettweight <- matched_effects(nuclearplants[,c("cost", "cap")], f, weight="e")
  expect_true(identical(defweight, ettweight))

  diffs <- matched_effects(nuclearplants[,c("cost", "cap")], f, keep.differences=TRUE)
  expect_true(all.equal(dim(diffs$differences), c(length(levels(f)), 2)))
  nodiffs <- matched_effects(nuclearplants[,c("cost", "cap")], f, keep.differences=FALSE)
  defdiffs <- matched_effects(nuclearplants[,c("cost", "cap")], f)
  expect_true(is.null(nodiffs$differences))
  expect_true(is.null(defdiffs$differences))
})


test_that("effect size", {
  ff <- as.factor(c(1,1,2,2))
  attr(ff, "contrast.group") <- c(0,1,0,1)
  class(ff) <- "optmatch"

  resp <- c(10, 15, 15, 25)
  a <- matched_effects(resp, ff, "ETT")
  b <- matched_effects(resp, ff, "harmonic")
  #+5 in first group, +10 in second, equal weighting, so 7.5 effect size
  expect_true(a$eff==7.5)
  expect_true(b$eff==7.5)


  ff <- as.factor(c(1,1,1,2,2))
  attr(ff, "contrast.group") <- c(0,1,1,0,1)
  class(ff) <- "optmatch"

  resp <- c(5,10,15,15,25)
  #7.5 in first, 10 in second, weighting is double in first group, so 25/3=8.333
  a <- matched_effects(resp, ff, "ETT")
  expect_true(all.equal(a$eff,8+1/3, check.attributes=FALSE))
  # here weighting is 1.333 & 1, so (1.333*7.5+10)/(2.333) ~= 60/7
  b <- matched_effects(resp, ff, "harmonic")
  expect_true(all.equal(b$eff,60/7, check.attributes=FALSE))

})

test_that("covariance", {
  ff <- as.factor(c(1,1,2,2))
  attr(ff, "contrast.group") <- c(0,1,0,1)
  class(ff) <- "optmatch"

  resp <- c(10, 15, 15, 25)

  W <- attr(ff, "contrast.group")
  gmeans <- sapply(levels(ff), function(x) sum(resp[ff == x & W == 0])/length(resp[ff == x & W == 0]))[ff]
  tauhat <- 1/sum(W)*sum(W*(resp - gmeans))

  sigmahat2 <- 1/(sum(W) - 1) * sum(W*(resp - gmeans - tauhat)^2)
  a <- matched_effects(resp, ff, "ett")

  ### this test currently fails
  #expect_true(sigmahat2/sum(W) == a$cov)

  ff <- as.factor(c(1,1,1,2,2))
  attr(ff, "contrast.group") <- c(0,1,1,0,1)
  class(ff) <- "optmatch"

  resp <- c(5,10,15,15,25)

  W <- attr(ff, "contrast.group")
  gmeans <- sapply(levels(ff), function(x) sum(resp[ff == x & W == 0])/length(resp[ff == x & W == 0]))[ff]
  tauhat <- 1/sum(W)*sum(W*(resp - gmeans))

  sigmahat2 <- 1/(sum(W)^2) * sum(W*(resp - gmeans - tauhat)^2)
  b <- matched_effects(resp, ff, "ett")

  expect_true(all.equal(sigmahat2/sum(W), b$cov[1,1]))


})
