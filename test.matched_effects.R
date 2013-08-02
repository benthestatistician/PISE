library(testthat)

context("matched_effects")

test_that("basic", {
  data(nuclearplants)

  f <- fullmatch(match_on(pr ~ cost, data=nuclearplants), data=nuclearplants)

  a <- matched_effects(nuclearplants$cost, f, weighting.scheme="h", keep.differences=FALSE)

  expect_equal(length(a), 2)
  expect_true(all(names(a) == c("effects", "covariance")))
  expect_true(is.numeric(a$effects))
  expect_equal(length(a$effects), 1)
  expect_true(is.matrix(a$covariance))
  expect_true(all(a$covariance[1,1] > 0))

  b <- matched_effects(nuclearplants[,c("cost","cap")], f, weighting.scheme="h", keep.differences=FALSE)

  expect_equal(length(b), 2)
  expect_true(all(names(b) == c("effects", "covariance")))
  expect_true(is.numeric(b$effects))
  expect_equal(length(b$effects), 2)
  expect_true(is.matrix(b$covariance))
  expect_true(all(a$covariance > 0))
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
  # here weighting is 1.333 & 1, so (1.333*7.5+10)/(2.333) ~= 60/7
  b <- matched_effects(resp, ff, "harmonic")
  expect_true(all.equal(b$eff,60/7, check.attributes=FALSE))



})
