library(testthat)
library(optmatch)
source("findMatchedSetSingletons.R")

context("findMatchedSetSingletons")

test_that("basic", {
  data(nuclearplants)

  ### Example without many-one matched sets
  f <- fullmatch(match_on(pr ~ cost, data=nuclearplants), data=nuclearplants)
  fmss0 <- findMatchedSetSingletons(f)
  
  expect_equal(names(fmss0), c("theOne.name", "theOne.position","isOneMany"))
  expect_equal(sapply(fmss0, class), c("character","integer","logical"))
### now have to check that it makes sense...
}
          )
