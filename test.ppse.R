source("ppse.R")
stopifnot(require("testthat"))
data(nuclearplants, package="optmatch")
aglm <- glm(pr~.-cost, data=nuclearplants, family=binomial)

expect_that(ppse(aglm), is_a("numeric"))

test_that("Gives same answer as original version", {
gpsc <- function(propensity.glm) {

  stopifnot(inherits(propensity.glm, "glm"))

  data <- model.matrix(propensity.glm)
  data <- data[,!colnames(data) == "(Intercept)", drop = FALSE]

  covx <- cov(data)

  # wrapping in a tryCatch as we observed some issues with bayesglm and sandwich.
  covb <- tryCatch(sandwich(propensity.glm),
                 error = function(e) { vcov(propensity.glm) })

  covb <- covb[,!colnames(covb) == "(Intercept)", drop = FALSE]
  covb <- covb[!rownames(covb) == "(Intercept)",, drop = FALSE]

  # could we ever get a covb with different dim than covx? is this worth checking?

  # calculate the correction for the expected difference in propensity scores
  ps <- predict(propensity.glm)
  rho.beta <- apply(data, 2, function(xi) { cor(xi, ps) })
  sd.x <- sqrt(diag(covx))

  srb <- sd.x * rho.beta
 
  sqrt(2 * sum(covx * covb) - 2 * (t(srb) %*% covb %*% srb)[1,1]) # indexing to turn a 1x1 matrix into a scalar
}
## here's the regression test
expect_equal(ppse(aglm), gpsc(aglm))
}
          )

test_that("glm.fit's non-update of weights at last stage",
          {
              expect_equal(FALSE,
                           isTRUE(all.equal(getglmQweights(aglm$linear.predictor,family=aglm$family),
                                            aglm$weights)))
              expect_true(all(abs(getglmQweights(aglm$linear.predictor,family=aglm$family) - aglm$weights) < 1e-5))
              
          }) #(documenting the reason not to just pull weights from the fitted model obj)

test_that("Agreement between ppse.qr & ppse.glm.",{
    expect_equal(ppse(aglm), ppse(aglm$qr, fitted.model=aglm))
    expect_equal(ppse(aglm), ppse(aglm$qr, data=model.frame(aglm), fitted.model=aglm))
    expect_error(ppse(aglm$qr, fitted.model=aglm), # Change this if/when further vcov support added
                 "only supports vcov as covariance extractor")
})
test_that("appropriately handles strata in formula", {
    aglm.s <- update(aglm, formula=update(formula(aglm), .~.-ne+strata(ne)))
    expect_true(ppse(aglm.s) < ppse(aglm))
} )

## (this test currently fails...)
test_that("deals with fitted cox models too",
          {
test1 <- list(time=c(4,3,1,1,2,2,3), 
              status=c(1,1,1,0,1,1,0), 
              x=c(0,2,1,1,1,0,0), 
              sex=c(0,0,0,0,1,1,1)) 
aph <- coxph(Surv(time, status) ~ x + strata(sex), test1)
expect_that(ppse(aph), is_a("numeric"))
} )


