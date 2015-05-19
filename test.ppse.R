data(nuclearplants, package="optmatch")
aglm <- glm(pr~.-cost, data=nuclearplants, family=binomial)
expect_that(ppse(aglm), is_a("numeric"))

## regression test
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

