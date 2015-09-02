source("ppse.R")
stopifnot(require("testthat"))
data(nuclearplants, package="optmatch")
aglm <- glm(pr~.-cost, data=nuclearplants, family=binomial)
expect_false(is.null(aglm$qr))
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

test_that("ppse() recognizes simplify argument",
          {
            expect_true(inherits(ppse(aglm, simplify=FALSE), "list"))
            expect_true(!inherits(ppse(aglm, simplify=TRUE), "list"))
            expect_true(inherits(ppse(aglm$qr, fitted.model=aglm, simplify=FALSE), "list"))
            expect_true(!inherits(ppse(aglm$qr, fitted.model=aglm, simplify=TRUE), "list"))
          })

test_that("glm.fit's non-update of weights at last stage",
### this documents the reason for `ppse.qr` not to just pull weights from the fitted model object
          {
              expect_false(isTRUE(all.equal(getglmQweights(aglm$linear.predictor,family=aglm$family),
                                            aglm$weights, check.attributes=FALSE)))
              expect_true(all(abs(getglmQweights(aglm$linear.predictor,family=aglm$family) - aglm$weights) < 1e-5))
          })

test_that("`getglmQweights()` mimics `glm.fit()`'s creation of quadratic weights", {
    aglm0 <- update(aglm, control=list(maxit=(aglm$iter-1)))
    expect_false(aglm0$converged)
    aglm1 <- update(aglm0, etastart=aglm0$linear.predictors, control=list(maxit=1))
    expect_true(aglm1$converged)
    expect_true(isTRUE(all.equal(getglmQweights(aglm0$linear.predictor,family=aglm0$family),
                                 aglm1$weights, check.attributes=FALSE))
                )
### But it doesn't work out if you take weights and coeff's from the same iteration (even the last one)
    expect_false(isTRUE(all.equal(getglmQweights(aglm1$linear.predictor,family=aglm1$family),
                                 aglm1$weights, check.attributes=FALSE))
                 )
})

test_that("updating a converged glm makes no discernible change to its linear predictors",
          {
              expect_true(isTRUE(all.equal(predict(aglm),predict(update(aglm,etastart=predict(aglm))) )))
              ## Numerically speaking, coefficients converge before the weight vector does.
              ## Not sure I should have been surprised to find this, but I was;
              ## so I'll digress briefly to document it.
              aglm2 <- update(aglm, etastart=aglm$linear.predictors)
              expect_true(isTRUE(all.equal(coef(aglm2), coef(aglm)))) #so coeff's have converged,
              expect_true(isTRUE(all.equal(predict(aglm2), predict(aglm)))) # as have the eta's...
              expect_false(isTRUE(all.equal(aglm2$weights, aglm$weights))) #but fitting weights have not
              expect_false(isTRUE(all.equal(vcov(aglm2), vcov(aglm)))) #nor has the information matrix.
              expect_true((aglm2$deviance-aglm$deviance)^2 < #Differences may show up in 
                          .Machine$double.eps)               # deviances -- but just barely.

          })

test_that("Crude alignment between QR decomp-based calcs and coeffs, cov-hats associated with fitted glms",
          {
              ## the `C_Cdqrls`-based coefficent estimates don't align precisely w/ QR-based reconstruction
              expect_error(ppse(aglm$qr, fitted.model=aglm,
                                tol.coeff.alignment=1e-7), "reported coefficients differ by up to")
              ## but with fuzzy goggles they do look the same
              expect_that(ppse(aglm, tol.coeff.alignment=1e-5), is_a("numeric"))
              ## On the other hand, nominal cov-hat's do line up to within numerical tolerance
              expect_true(isTRUE(all.equal(qr.R(stats:::qr.lm(aglm)) %*%
                                           vcov(aglm) %*% # this is basically `chol2inv(qr.R(aglm$qr))`, via summary.glm
                                           t(qr.R(stats:::qr.lm(aglm))),
                                           diag(stats:::qr.lm(aglm)$rank), # under the Q basis, C-hat is the identity
                                           check.attributes=FALSE))
                          )
          })


test_that("Crude agreement between ppse.qr & ppse.glm",{
    expect_true(ppse(aglm)!=ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=F))
    expect_true(ppse(aglm)!=ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=T))
    expect_true(abs(ppse(aglm) - ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=F))<1e-5)
    expect_true(abs(ppse(aglm) - ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=T))<1e-5)
    expect_error(ppse(aglm$qr, covariance.extractor=sandwich, fitted.model=aglm), # Change this if/when further vcov support added
                 "covariance.extractor")
})

test_that("appropriately handles NA coefs",
          {
              aglm.alt <- update(aglm, formula=update(formula(aglm), .~.+factor(ne)))
              expect_true(any(is.na(coef(aglm.alt))))
              expect_equal(ppse(aglm), ppse(aglm.alt))
              expect_equal(ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=F),
                           ppse(aglm.alt$qr, fitted.model=aglm.alt, coeffs.from.fitted.model=F))
              expect_equal(ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=T),
                           ppse(aglm.alt$qr, fitted.model=aglm.alt, coeffs.from.fitted.model=T))
           })

test_that("appropriately handles strata in formula", {
    expect_true(require("survival"))
    aglm.s <- update(aglm, formula=update(formula(aglm),
                                        #.~.-ne+survival:::strata(ne))) # ppse doesn't recognize now;
                               .~.-ne+strata(ne))) #thus the `require("survival")` above.  Oughta fix that...
    expect_true(ppse(aglm.s) < ppse(aglm))
} )

test_that("redo_qr",
          {
              expect_true(abs(ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=F) -
                                  ppse(redo_qr(aglm, LAPACK=F), fitted.model=aglm, coeffs.from.fitted.model=F))
                          < 1e-5)
              expect_true(abs(ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=F) -
                                  ppse(redo_qr(aglm, LAPACK=T), fitted.model=aglm, coeffs.from.fitted.model=F))
                          < 1e-5)
              expect_true(abs(ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=T) -
                                  ppse(redo_qr(aglm, LAPACK=F), fitted.model=aglm, coeffs.from.fitted.model=F))
                          < 1e-5)
              expect_true(abs(ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=T) -
                                  ppse(redo_qr(aglm, LAPACK=T), fitted.model=aglm, coeffs.from.fitted.model=F))
                          < 1e-5)

          })
## (this test currently fails...)
##test_that("deals with fitted cox models too",
##          {
##test1 <- list(time=c(4,3,1,1,2,2,3), 
##              status=c(1,1,1,0,1,1,0), 
##              x=c(0,2,1,1,1,0,0), 
##              sex=c(0,0,0,0,1,1,1)) 
##aph <- coxph(Surv(time, status) ~ x + strata(sex), test1)
##expect_that(ppse(aph), is_a("numeric"))
##} )


