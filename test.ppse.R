source("ppse.R")
stopifnot(require("testthat"))
data(nuclearplants, package="optmatch")
aglm <- glm(pr~.-cost, data=nuclearplants, family=binomial)
expect_false(is.null(aglm$qr))
expect_that(ppse(aglm), is_a("numeric"))

test_that("Gives same answer as original version", {
gpsc <- function(propensity.glm, covariance.extractor=vcov) {

  stopifnot(inherits(propensity.glm, "glm"))

  data <- model.matrix(propensity.glm)
  data <- data[,!colnames(data) == "(Intercept)", drop = FALSE]

  covx <- cov(data)

  covb <- covariance.extractor(propensity.glm)

  covb <- covb[,!colnames(covb) == "(Intercept)", drop = FALSE]
  covb <- covb[!rownames(covb) == "(Intercept)",, drop = FALSE]

  # calculate the correction for the expected difference in propensity scores
  ps <- predict(propensity.glm)
  rho.beta <- apply(data, 2, function(xi) { cor(xi, ps) })
  sd.x <- sqrt(diag(covx))

  srb <- sd.x * rho.beta
 
  sqrt(2 * sum(covx * covb) - 2 * (t(srb) %*% covb %*% srb)[1,1]) 
}
## here are the regression tests
expect_equal(ppse_notstabilized(aglm), gpsc(aglm))
expect_true( abs(ppse(aglm) - gpsc(aglm)) < 1e-5 )
} )
test_that("ppse() recognizes simplify argument",
          {
            expect_true(inherits(ppse_notstabilized(aglm, simplify=FALSE), "list"))
            expect_true(!inherits(ppse_notstabilized(aglm, simplify=TRUE), "list"))
            expect_true(inherits(ppse(aglm, simplify=FALSE), "list"))
            expect_true(!inherits(ppse(aglm, simplify=TRUE), "list"))
          })

test_that("glm.fit's non-update of weights at last stage",
### this documents the reason for `ppse.qr` not to just pull weights from the fitted model object
          {
              expect_false(isTRUE(all.equal(getglmQweights(aglm$linear.predictor,family=aglm$family),
                                            aglm$weights, check.attributes=FALSE)))
              expect_true(all(abs(getglmQweights(aglm$linear.predictor,family=aglm$family) - aglm$weights) < 1e-5))
          })

test_that("`getglmQweights()` mimics `glm.fit()`'s creation of quadratic weights", {
    aglm0 <- suppressWarnings(update(aglm, control=list(maxit=(aglm$iter-1))))
    expect_false(aglm0$converged)
    aglm1 <- suppressWarnings(update(aglm0, etastart=aglm0$linear.predictors, control=list(maxit=1)))
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
              expect_error(ppse(aglm,
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


test_that("Crude agreement between ppse__notstabilized & ppse_via_qr",{
    expect_true(ppse_notstabilized(aglm)!=ppse(aglm, coeffs.from.fitted.model=F))
    expect_true(ppse_notstabilized(aglm)!=ppse(aglm, coeffs.from.fitted.model=T))
    expect_true(abs(ppse_notstabilized(aglm) - ppse(aglm, coeffs.from.fitted.model=F))<1e-5)
    expect_true(abs(ppse_notstabilized(aglm) - ppse(aglm, coeffs.from.fitted.model=T))<1e-5)
    expect_true(abs(ppse_notstabilized(aglm, covariance.estimator="sandwich") -
                        ppse(aglm, covariance.estimator="sandwich")) <1e-5)
})

test_that("Numerical stabilization from ppse.qr",{
    ## Josh E example in which ppse.glm doesn't work.
    set.seed(201510)
    n <- 100
    d <- data.frame(y=sample(0:1, n, TRUE),
                    x1=rnorm(n),
                    x2=rnorm(n))
    d$x3 <- d$x1+d$x2+rnorm(n,sd=.0000001)
    mod1 <- glm(y~., data=d, family=binomial)
    ppse1_nostab <- ppse_notstabilized(mod1, "sandwich", simplify=FALSE)
    expect_true(sum(ppse1_nostab$cov.beta *
                    makeSperp(ppse1_nostab$cov.X, ppse1_nostab$betahat))<0)
    ppse1_qr <- ppse(mod1, "sandwich", simplify=FALSE)
    expect_true(sum(ppse1_qr$cov.beta *
                    makeSperp(ppse1_qr$cov.X, ppse1_qr$betahat)
                    )>0)
})
test_that("appropriately handles NA coefs",
          {
              aglm.alt <- update(aglm, formula=update(formula(aglm), .~.+factor(ne)))
              expect_true(any(is.na(coef(aglm.alt))))
              expect_equal(ppse_notstabilized(aglm), ppse_notstabilized(aglm.alt))
              expect_equal(ppse(aglm, coeffs.from.fitted.model=F),
                           ppse(aglm.alt, coeffs.from.fitted.model=F))
              expect_equal(ppse(aglm, coeffs.from.fitted.model=T),
                           ppse(aglm.alt, coeffs.from.fitted.model=T))
           })

test_that("ppse_notstabilized appropriately handles strata in formula", {
    expect_true(require("survival"))
    aglm.s <- update(aglm, formula=update(formula(aglm),
                                        #.~.-ne+survival:::strata(ne))) # ppse doesn't recognize now;
                               .~.-ne+strata(ne))) #thus the `require("survival")` above.  Oughta fix that...
    expect_true(ppse_notstabilized(aglm.s) < ppse_notstabilized(aglm))
    expect_equal(ppse_notstabilized(aglm.s, terms.to.sweep.out=NULL),
                 ppse_notstabilized(aglm))
} )

test_that("ppse() handles strata() terms *if* placed next to Intercept", {
    expect_true(require("survival"))
    aglm.s <- update(aglm, formula=update(formula(aglm),
                               .~.-ne+strata(ne))) 
    expect_warning(ppse(aglm.s), "ignored")
    expect_equal(ppse(aglm.s), ppse(aglm) )
    aglm.s2 <- update(aglm, formula=update(formula(aglm),
                               .~strata(ne)+.-ne)) 
    expect_equal(ppse(aglm.s2, terms.to.sweep.out=NULL),
                 ppse(aglm) )
    expect_true(ppse(aglm.s2) < ppse(aglm))
} )

test_that("accepts arm::bayesglm fits", {
    expect_true(require("arm"))
    bglm <- bayesglm(pr~.-cost, data=nuclearplants, family=binomial)
    expect_false(is.null(bglm$qr))
    expect_that(ppse_notstabilized(bglm), is_a("numeric"))
    expect_true(inherits(ppse_notstabilized(bglm, simplify=FALSE), "list"))
    expect_error(ppse(bglm, # have yet to teach it how to rebuild coefs 
                      coeffs.from.fitted.model=FALSE)) # from QR
    expect_true(inherits(ppse(bglm, simplify=FALSE,
                              coeffs.from.fitted.model=TRUE), ## with this it ought to go through
                         "list"))
    expect_equal(ppse(bglm), ppse(bglm, coeffs.from.fitted.model=TRUE))
    expect_true(inherits(ppse(bglm, covariance.estimator="sandwich",
                              simplify=FALSE,
                              coeffs.from.fitted.model=TRUE), 
                         "list"))
})




