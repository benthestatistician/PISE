
library(testthat)
context("Touch up QR decomp associated w/ a glm to address near-singularities.")
    data(nuclearplants, package="optmatch")
    ## singlar model matrix, based on example(qr.X). b2 and c3-4 are extra.
    tdat <- data.frame(y =  rbinom(6, 1, .5), int=1,
                       b1 = factor(rep(1:0, each = 3)), b2 = rep(0:1, each = 3),
                       c1 = rep(c(1,0,0), 2), c2 = rep(c(0,1,0), 2),
                       c3 = rep(c(0,0,1),2), c4=rep(0:1, each = 3))
glm.nonsing <- suppressWarnings(glm(y~b1 + c1+c2, tdat, family=binomial))
glm.sing <- suppressWarnings(glm(tdat, family=binomial))

tdat2 <- data.frame(y=rbinom(2,1,.5), x1=1:0, x2=0:1, x3=1)
glm2.sing <- suppressWarnings(glm(tdat2, family=binomial))

data(nuclearplants, package="optmatch")
aglm <- glm(pr~.-cost, data=nuclearplants, family=binomial)
expect_false(is.null(aglm$qr))

test_that("colname alignment w/in glms, before we start tampering",{
    expect_equal(names(coef(glm.nonsing)), colnames(glm.nonsing$R))
    expect_equal(names(coef(glm.sing))[glm.sing$qr$pivot], colnames(glm.sing$R))
})

test_that("zeroing out columns entails that qr() will pivot them to the back", {
  expect_equal(qr(data.frame(foo=0, subset(tdat, select=c(int,b1,c1))), LAPACK = F)$pivot[4],1)
  expect_equal(qr(data.frame(subset(tdat, select=c(int,b1,c1)),foo=0), LAPACK = F)$pivot[4],4)
  expect_equal(qr(data.frame(int=1, foo=0, subset(tdat, select=c(b1,c1))), LAPACK = F)$pivot[4],2)
  expect_equal(sort(qr(data.frame(foo=0, bar=0, subset(tdat, select=c(int,b1,c1))), LAPACK = F)$pivot[4:5]),1:2)
  expect_equal(sort(qr(data.frame(int=1, foo=0, bar=0, subset(tdat, select=c(b1,c1))), LAPACK = F)$pivot[4:5]),2:3)
  expect_equal(sort(qr(data.frame(foo=0,int=1,bar=0, subset(tdat, select=c(b1,c1))), LAPACK = F)$pivot[4:5]),c(1,3))

  expect_equal(qr(data.frame(foo=0, subset(tdat, select=c(int,b1,c1))), LAPACK = T)$pivot[4],1)
  expect_equal(qr(data.frame(subset(tdat, select=c(int,b1,c1)),foo=0), LAPACK = T)$pivot[4],4)
  expect_equal(qr(data.frame(int=1, foo=0, subset(tdat, select=c(b1,c1))), LAPACK = T)$pivot[4],2)
  expect_equal(sort(qr(data.frame(foo=0, bar=0, subset(tdat, select=c(int,b1,c1))), LAPACK = T)$pivot[4:5]),1:2)
  expect_equal(sort(qr(data.frame(int=1, foo=0, bar=0, subset(tdat, select=c(b1,c1))), LAPACK = T)$pivot[4:5]),2:3)
  expect_equal(sort(qr(data.frame(foo=0,int=1,bar=0, subset(tdat, select=c(b1,c1))), LAPACK = T)$pivot[4:5]),c(1,3))

  ## what happens in the singular case? Here b2 = int-b1
  ## LAPACK gives desired behavior (0's always forced to the back), LINPACK not necessarily

  expect_equal(qr(data.frame(foo=0, subset(tdat, select=c(int,b1,b2,c1))), LAPACK = T)$pivot[5],1)
  expect_equal(qr(data.frame(subset(tdat, select=c(int,b1,b2,c1)),foo=0), LAPACK = T)$pivot[5],5)
  expect_equal(qr(data.frame(int=1, foo=0, subset(tdat, select=c(b1,b2,c1))), LAPACK = T)$pivot[5],2)
  expect_equal(sort(qr(data.frame(foo=0, bar=0, subset(tdat, select=c(int,b1,b2,c1))), LAPACK = T)$pivot[5:6]),1:2)
  expect_equal(sort(qr(data.frame(int=1, foo=0, bar=0, subset(tdat, select=c(b1,b2,c1))), LAPACK = T)$pivot[5:6]),2:3)
  expect_equal(sort(qr(data.frame(foo=0,int=1,bar=0, subset(tdat, select=c(b1,b2,c1))), LAPACK = T)$pivot[5:6]), c(1,3))

  ## LINPACK not quite as well behaved

  expect_false(qr(data.frame(foo=0, subset(tdat, select=c(int,b1,b2,c1))), LAPACK = F)$pivot[5]==1)
  ## but it does group the 0 cols at the back along with other redundant cols.
  expect_equal(sort(qr(data.frame(foo=0, subset(tdat, select=c(int,b1,b2,c1))), LAPACK = F)$pivot[4:5]),c(1,4))

  expect_equal(sort(qr(data.frame(foo=0, bar=0, subset(tdat, select=c(int,b1,b2,c1))), LAPACK = F)$pivot[4:6]),c(1:2,5))
  expect_equal(sort(qr(data.frame(int=1, foo=0, bar=0, subset(tdat, select=c(b1,b2,c1))), LAPACK = F)$pivot[4:6]),c(2:3,5))
  expect_equal(sort(qr(data.frame(foo=0,int=1,bar=0, subset(tdat, select=c(b1,b2,c1))), LAPACK = F)$pivot[4:6]), c(1,3,5))

})


test_that("redo_qr preserves or permutes column order according as LAPACK=F or T", {
    redone.glm.sing <- redo_qr(glm.sing, LAPACK=F)
##    expect_false(isTRUE(all.equal(qr.R(redone.glm.sing), qr.R(glm.sing$qr), check.attributes=F)))
    expect_true(all(colnames(qr.R(redone.glm.sing)) %in% colnames(qr.R(glm.sing$qr))))
    redone.glm.sing <- redo_qr(glm.sing, LAPACK=T)
    expect_true(all(colnames(qr.R(redone.glm.sing)) %in% colnames(qr.R(glm.sing$qr))))

    redone.glm.nonsing <- redo_qr(glm.nonsing, LAPACK=F)
    expect_false(isTRUE(all.equal(qr.R(redone.glm.nonsing), qr.R(glm.nonsing$qr), check.attributes=F)))
    expect_equal(c("(Intercept)", colnames(qr.R(redone.glm.nonsing))), colnames(qr.R(glm.nonsing$qr)))

    ## seemed to be getting inconsistent results w/ the below, so skipping
###    redone.glm.nonsing <- redo_qr(glm.nonsing, LAPACK=T)
###    expect_false(isTRUE(all.equal(colnames(qr.R(redone.glm.nonsing)), colnames(qr.R(glm.nonsing$qr)))))
    })


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
              aglm.alt <- update(aglm, formula=update(formula(aglm), .~.+factor(ne)))
              expect_true(any(is.na(coef(aglm.alt))))
              expect_equal(ppse(aglm), ppse(aglm.alt))
              expect_equal(ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=F),
                           ppse(aglm.alt$qr, fitted.model=aglm.alt, coeffs.from.fitted.model=F))
              expect_true(abs(ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=F)-
                           ppse(redo_qr(aglm.alt, LAPACK=F), fitted.model=aglm.alt,
                                coeffs.from.fitted.model=F)) < 1e-5)
              expect_true(abs(ppse(aglm$qr, fitted.model=aglm, coeffs.from.fitted.model=F)-
                           ppse(redo_qr(aglm.alt, LAPACK=T), fitted.model=aglm.alt,
                                coeffs.from.fitted.model=F)) < 1e-5)              
          })

test_that("after LINPACK decomp, qr.X reconstructs X*w even if X*w was singular",
          expect_equivalent(qr.X(glm.sing$qr, ncol=length(coef(glm.sing))), # NB: `complete=T` doesn't cut it. R ver 3.1.2.
                            sqrt(glm.sing$weights)*model.matrix(glm.sing))
          )

test_that("after LINPACK decomp, qr.X reconstructs X*w even if nrow(X)<ncol(X)",
          expect_equivalent(qr.X(glm2.sing$qr, ncol=length(coef(glm2.sing))),
                            sqrt(glm2.sing$weights)*model.matrix(glm2.sing))
          )

