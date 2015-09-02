
library(testthat)
context("Touch up QR decomp associated w/ a glm to address near-singularities.")
    data(nuclearplants, package="optmatch")
    ## singlar model matrix, based on example(qr.X). b2 and c3-4 are extra.
    tdat <- data.frame(y =  rbinom(6, 1, .5), 
                       b1 = factor(rep(1:0, each = 3)), b2 = rep(0:1, each = 3),
                       c1 = rep(c(1,0,0), 2), c2 = rep(c(0,1,0), 2),
                       c3 = rep(c(0,0,1),2), c4=rep(0:1, each = 3))
glm.nonsing <- suppressWarnings(glm(y~b1 + c1+c2, tdat, family=binomial))
glm.sing <- suppressWarnings(glm(tdat, family=binomial))

tdat2 <- data.frame(y=rbinom(2,1,.5), x1=1:0, x2=0:1, x3=1)
glm2.sing <- suppressWarnings(glm(tdat2, family=binomial))

test_that("colname alignment w/in glms, before we start tampering",{
    expect_equal(names(coef(glm.nonsing)), colnames(glm.nonsing$R))
    expect_equal(names(coef(glm.sing))[glm.sing$qr$pivot], colnames(glm.sing$R))
})

test_that("check that redo_qr w/ LAPACK=F is a non-op", {
    redone.glm.nonsing <- redo_qr(glm.nonsing, LAPACK=F)
    expect_equivalent(redone.glm.nonsing, glm.nonsing)
    expect_equal(colnames(redone.glm.nonsing$R), colnames(glm.nonsing$R))
                                                         })

test_that("after LINPACK decomp, qr.X reconstructs X*w even if X*w was singular",
          expect_equivalent(qr.X(glm.sing$qr, ncol=length(coef(glm.sing))), # NB: `complete=T` doesn't cut it. R ver 3.1.2.
                            sqrt(glm.sing$weights)*model.matrix(glm.sing))
          )

test_that("after LINPACK decomp, qr.X reconstructs X*w even if nrow(X)<ncol(X)",
          expect_equivalent(qr.X(glm2.sing$qr, ncol=length(coef(glm2.sing))),
                            sqrt(glm2.sing$weights)*model.matrix(glm2.sing))
          )

test_that("reconstruction of singular & nonsingular X*w after the qr redo", {
    redone.glm.sing <- redo_qr(glm.sing, LAPACK=T)
    redone.glm.nonsing <- redo_qr(glm.nonsing, LAPACK=T)
    expect_equal(qr.X(redone.glm.sing$qr, ncol=length(coef(glm.sing))),
                 qr.X(glm.sing$qr, ncol=length(coef(glm.sing))))
    expect_equal(qr.X(redone.glm.nonsing$qr, ncol=length(coef(glm.nonsing))),
                 qr.X(glm.nonsing$qr, ncol=length(coef(glm.nonsing))))
} )

test_that("redo_qr tweaks the QR if X*w was singular", # something more descriptive of what we hope for?
          {
    expect_true(setequal(colnames(redone.glm.sing$R), colnames(glm.sing$R)))
    expect_false(all(colnames(redone.glm.sing$R)==colnames(glm.sing$R))) 
} )


ppty1 <- glm(pr ~ . - (pr + cost), family = binomial(), data = np.m)
    np.m <- transform(nuclearplants,
                      fuzzydate=date+rnorm(nrow(nuclearplants), sd=.Machine$double.eps^.5)
                      )
    ppty2 <- glm(pr ~ . - (pr + cost), family = binomial(), data = np.m)
    
