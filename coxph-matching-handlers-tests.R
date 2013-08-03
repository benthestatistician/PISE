library(survival)
library(optmatch)
library(testthat)
source("coxph-matching-handlers.R")
context("Helpers for risk-set matching using survival:coxph")


example(cluster)

### need tests for exactMatch() applied to a Surv

matches <- models <- matchons <- list()


### Adapting some tests from coxph() examples
   test1 <- list(time=c(4,3,1,1,2,2,3), 
                   status=c(1,1,1,0,1,1,0), 
                   x=c(0,2,1,1,1,0,0), 
                   sex=c(0,0,0,0,1,1,1)) 
test2 <- within(test1, x <- c(NA, x[-1]))
test2 <- as.data.frame(test2)
row.names(test2)

expect_equal(ncol(convertSurv2to3(with(test1, Surv(time, status)))),3) 

matchons$test1a <- exactMatch(with(test1, Surv(time, status)))
expect_is(matchons$test1a, "InfinitySparseMatrix")

matchons$test2a <- exactMatch(with(na.omit(test2), Surv(time, status)), row.names=rownames(na.omit(test2))) 
expect_is(matchons$test2a, "InfinitySparseMatrix")
expect_identical(as.matrix(matchons$test1a)[-1,], as.matrix(matchons$test2a))

matchons$test1b <- match_on( coxph(Surv(time, status) ~ x, test1) )
expect_is(matchons$test1b, "InfinitySparseMatrix")
expect_equal(dim(matchons$test1a), dim(matchons$test2a) + 1:0) # Got to here

### With strata(), match_on() should only entertain within-stratum pairings
matchons$test1c <- match_on( coxph(Surv(time, status) ~ x + strata(sex), test1) )

### Need a test for standardization.scale, with and witout a strata term in the model.

### Confirm handling of within= arg
em.sex <- exactMatch(status~sex, data=test1)
matchons$test1d <- match_on( coxph(Surv(time, status) ~ x, test1), within=em.sex )
with(matchons, expect_equal(test1c < Inf, test1d < Inf))
### Now adapting from cluster() examples
###should throw warning, cluster() or frailty() term detected, maybe you wanted to give a "units" argument
test_that(matchons$cluster1 <- match_on(coxph(Surv(time, status) ~ rx + cluster(litter), rats)), gives_warning())
test_that(matchons$frailty1 <- match_on(coxph(Surv(time, status) ~ rx + frailty(litter), rats)), gives_warning())

### Test of ISM method for aggregate -- to be written


### Calipers.  match_on's caliper= option not working at present; so use caliper()
### optmatch doesn't seem to unit test equivalence of use of caliper= arg w/ results of caliper(),
### so I'm suspicious that the problem may be there.
models$heart1 <- coxph(Surv(start, stop, event) ~ age + year + transplant + cluster(id), data=heart)
matchons$heart1 <- match_on(models$heart1)
coxps <- predict(models$heart1)
names(coxps) <- row.names(heart)
w1 <- with(heart, exactMatch(Surv(start, stop, event)))
expect_equal(as.matrix(matchons$heart1),
             as.matrix(match_on(coxps, z=heart$event,
                                       within=w1)))
matchons$heart.cal1 <- caliper(matchons$heart1, 1)
table(apply(as.matrix(matchons$heart.cal1), 1, function(x) any(is.finite(x))))
table(apply(as.matrix(matchons$heart.cal1), 2, function(x) any(is.finite(x))))
expect_that(with(matchons, all(as.matrix(heart1)[is.finite(as.matrix(heart.cal1))] <= 1)), is_true())

matchons$heart1.w.cal0 <- match_on(coxps, z=heart$event, within=w1, caliper=1)
expect_that(with(matchons, all(as.matrix(heart1.w.cal0)[is.finite(as.matrix(heart1.w.cal0))] <= 1)), is_true())

matchons$heart1.w.cal1 <- match_on(models$heart1, caliper=1) # bombs here.  how come?  bug appears to be local
expect_that(with(matchons, all(as.matrix(heart1.w.cal1)[is.finite(as.matrix(heart1.w.cal1))] <= 1)), is_true())

summary(with(matchons, as.matrix(heart1.w.cal1)[is.finite(as.matrix(heart1.w.cal1))]))

### Does it balance?
matchons$heart1.w.cal1 <- with(matchons, heart1 + heart.cal1)
matches$heart1.w.cal1 <- fullmatch(matchons$heart1.w.cal1, data=heart)
expect_equal(as.logical(heart$event), attr(matches$heart1.w.cal1, "contrast.group"))

matchons$heart1.w.cal05 <- matchons$heart1 + caliper(matchons$heart1, .5)
matches$heart1.w.cal05 <- fullmatch(matchons$heart1.w.cal05, data=heart)
### summary.optmatch not equipped to handle wacky coxph nonstorage of data

library(RItools)
xBalance(event~age + year + transplant,
         strata=list(unstrat=NULL,
           cal0.5=~matches$heart1.w.cal05,
           cal1.0=~matches$heart1.w.cal1),
         data=heart,
         report = c('adj.means', 'z.scores', 'chisquare.test'))
