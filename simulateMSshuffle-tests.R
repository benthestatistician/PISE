source("simulateMSshuffle.R")
source("sumStatStuff.R")
source("cumulant.R")
require('optmatch')
require('RItools')
data(nuclearplants)

test <- function(t, m = "Error!") {
  if (!t) {
    stop(m)
  }
}
test.eq <- function(x,y, digits=12, m = "Error!") {
  if (any(abs(x-y)>= 10^(-digits))) { stop(m) }
}


shouldError <- function(expr, msg = "Exception should be thrown") {
  r <- try(expr, silent = T)
  if (!inherits(r, "try-error")) {
    stop(msg)
  }
}


smstest <- list()
aGlm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
smstest$ps3 <- mdist(aGlm, structure.fmla=~pt)
fm3 <- fullmatch(smstest$ps3)
summary(fm3)

### Tests of makeAssignmentRecipe
smstest$ar1 <- makeAssignmentRecipe(aGlm$linear.predictors, fm3)
smstest$ar2a <- makeAssignmentRecipe(aGlm, fm3)
test(identical(smstest$ar1, smstest$ar2a))
names(attributes(smstest$ar1))
test(with(smstest, !any(is.na(as.integer(names(ar1[[1]]))))))
(smstest$notamanyone <- which.min(attr(smstest$ar1,"manyone")))
with(smstest,
     test.eq(aGlm$linear[as.integer(names(ar1[[notamanyone]]))],
        ar1[[notamanyone]] +
        attr(ar1, "matched.set.mean.linear.propensities")[notamanyone]
             )
     )

(smstest$amanyone <- which.max(attr(smstest$ar1,"manyone")))
with(smstest,
     test.eq(aGlm$linear[as.integer(names(ar1[[amanyone]]))],
        -ar1[[amanyone]] +
        attr(ar1, "matched.set.mean.linear.propensities")[amanyone]
             )
     )

smstest$ar2b <- makeAssignmentRecipe(aGlm, fm3, truncate.linear.ps.at=1)
smstest$ar0 <- makeAssignmentRecipe(NULL, fm3)

### Tests of msShuffler
(smstest$shuff0 <- msShuffler(smstest$ar0))
msShuffler(fm3)
(smstest$shuff1 <- msShuffler(smstest$ar1))
(smstest$shuff2b <- msShuffler(smstest$ar2b))
smstest$shuff2 <- msShuffler(smstest$ar1, shuffles=100)

### Tests of tilt.from.permdist
test.eq(tilt.from.permdist(smstest$shuff0), 1)
tilt.from.permdist(smstest$shuff1)
tilt.from.permdist(smstest$shuff2b) # should be less extreme than previous
### rdist
mystat <- function(z, tes=seq(-100, 50, by=50))
{ n <- nrow(nuclearplants)
  R <- with(nuclearplants,
            resid(lm(cost - pr*matrix(tes, n, length(tes), byrow=TRUE,
                        dimnames=list(rownames(nuclearplants),
                        gsub("-", "m", paste("te",tes,sep=".")) )
                             ) ~ date+cap+pt)
                  )
            )
                  coef(lm(R ~z+fm3, data=nuclearplants))['z',]
}

(rdObj1 <- rdist(smstest$shuff1,statistic.=mystat))

  R <- with(nuclearplants,
            resid(lm(cost - pr*matrix(seq(-100, 50, by=50),
                             nrow(nuclearplants), 4, byrow=TRUE,
                             dimnames=list(rownames(nuclearplants),
                               gsub("-", "m", paste("te",seq(-100, 50, by=50),sep="."))
                               )
                                      )
                     ~date+cap+pt))
            )

(xb1 <- xBalance(pr~ te.m100+ te.m50+te.0+te.50,
                 strata=fm3, data=data.frame(nuclearplants, R),
                 report=c('adj.mean.diffs','p.values')) )

estimateMean(rdObj1, FUN=function(x) rbind(EX=x, EX2=x^2, pv.2s=abs(x)>=abs(xb1$results[, 'adj.diff',])))

estimateMean(rdObj1, FUN=function(x) rbind(EX=x, EX2=x^2, pv.2s=abs(x)>=abs(xb1$results[, 'adj.diff',])), propensity=0)

rdObj2 <- rdist(smstest$shuff2,statistic.=mystat)
estimateSD(rdObj2, FUN=function(x) x, propensity=0) # Hmm, errors far too small -- fix someday?  

sapply(as.data.frame(R), function(x) sqrt(sumStatVar(x, attr(fm3, "contrast.group"), fm3))/
  (attr(stratumStructure(fm3),"comparable.num.matched.pairs")/2) )

estimateMean(rdObj2, FUN=function(x) rbind(EX=x, pv.2s=abs(x)>=abs(xb1$results[, 'adj.diff',])))
Eadjmeandiff(as.data.frame(R), fm3, aGlm$linear.predictors)


### Some trickery to get full permutation distribution:
smstest$shuff.all <- as.data.frame(t(expand.grid(lapply(smstest$ar0, function(x) as.integer(names(x)) ))))
names(smstest$shuff.all) <- 1L:length(smstest$shuff.all)

attributes(smstest$shuff.all) <- c(attributes(smstest$shuff.all),
                                   attributes(msShuffler(fm3))[c('class', 'names.all.units', 'contrast.group.template', 'manyone', 'mssize', 'contrast.group', 'positions.of.matched.units', 'matched.set.mean.linear.propensities', 'linear.propensities', 'matched.set.linear.propensities')])
rdObj.all <- rdist(smstest$shuff.all , statistic.=mystat)

estimateMean(rdObj.all, FUN=function(x) rbind(EX=x, pv.2s=abs(x)>=abs(xb1$results[, 'adj.diff',])))[[1]]
estimateMean(rdObj2, FUN=function(x) rbind(EX=x, pv.2s=abs(x)>=abs(xb1$results[, 'adj.diff',])), propensity=0)

all.equal(estimateMean(rdObj.all, FUN=function(x) x^2)[[1]], 
sapply(as.data.frame(R), function(x) sumStatVar(x, attr(fm3, "contrast.group"), fm3))/
  (attr(stratumStructure(fm3),"comparable.num.matched.pairs")/2)^2)

all.equal(estimateSD(rdObj.all)[[1]]*sqrt((length(rdObj.all$shuff)-1)/length(rdObj.all$shuff) ),
          sapply(as.data.frame(R), function(x) sqrt(sumStatVar(x, attr(fm3, "contrast.group"), fm3))/
          (attr(stratumStructure(fm3),"comparable.num.matched.pairs")/2))
          )
