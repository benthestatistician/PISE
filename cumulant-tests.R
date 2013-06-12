require("optmatch")
data(nuclearplants)
source("sumStatStuff.R")
source("cumulant.R")
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

getH <- function(thematch)
  {
stopifnot(class(thematch)[1]=="optmatch")
thetab <- unclass(table(thematch, attr(thematch, "contrast.group")))
theHes <- apply(thetab, 1, function(x) 1/(1/x[1] + 1/x[2]))
sum(theHes)
  }

aps <- glm(pr ~ t1 + t2 + pt, family=binomial(), data=nuclearplants)
pairm <- pairmatch(mdist(aps))
fullm <- fullmatch(mdist(aps))

### Add tests of MScumfct and related

pair1 <- pairm=='m.1'
pair1 <- !is.na(pair1) &pair1
adf1 <- data.frame(TrTx=attr(pairm, "contrast.group")[pair1],
                             pSlIn=c(0,0), nuclearplants[pair1,  c("t1", "t2")])
test.eq(MScumfct(adf1, t=0, derivs=1), c(0,0))
test.eq(MScumfct(adf1, t=0, derivs=2), pair1var <- diag(var(nuclearplants[pair1,  c("t1", "t2")]))/2)

### Proper handling of throwaway case that m is 1
test.eq(MScumfct(adf1[1,], t=0, derivs=1), c(0,0))
test.eq(MScumfct(adf1[1,], t=0:1, derivs=1), data.frame(t1=c(0,0),t2=c(0,0)))

test(all.equal(pad.bottom.rows(adf=adf1, padto=3), data.frame(TrTx=c(T,F,F), rbind(adf1[-1], c(-Inf,0,0)))))
test(all.equal(pad.bottom.rows(adf=adf1, padto=4), data.frame(TrTx=c(T,F,F,F), rbind(adf1[-1], c(-Inf,0,0), c(-Inf,0,0)))))

test.eq(MScumfct(data.frame(TrTx=c(T,F,T), rbind(adf1[-1], c(Inf,0,0))), t=0, derivs=1), c(0,0))
test.eq(MScumfct(data.frame(TrTx=c(T,F,F), rbind(adf1[-1], c(-Inf,0,0))), t=0, derivs=1), c(0,0))
test.eq(MScumfct(data.frame(TrTx=c(T,F,T), rbind(adf1[-1], c(Inf,0,0))), t=0, derivs=2), pair1var)
test.eq(MScumfct(data.frame(TrTx=c(T,F,F), rbind(adf1[-1], c(-Inf,0,0))), t=0, derivs=2), pair1var)

### Testing MScumfct with list of pre-formed cumulant expressions
MSClists <- list()
MSClists$a <- MS.cumfct.expr.factory(m=c(2,5,10), derivs=1)
test.eq(MScumfct(adf1, t=0, derivs=MSClists$a), c(0,0))
test.eq(MScumfct(pad.bottom.rows(adf1, padto=4), t=0, derivs=MSClists$a), c(0,0))

MSClists$b <- MS.cumfct.expr.factory(m=c(2,5,10), derivs=2)
test.eq(MScumfct(adf1, t=0, derivs=MSClists$b), pair1var)
test.eq(MScumfct(pad.bottom.rows(adf1, padto=4), t=0, derivs=MSClists$b), pair1var)

### basic operation
cumulantfct(t=0, vars=nuclearplants$t1, thematch=pairm, pscore.linear=NULL, derivs=0)
test.eq(cumulantfct(t=0, vars=nuclearplants$t1, thematch=pairm, pscore.linear=NULL, derivs=1)[1,2],0)
test.eq(cumulantfct(t=0, vars=nuclearplants$t1, thematch=pairm, pscore.linear=NULL, derivs=2)[1,2],
     sumStatVar(nuclearplants$t1, attr(pairm, "contrast.group"), pairm)/(getH(pairm)^2)
     )

test.eq(cumulantfct(t=0, vars=nuclearplants$t1, thematch=fullm, pscore.linear=NULL, derivs=1)[1,2],0)

test.eq(cumulantfct(t=0, vars=nuclearplants$t1, thematch=fullm, pscore.linear=NULL, derivs=2)[1,2],
     sumStatVar(nuclearplants$t1, attr(fullm, "contrast.group"), fullm)/(getH(fullm)^2)
     )

### Not quite ready to take multiple derivs at once, but on its way.  Have to get the following 
### to spit it out in a proper form:
### cumulantfct(t=0, vars=nuclearplants$t1, thematch=fullm, pscore.linear=NULL, derivs=1:2)
###cumulantfct(t=0:1, vars=nuclearplants$t1, thematch=fullm, pscore.linear=NULL, derivs=1:2)
### cumulantfct(t=0:1, vars=nuclearplants[c('t1','t2')], thematch=fullm, pscore.linear=NULL, derivs=1:2)

### Takes data frame vars argument
test.eq(cumulantfct(t=0, vars=nuclearplants[c('t1','t2')], thematch=pairm, pscore.linear=NULL, derivs=1)[1,],rep(0, 2))
test.eq(cumulantfct(t=0, vars=nuclearplants[c('t1','t2')], thematch=pairm, pscore.linear=NULL, derivs=2)[1:2,2],
     c(sumStatVar(nuclearplants$t1, attr(pairm, "contrast.group"), pairm),
       sumStatVar(nuclearplants$t2, attr(pairm, "contrast.group"), pairm))/
        (getH(pairm)^2)
     )

### Takes vector t arguments
test.eq(cumulantfct(t=((-1:1)/2), vars=nuclearplants$t1, thematch=pairm, pscore.linear=NULL, derivs=1)[2,2],0)
test.eq(cumulantfct(t=((-1:1)/2), vars=nuclearplants[c('t1','t2')], thematch=pairm, pscore.linear=NULL, derivs=1)[2,2:3],rep(0,2))

### With non-null pscore argument, reproduces expectations calculated independently
test.eq(cumulantfct(t=0, vars=nuclearplants[c('t1','t2')], thematch=pairm, pscore.linear=aps$lin, derivs=1)[1:2,2],
        Eadjmeandiff(vars=nuclearplants[c('t1','t2')], thematch=pairm, pscore.linear=aps$lin) )

test.eq(cumulantfct(t=0, vars=nuclearplants[c('t1','t2')], thematch=fullm, pscore.linear=aps$lin, derivs=1)[1:2,2],
        Eadjmeandiff(vars=nuclearplants[c('t1','t2')], thematch=fullm, pscore.linear=aps$lin) )


### With matched pairs, untilted variance exceeds tilted variance

test(cumulantfct(t=0, vars=nuclearplants$t1, thematch=pairm, pscore.linear=NULL, derivs=2)[1,2]>
cumulantfct(t=0, vars=nuclearplants$t1, thematch=pairm, pscore.linear=aps$lin, derivs=2)[1,2])

### VSK function
VSKres.new <- list()
VSKres.new$pair.ps0.t1 <- VSK(NULL, pairm, nuclearplants$t1)
VSKres.new$pair.aps.t1 <- VSK(aps, pairm, nuclearplants$t1)

VSKres.new$pair.ps0.t1t2 <- VSK(NULL, pairm, nuclearplants[c('t1', 't2')])

VSKres.new$full.ps0.t1 <- VSK(NULL, fullm, nuclearplants$t1)
VSKres.new$full.aps.t1 <-VSK(aps, fullm, nuclearplants$t1)

load("cumulant-tests-results.RData")
for (nm in names(VSKres.new))
  {
    cat(paste("Testing ", nm, "...\n"))
    test.eq(VSKres.new[[nm]], VSKres.old[[nm]])
  }

VSKres.old <- VSKres.new
save(VSKres.old, file="cumulant-tests-results.RData")
