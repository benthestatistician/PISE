source("firstThreeMoments.R")
source("SparseM-based-tools.R")
require("optmatch")
data(nuclearplants)

aps <- glm(pr ~ t1 + t2 + pt, family=binomial(), data=nuclearplants)
pairm <- pairmatch(mdist(aps), data=nuclearplants)
fullm <- fullmatch(mdist(aps), data=nuclearplants)

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

load("cumulant-tests-results.RData")
first3Mres <- list()
first3Mres$pair.ps0 <- first3M(nuclearplants[c('t1', 't2')], pairm)
first3Mres$pair.aps <- first3M(nuclearplants[c('t1', 't2')], pairm, predict(aps))
first3Mres$full.ps0 <- first3M(nuclearplants[c('t1', 't2')], full)
first3Mres$full.aps <- first3M(nuclearplants[c('t1', 't2')], full, predict(aps))

for (nm in names(VSKres.new))
  {
    cat(paste("Testing ", nm, "...\n"))
    test.eq(first3Mres[[nm]],# transpose me?
            VSKres.old[[nm]][c("Exp", "Var", "Skew")])
  }
