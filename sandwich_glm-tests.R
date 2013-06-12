source("sandwich_glm.R")
require(optmatch)
require(SparseM)
data(nuclearplants)
aPS <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
dim(GetInvObservedInfo(aPS))
as.matrix(SparseMMFromFactor(as.factor(ifelse(as.logical(nuclearplants$pt), 'yespt', 'nopt'))))
dim(EmpiricalScoreCov(aPS))

sqrt(diag(vcov(aPS)))
sqrt(diag(SandwichVcov(aPS)))

fm1 <- fullmatch(mdist(aPS))
sqrt(diag(SandwichVcov(aPS, fm1)))

fm2 <- fullmatch(mdist(aPS) + caliper(.2, aPS ) )
sqrt(diag(SandwichVcov(aPS, fm2)))
