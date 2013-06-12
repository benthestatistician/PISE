require(optmatch)
data(nuclearplants)
source("matchedDifferencesFromOptmatch.R")
source("sandwich_glm.R")
pm1 <- pairmatch( mdist(pr ~ t1 + t2, data =nuclearplants) )

MDMp1 <- MatchedDiffMaker(pm1, type="dense")
stopifnot(round(coef(lm(cost~pr+pm1, data=nuclearplants))["pr"]-mean(MDMp1$diffMaker%*%nuclearplants$cost), 10)==0)

tm1 <- pairmatch( mdist(pr ~ t1 + t2, data =nuclearplants), controls=2 )
MDMt1 <- MatchedDiffMaker(tm1, type="dense")
stopifnot(round(coef(lm(cost~pr+tm1, data=nuclearplants))["pr"]-mean(MDMt1$diffMaker%*%nuclearplants$cost), 10)==0)

fm1 <- fullmatch(mdist(pr~t1+t2, data=nuclearplants))
MDMf1 <- MatchedDiffMaker(fm1, type="dense")

stopifnot(round(coef(lm(cost~pr+fm1, data=nuclearplants))["pr"]-
with(MDMf1, sum(diffMaker%*%nuclearplants$cost /
     (1+rep(theOthers.HowMany, theOthers.HowMany)) # size of matched set corresp. to matched diff
                )/
     sum(theOthers.HowMany/(theOthers.HowMany+1)) # Half the sum of harmonic means of (\#Tx, \#Ctl)
     ),
                10)==0
          )

aGlm <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
fm2 <- fullmatch(caliper(.2, aGlm ) )
stratumStructure(fm2)
MDMf2 <- MatchedDiffMaker(fm2, type='dense')

stopifnot(require("SparseM"))
sMDMp1 <- MatchedDiffMaker(pm1, type="SparseM-csr")
sMDMt1 <- MatchedDiffMaker(tm1, type="SparseM-csr")
sMDMf1 <- MatchedDiffMaker(fm1, type="SparseM-csr")
sMDMf2 <- MatchedDiffMaker(fm2, type="SparseM-csr")

stopifnot(with(nuclearplants,
               all.equal(MDMp1$diffMaker%*%cost,
                         sMDMp1$diffMaker%*%cost)),
          with(nuclearplants,
               all.equal(MDMt1$diffMaker%*%cost,
                         sMDMt1$diffMaker%*%cost)),
          with(nuclearplants,
               all.equal(MDMf1$diffMaker%*%cost,
                         sMDMf1$diffMaker%*%cost)),
          with(nuclearplants,
               all.equal(MDMf2$diffMaker%*%cost,
                         sMDMf2$diffMaker%*%cost))
          )

MatchedDiffTable(nuclearplants$cost, MDMf2)
stopifnot(all.equal(with(nuclearplants, MatchedDiffTable(cost, MDMf2)), MatchedDiffTable(nuclearplants['cost'], MDMf2)),
          all.equal(with(nuclearplants, MatchedDiffTable(cost, MDMf2)), subset(MatchedDiffTable(nuclearplants[c('cost', 'cap')], MDMf2), select=-cap)),
          all.equal(MatchedDiffTable(nuclearplants['cost'], MDMf2), MatchedDiffTable(nuclearplants['cost'], sMDMf2))
          )

apsd <- PSdiffs(aGlm, fm2)
stopifnot(all(abs(predict(aGlm)[as.character(apsd$theTx)] - predict(aGlm)[as.character(apsd$theCtl)] - apsd$diff)<.Machine$double.eps^.5))


apsd.c <- PSdiffs(aGlm, fm2, covariance.extractor=SandwichVcov)
all.equal(apsd[c("theTx","theCtl","diff")], apsd.c[c("theTx","theCtl","diff")])
