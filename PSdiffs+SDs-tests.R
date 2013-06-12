require(optmatch)
data(nuclearplants)
source("matchedDifferencesFromOptmatch.R")
source("PSdiffs+SDs.R")

### Cross-checking getPSdiffs (assisted by optmatch:matched.distances) against PSdiffs
aPS <- glm(pr~.-(pr+cost), family=binomial(), data=nuclearplants)
fm2 <- fullmatch(caliper(.2, aPS ) )
fm2psdiffs <- PSdiffs(aPS, fm2)
fm2psdiffs$diff <- abs(fm2psdiffs$diff)

PSdiffMats <- getPSdiffs(aPS)

fm2psdiffs.c1 <- 
do.call(rbind, 
        lapply(matched.distances(fm2, PSdiffMats$differences, preserve.unit.names = TRUE),
               function(mat) {
                 n <- length(mat)
                 if (nrow(mat)>1)
                   {
                     theTx <- dimnames(mat)[[1]]
                     theCtl <- rep(dimnames(mat)[[2]], n)
                   } else
                 {
                   theTx <- rep(dimnames(mat)[[1]], n)           
                   theCtl <- dimnames(mat)[[2]]
                 }
                 data.frame(theTx, theCtl, diff=abs(as.vector(mat)))
               }
               )
        )
stopifnot(isTRUE(all.equal(as.character(fm2psdiffs[[1]]), as.character(fm2psdiffs.c1[[1]]))),
          isTRUE(all.equal(as.character(fm2psdiffs[[2]]), as.character(fm2psdiffs.c1[[2]]))),
          isTRUE(all.equal(fm2psdiffs[[3]], fm2psdiffs.c1[[3]]))
          )

fm2psdiffs.c2 <- 
do.call(rbind, 
        lapply(matched.distances(fm2, PSdiffMats$SEs, preserve.unit.names = TRUE),
               function(mat) {
                 n <- length(mat)
                 if (nrow(mat)>1)
                   {
                     theTx <- dimnames(mat)[[1]]
                     theCtl <- rep(dimnames(mat)[[2]], n)
                   } else
                 {
                   theTx <- rep(dimnames(mat)[[1]], n)           
                   theCtl <- dimnames(mat)[[2]]
                 }
                 data.frame(theTx, theCtl, diff=abs(as.vector(mat)))
               }
               )
        )

stopifnot(isTRUE(all.equal(as.character(fm2psdiffs[[1]]), as.character(fm2psdiffs.c2[[1]]))),
          isTRUE(all.equal(as.character(fm2psdiffs[[2]]), as.character(fm2psdiffs.c2[[2]]))),
          isTRUE(all.equal(fm2psdiffs[[4]], fm2psdiffs.c2[[3]]))
          )
### End cross-check of getPSdiffs against PSdiffs

## Now cross-checking PSSEdist against getPSdiffs
names(dimnames(PSdiffMats$SEs)) <- NULL
stopifnot(isTRUE(all.equal(PSdiffMats$SEs, PSSEdist(aPS)[[1]])))
## End cross-check of PSSEdist against getPSdiffs
