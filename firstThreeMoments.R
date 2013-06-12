### requires optmatch version 0.8 or later, due to matched() bug in versions 0.8_{1,2}
first3M <- function(x, thematch, ...) UseMethod("first3M")
first3M.formula <- function(x, data, subset,...)
  {
### Cribbing from match_on:
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("x", "data", "subset"), # maybe later add "na.action"
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    names(mf)[names(mf) == "x"] <- "formula"
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())  
    if (dim(mf)[2] < 2) {
    stop("Formula must have a right hand side with at least one variable.")
  }
    x <- subset(model.matrix(x, mf), T, -1) # drop the intercept
    first3M(x,...)
  }
first3M.default <- function(x, thematch, pscore.linear=NULL,...) # formula, ... methods reshape and send here
  {
  stopifnot(inherits(thematch, "optmatch"),
            is.null(pscore.linear) || is.numeric(pscore.linear) || is.data.frame(pscore.linear),
            !is.data.frame(pscore.linear) || all(sapply(pscore.linear, is.numeric)),
            nrow(data) == {if (is.data.frame(pscore.linear))
                             nrow(pscore.linear) else length(pscore.linear)},
            is.null(pscore.linear) || all(complete.cases(pscore.linear)),
            (is.data.frame(x) && all(sapply(x, is.numeric)) ) || is.numeric(x),
            require("SparseM")
            )
  if (!is.null(pscore.linear) && all(pscore.linear>=0))
    warning("pscore.linear is all positive. Are these really linear scores, ie logits of conditional probabilities?")
  if (is.null(pscore.linear)) pscore.linear <- rep(0, length(thematch))
  ismatched <- !is.na(thematch) 
  thematch <- thematch[ismatched, drop=TRUE]
  designmat <- SparseMMFromFactor(thematch)

### Group-mean-centering the xes and the pscores.
### Given titling, this isn't enough to center the moments,
### but it should numerically stabilize the tilting and it's
### first of 2 necessary steps for handling many-one matched sets
  
  pscore.adj <- resid(SparseM:::slm.fit(designmat, subset(pscore.linear, ismatched)))
  x.adj <- resid(SparseM:::slm.fit(designmat, as.matrix(subset(x, ismatched))))  

### Second of 2 necessary preparations for ignoring distinction between
### many-one and one-many matched sets. 
  cfactor <- ifelse(identify_manyone_matchedsets(thematch), -1, 1)
  pscore.adj <- pscore.adj * cfactor
  x.adj <- x.adj * cfactor

### Tilting
  ps.wts <- exp(pscore.adj)
  rm(pscore.adj)
  ms.sum.maker <- stratum_summing_SparseM(thematch)

### Harmonic overall divisor, ie half the sum of harmonic means of matched set sizes
  h <- sum( (1 + 1/(matchedsetsize(thematch)-1))^(-1) )
### Now define a function of a single vector of PS-based weights to compute the various moments.
  momentmaker <- function(pswts)
    {
        ms.wtsums <- ms.sum.maker %*% pswts
        ms.wtsums <- drop(ms.wtsums)
        m1 <- ms.sum.maker %*% (x.adj * pswts)
        m1 <- as.matrix(m1)
        m1 <- m1/ms.wtsums
        m2 <- ms.sum.maker %*% (x.adj^2 * pswts)
        m2 <- as.matrix(m2)
        m2 <- m2/ms.wtsums
        m3 <- ms.sum.maker %*% (x.adj^3 * pswts)
        m3 <- as.matrix(m3)
        m3 <- m3/ms.wtsums
        
        c2 <- m2 - m1^2
        c3 <- m3 -3*m1*m2 + 2*m1^3

        cSc2 <- colSums(c2)
        rbind(Exp=colSums(m1)/h,Var=cSc2/h^2,Skew=colSums(c3)/cSc2^1.5)
    }

  ps.wts <- as.data.frame(ps.wts)
  sapply(ps.wts, momentmaker, simplify = "array")
}

identify_manyone_matchedsets <- function(thematch, expanded=TRUE)
  {
    thetx <- attr(thematch, "contrast.group")
    manyone <- tapply(thetx, thematch, sum) > 1
    if (expanded) unsplit(manyone, thematch) else manyone
  }

### Not currently using the below, but here it is just in case:
matchedsetsize <- function(thematch, expanded=TRUE, logged=FALSE)
  {
stopifnot(inherits(thematch, "optmatch"),
          !is.null(attr(thematch, "contrast.group"))
          )
thetab <- table(thematch)
if (any(thetab<2)) stop("? Found matched sets of size 1 or less.  Please fix.")
if (logged) thetab <- log(thetab)
if (expanded) unsplit(thetab, thematch) else thetab
}
