cumulantfct <- function(t, vars, thematch, pscore.linear,
                        derivs=0 # NB: poss. to pass expressions in lieu of derivs,
###                                using MS.cumfct.expr.factory(),
###                                although not presently (Rev116) tested.
###                                Multiple numeric derivs not yet working but see tests (R117).
                        )
{
  stopifnot(inherits(thematch, "optmatch"), is.data.frame(vars) || is.numeric(vars), 
            is.null(pscore.linear) || is.numeric(pscore.linear), 
            if (is.data.frame(vars)) nrow(vars)==length(thematch) else length(vars)==length(thematch),
            length(thematch)=={if (is.data.frame(vars)) nrow(vars) else length(vars)},
            all(!is.na(vars)),
            is.null(pscore.linear) || all(!is.na(pscore.linear)),
            is.null(pscore.linear) || length(pscore.linear)==length(thematch)
            )
	if (!is.null(pscore.linear) && all(pscore.linear>=0)) warning("pscore.linear is all positive. Are these really linear scores, ie logits of conditional probabilities?")
	if (is.numeric(vars)) vars <- data.frame(v=vars)
	if (any(c("TrTx","pSlIn") %in% names(vars))) stop("please use variable names other than 'TrTx' or 'pSlIn'.")
	if (is.null(pscore.linear)) pscore.linear <- rep(0, length(thematch))

  mtd <- matched(thematch)
  if (!all(mtd)) {
    thematch <- thematch[mtd, drop=TRUE]
  }
  thescale <- sum(RItools:::harmonic(data.frame(Tx.grp=attr(thematch, "contrast.group"),
                                                stratum.code=thematch)
                                     ) # H-means of tx & ctl group sizes, by matched set
                  )/2 # half the sum of H-means
  thedfs <- split(data.frame(TrTx=attr(thematch, "contrast.group"),
                             pSlIn=pscore.linear[mtd], vars[mtd,]), 
                  thematch)
  all.ms <- sapply(thedfs, nrow)

  wking.MS.cumfct.exprs <- if (is.numeric(derivs))
    {
      wking.ms <- c(2:8, seq(10, ceiling(max(all.ms)/10)*10, by=10))
      reduced.ms <- wking.ms[1 + findInterval(all.ms, wking.ms + 1)]
      MS.cumfct.expr.factory(reduced.ms, derivs)
    } else {
      stopifnot(is.list(derivs), all(substring(names(derivs),1,5)=="m.is."),
                max(as.numeric(substring(names(derivs),6))) >= max(all.ms) )
      derivs
    }
  
  mscums <- lapply(thedfs, MScumfct, t=t, derivs=wking.MS.cumfct.exprs)

  ans <- mscums[[1]]
  if (length(mscums)-1) for (ii in 2:length(mscums)) ans <- ans + mscums[[ii]]
  
  ans <- ans/(thescale^attr(wking.MS.cumfct.exprs, "which.derivatives"))
  cbind(t=t, ans)
	}
	
MScumfct <- function(adf, t, derivs=0)
{
  stopifnot(is.data.frame(adf), length(adf)>=3,
            all(c("TrTx","pSlIn") %in% names(adf)),
            is.logical(adf$TrTx), is.numeric(t), is.numeric(derivs) || is.list(derivs),
            !is.list(derivs) || (!is.null(names(derivs)) &&
                                 all(!is.na(availsizes <- as.numeric(substring(names(derivs),6)) )))
            )
ifx <- is.finite(adf[["pSlIn"]])
cnms <- names(adf)!=c("TrTx")
adf[ifx,cnms] <- scale(adf[ifx,cnms], scale=F, center=T)
#ix <- !(names(adf)%in%c("TrTx","pSlIn"))
#adf[ix] <- sapply(adf[ix], function(x) scale(x, scale=F, center=T) )
if (is.numeric(derivs))
  MSCumfct(adf, t, derivs=as.integer(derivs)) else {
    stopifnot(is.finite(m <- min(availsizes[availsizes>=nrow(adf)]) ))
    evalMScumfct(pad.bottom.rows(adf=adf, padto=m),t,cfexpr=derivs[[paste("m.is.",as.integer(m),sep="")]])
  }
  }

MS.cumfct.expr.factory <- function(m, derivs)
  {
DD <- function(expr,name, order = 1) {
  if(order == 1) D(expr,name)
  else DD(D(expr, name), name, order - 1)
}
make.MS.cumfct.expr <- function(m)
  {
    vnms <- paste("v",1:m, sep="")
    psnms <- paste("ps",1:m, sep="")
    thesummands <- paste("exp(t*", vnms, " + ", psnms, ")", sep="")
    exptext <- paste("log(", paste(thesummands, collapse=" + "), ")", sep=" ")
    expr <- parse(text=exptext, srcfile=NULL)
    if (derivs==0) return(expr)
    
    if (length(derivs)==1) return(DD(expr, "t", order=derivs)) else
    {
      expr <- lapply(derivs, function(d) DD(expr, "t", order=d))
      exptexts <- lapply(expr, function(x) deparse(x, control=c("keepInteger", "keepNA")) )
      exptexts[1:(length(exptexts)-1)] <-
        lapply(exptexts[1:(length(exptexts)-1)], function(charvec) {
          lcv <- length(charvec)
          charvec[lcv] <- paste(charvec[lcv], ",")
          charvec
        } )
      Exptext <- c("c(" , unlist(exptexts), ")")
      return(parse(text=Exptext, srcfile=NULL))
    }
}
m <- setdiff(m, 1)
ans <- c(list(rep(parse(text="rep(0, length(t))",srcfile=NULL),
                  length(derivs))), #If m is 1 then cumfct is automatically 0 (& its derivs vanish)
         lapply(m, make.MS.cumfct.expr))
names(ans) <- paste("m.is.", c(1L, as.integer(m)), sep="")
attr(ans, "which.derivatives") <- derivs 
ans
}

makeVals <- function(vnm, adf, t) { # Can this be done more directly using substitute()?
  m <- nrow(adf)
  thevals <- unlist(adf[c(vnm, "pSlIn")])
  names(thevals) <- c(paste("v",1:m, sep=""), paste("ps",1:m, sep=""))
  c(list(t=t), as.list(thevals))
}

MSCumfct <- function(adf, t, derivs)
  {
m <- nrow(adf)
exprlst <- MS.cumfct.expr.factory(m=m, derivs=derivs)

trtx <- which(names(adf)=="TrTx")
if (sum(adf$TrTx)>=2) adf[-trtx] <- sapply(adf[-trtx], function(x) -1*x)

thevars <- setdiff(names(adf), c("TrTx","pSlIn"))
sapply(thevars, function(nm, adf, t) eval(exprlst[[paste("m.is.",as.integer(m),sep="")]], makeVals(nm, adf, t)), adf=adf, t=t)

}

evalMScumfct <- function(adf, t, cfexpr)
  {
### would be nice to have a check that cfexpr has as many free variables as adf has rows.
### (the check that there are no more free variables shoule be automatic, as otherwise this will fail.)

trtx <- which(names(adf)=="TrTx")
if (sum(adf$TrTx)>=2) adf[-trtx] <- sapply(adf[-trtx], function(x) -1*x)

thevars <- setdiff(names(adf), c("TrTx","pSlIn"))
sapply(thevars, function(nm, adf, t) eval(cfexpr, makeVals(nm, adf=adf, t=t)), adf=adf, t=t)
  }
pad.bottom.rows <- function(adf, padto)
{
stopifnot(is.data.frame(adf), all.equal(c("TrTx","pSlIn"), names(adf)[1:2]), is.numeric(padto))
nvars <- length(adf)-2
xtrarows <- max(0,padto -nrow(adf))
if (!xtrarows) return(adf)

newdf <- data.frame(TrTx=rep(FALSE, xtrarows), pSlIn=rep(-Inf, xtrarows))
if (nvars) 
	{
	varnms <- setdiff(names(adf), c("TrTx","pSlIn"))
	newdf <- data.frame(newdf, matrix(0,xtrarows,nvars))
	names(newdf) <- c("TrTx","pSlIn", varnms)
}
rbind(adf, newdf)
}
VSK <- function(pscore, thematch, vars)
  {
    if (is.null(pscore)||is.numeric(pscore)) pscore.linear <- pscore
    if (inherits(pscore, "glm")) pscore.linear <- pscore$linear.predictors

K1=cumulantfct(t=0,vars=vars,thematch=thematch,pscore.linear=pscore.linear, derivs=1)[,-1]
K2=cumulantfct(t=0,vars=vars,thematch=thematch,pscore.linear=pscore.linear, derivs=2)[,-1]
K3=cumulantfct(t=0,vars=vars,thematch=thematch,pscore.linear=pscore.linear, derivs=3)[,-1]
K4=cumulantfct(t=0,vars=vars,thematch=thematch,pscore.linear=pscore.linear, derivs=4)[,-1]

    if (is.null(dim(vars)))
      {return(c(Exp=K1,Var=K2,Skew=K3/((K2)^(3/2)),Kurt=K4/(K2^2)))
     } else return(data.frame(Exp=K1, Var=K2, Skew=K3/((K2)^(3/2)),
                              Kurt=K4/(K2^2), row.names=names(vars)))
}

Fhat <- function(x, themean, thevar, theK3=0, theK4=0)
 {
y <- (x-themean)/sqrt(thevar)
firsto <- pnorm(y)
secondo <- dnorm(y)*(theK3/6)*(y^2-1)
thirdo <-
  dnorm(y)*( (theK4/24)*(y^3 - 3*y) +
            (theK3^2/72)*(y^5 - 10*y^3 + 15*y)
            )
ans <- data.frame(x=x, first=firsto, second=(firsto-secondo), third=(firsto-secondo-thirdo))
if (length(x)==1) {ans <- drop(ans)} else row.names(ans) <- names(x)
ans
 }
coverage <- function(pscore, thematch, vars, alpha=.05)
  {
Ks.notilt <- VSK(NULL, thematch, vars)
Ks.tilt <- VSK(pscore, thematch, vars)
upper <- with(as.list(Ks.notilt), Exp + qnorm(1-alpha/2)*sqrt(Var))
lower <- with(as.list(Ks.notilt), Exp + qnorm(alpha/2)*sqrt(Var))

cov.notilt <- with(as.list(Ks.notilt),
                   Fhat(upper, Exp, Var, Skew, Kurt)[-1] - Fhat(lower, Exp, Var, Skew, Kurt)[-1])

cov.tilt <- with(as.list(Ks.tilt),
                 Fhat(upper, Exp, Var, Skew, Kurt)[-1] - Fhat(lower, Exp, Var, Skew, Kurt)[-1])

rbind(P=c(Ks.notilt,cov.notilt), Q=c(Ks.tilt,cov.tilt))
}
contig.EQ <- function(pscore, thematch, vars)
  {
pscore <- groupCenter(pscore, thematch)
mtd <- matched(thematch)
theweights <- unsplit(lapply(table(thematch[mtd,drop=TRUE]), function(m) rep(1/m,m)),
                      thematch[mtd, drop=TRUE])
if (is.data.frame(pscore) && length(pscore)>1) stop("contig.EQ accepts 1 pscore at a time.")
if (is.data.frame(pscore)) pscore <- pscore[[1]]
vars <- vars[mtd,]
thescale <- sum(RItools:::harmonic(data.frame(Tx.grp=attr(thematch, "contrast.group"),
                                              stratum.code=thematch)
                                   ) # H-means of tx & ctl group sizes, by matched set
                )/2 # half the sum of H-means
sapply(theweights*pscore*vars, sum)/thescale
  }

pvals.PandQ <- function(pscore, thematch, vars, side="two", samples=0)
  {
    require("RItools")
    if (side!="two") stop("Only 2-sided p-values implemented at present.")
	stopifnot(is.numeric(samples) || is.logical(samples) || is.data.frame(samples) || is.null(samples),
	!is.data.frame(samples) || nrow(samples)==length(thematch))
	
	if (is.logical(samples)&&samples) samples <- 1
	
    Ks.notilt <- VSK(NULL, thematch, vars)
    Ks.tilt <- VSK(pscore, thematch, vars)
    newdf <- data.frame(TrTx=attr(thematch, "contrast.group"), vars)
    thediffs <- xBalance(fmla=formula(newdf), strata=thematch,
                         data=newdf, report="adj.mean.diffs")[[1]][,'adj.diff',,drop=T]
    thediffs <- abs(thediffs)
p.notilt <- with(as.list(Ks.notilt),
                 Fhat(-thediffs, Exp, Var, Skew, Kurt)[-1] + 1 - Fhat(thediffs, Exp, Var, Skew, Kurt)[-1])

p.tilt <- with(as.list(Ks.tilt),
               Fhat(-thediffs, Exp, Var, Skew, Kurt)[-1] + 1 - Fhat(thediffs, Exp, Var, Skew, Kurt)[-1])
p.sd <- sqrt(Ks.notilt['Var'])
ans <- list(P=cbind(`~Exp/SD.P`=0,
              Ks.notilt[,"Exp", drop=TRUE]/p.sd,
              Ks.notilt[,2:4], p.notilt),
            Qhat=cbind(contig.EQ(pscore, thematch, vars)/p.sd,
              Ks.tilt[,"Exp", drop=TRUE]/p.sd,
              Ks.tilt[,2:4], p.tilt))
    names(ans[["P"]])[1:2] <- names(ans[["Qhat"]])[1:2]<- c('~Exp/SD.P', 'Exp/SD.P')
    if (is.numeric(samples) && samples && inherits(pscore, "glm"))
      {
      morepscores <- samplePropensities(pscore, samples=samples)
      names(morepscores) <- paste("Qdraw", 1:samples, sep="")
	  }
	  
	if (is.data.frame(samples)) morepscores <- samples
	if (exists('morepscores'))
	  {
      ans2 <- lapply(morepscores, function(theps) {
           theKs <- VSK(theps, thematch, vars)
           thep <- with(as.list(theKs),
               Fhat(-thediffs, Exp, Var, Skew, Kurt)[-1] + 1 - Fhat(thediffs, Exp, Var, Skew, Kurt)[-1])
           A <- cbind(contig.EQ(theps, thematch, vars)/p.sd,
                      theKs[,"Exp", drop=TRUE]/p.sd,
                      theKs[,2:4], thep)
           names(A)[1:2] <-  c('~Exp/SD.P', 'Exp/SD.P')
           A
      }
                     )
      return(c(ans, ans2))
    } else return(ans)
}

pvals_cumbased <- function(pscore, thematch, vars, side="two", samples=0)
  {
    require("RItools")
    if (side!="two") stop("Only 2-sided p-values implemented at present.")
	stopifnot(is.numeric(samples) || is.logical(samples) || is.data.frame(samples) || is.matrix(samples) || is.null(samples),
                  !is.data.frame(samples) || nrow(samples)==length(thematch),
                  !is.matrix(samples) || nrow(samples)==sum(!is.na(coef(pscore)))
                  )
	
	if (is.logical(samples)&&samples) samples <- 1
	if (is.matrix(samples)) samples <- reconstructPropensities(samples, pscore, strata=thematch)
    Ks.notilt <- VSK(NULL, thematch, vars)
    Ks.tilt <- VSK(pscore, thematch, vars)
    newdf <- data.frame(TrTx=attr(thematch, "contrast.group"), vars)
    thediffs <- xBalance(fmla=formula(newdf), strata=thematch,
                         data=newdf, report="adj.mean.diffs")[[1]][,'adj.diff',,drop=T]
    thediffs <- abs(thediffs)
p.notilt <- with(as.list(Ks.notilt),
                 Fhat(-thediffs, Exp, Var, Skew, Kurt)[c("first", "third")] + 1 - Fhat(thediffs, Exp, Var, Skew, Kurt)[c("first", "third")])

p.tilt <- with(as.list(Ks.tilt),
               Fhat(-thediffs, Exp, Var, Skew, Kurt)[c("first", "third")] + 1 - Fhat(thediffs, Exp, Var, Skew, Kurt)[c("first", "third")])
p.sd <- sqrt(Ks.notilt['Var'])
ans <- list(P=cbind(#`~Exp/SD.P`=0,
              Ks.notilt[,"Exp", drop=TRUE]/p.sd,
              rep(1, length(p.sd)),
              Ks.notilt[,3:4], p.notilt),
            Qhat=cbind(#contig.EQ(pscore, thematch, vars)/p.sd,
              Ks.tilt[,"Exp", drop=TRUE]/p.sd,
              sqrt(Ks.tilt[,"Var", drop=TRUE])/p.sd,
              Ks.tilt[,3:4], p.tilt))
    names(ans[["P"]])[1:2] <- names(ans[["Qhat"]])[1:2]<- c('Exp/SD.P', 'SD/SD.P')
    if (is.numeric(samples) && samples && inherits(pscore, "glm"))
      {
      morepscores <- samplePropensities(pscore, samples=samples)
      names(morepscores) <- paste("Qdraw", 1:samples, sep="")
	  }
	  
	if (is.data.frame(samples)) morepscores <- samples
	if (exists('morepscores'))
	  {
      ans2 <- lapply(morepscores, function(theps) {
           theKs <- VSK(theps, thematch, vars)
           thep <- with(as.list(theKs),
               Fhat(-thediffs, Exp, Var, Skew, Kurt)[c("first", "third")] + 1 - Fhat(thediffs, Exp, Var, Skew, Kurt)[c("first", "third")])
           A <- cbind(#contig.EQ(theps, thematch, vars)/p.sd,
                      theKs[,"Exp", drop=TRUE]/p.sd,
                      sqrt(theKs[,"Var",drop=TRUE])/p.sd,
                      theKs[,3:4], thep)
           names(A)[1:2] <-  c('Exp/SD.P', 'SD/SD.P')
           A
      }
                     )
      return(c(ans, ans2))
    } else return(ans)
}


## For checking purposes:
Eadjmeandiff <- function(vars, thematch, pscore.linear)
  {
  stopifnot(inherits(thematch, "optmatch"), is.data.frame(vars) || is.numeric(vars), 
            is.null(pscore.linear) || is.numeric(pscore.linear), 
            if (is.data.frame(vars)) nrow(vars)==length(thematch) else length(vars)==length(thematch),
            length(thematch)=={if (is.data.frame(vars)) nrow(vars) else length(vars)},
            all(!is.na(vars)),
            is.null(pscore.linear) || all(!is.na(pscore.linear)),
            is.null(pscore.linear) || length(pscore.linear)==length(thematch)
            )
	if (!is.null(pscore.linear) && all(pscore.linear>=0)) warning("pscore.linear is all positive. Are these really linear scores, ie logits of conditional probabilities?")
	if (is.numeric(vars)) vars <- data.frame(v=vars)
	if (any(c("TrTx","pSlIn") %in% names(vars))) stop("please use variable names other than 'TrTx' or 'pSlIn'.")
	if (is.null(pscore.linear)) pscore.linear <- rep(0, length(thematch))

  mtd <- matched(thematch)
  if (!all(mtd)) {
    thematch <- thematch[mtd, drop=TRUE]
  }
 thedfs <- split(data.frame(pSlIn=pscore.linear[mtd], vars[mtd,]), 
                  thematch)
 thedfs <- lapply(thedfs, function(df) as.data.frame(scale(df, center=TRUE, scale=FALSE)))
  
  manyones <- tapply(attr(thematch, "contrast.group"), thematch, function(x) sum(x)>1)

if (any(manyones)) thedfs[manyones] <- lapply(thedfs[manyones], function(df) -1*df)

msmeans <- lapply(thedfs, function(df)
         {
           exp.lps <- exp(df[["pSlIn"]])
           sapply(df[(names(df)!="pSlIn")], function(v) sum(v*exp.lps))/
             sum(exp.lps)
         }
                  )
  ans <- msmeans[[1]]
   if (length(msmeans)>1) for (ii in 2:length(msmeans)) ans <- ans + msmeans[[ii]]
  
    thescale <- sum(RItools:::harmonic(data.frame(Tx.grp=attr(thematch, "contrast.group"),
                                                stratum.code=thematch)
                                     ) # H-means of tx & ctl group sizes, by matched set
                  )/2 # half the sum of H-means

  ans/thescale
  }
