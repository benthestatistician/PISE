makeXandCombdiffX <- function(fmla, strata, data, na.rm=FALSE) #, impfn=median
  {
stopifnot(inherits(fmla, "formula"),
          inherits(strata, "factor") || inherits(strata, "formula"),
          is.null(data) || is.data.frame(data)#,is.function(impfn)
          )


### NA Handling, copied from xBalance ##

  if (na.rm==TRUE)
    {
      tfmla <- terms.formula(ufmla,dat=data, keep.order=TRUE)
    } else
  {
    data <- RItools:::naImpute(fmla,data)#,impfn
    tfmla <- attr(data, 'terms')
  }
### End NA handling ###

mf <- model.frame(tfmla, data)
###Extract the treatment var
  if (!attr(tfmla, "response")>0)
    stop("fmla must specify a treatment group variable")

  zz <- eval(tfmla[[2]], mf, parent.frame())
  if (!is.numeric(zz) & !is.logical(zz))
    stop("LHS of fmla should be logical or numeric")
  if (any(is.na(zz))) stop('NAs on LHS of fmla not allowed.')
### End extract treatment var


mm1 <- model.matrix(tfmla, mf)
vnms <- colnames(mm1)
vnms <- vnms[vnms!="(Intercept)"]
mm1 <- mm1[,vnms]

if (inherits(strata, "formula"))
  {
  ss <- eval(attr(terms(strata), "variables"), data)
    cenv <- environment()
    parent.env(cenv) <- environment(strata)
    environment(strata) <- cenv
  thelm1 <- lm(update(strata, mm1~.), data=data)
  } else thelm1 <- lm(mm1 ~ strata)

xmat <- resid(thelm1)
if (ONEVAR <- is.null(dim(xmat))) xmat <- matrix(xmat, length(xmat),1,dimnames=list(names(xmat), vnms))
###zz <- zz[-thelm1$na.action] not necessary after all
thelm2 <- update(thelm1, ~.+zz)
combdiff <- if (ONEVAR) {coef(thelm2)['zz']} else coef(thelm2)['zz',]
names(combdiff) <- vnms

if (inherits(strata, "formula")) {
  strata <- if (length(ss)-1) interaction(ss, drop=TRUE) else as.factor(ss[[1]])
  }

list(xmat=xmat, combdiff=combdiff, strata=strata)
}

oneOverStratSize <- function(strata)
  {ts <- table(strata)
   mtd <- if (inherits(strata, "optmatch")) matched(strata) else !is.na(strata)
   strata <- strata[mtd,drop=TRUE]
    ans <- unsplit(1/ts, strata, drop=FALSE)
 names(ans) <- names(strata)
 ans
 }

getVarCombdiffX <- function(xmat, thematch=NULL, tol=sqrt(.Machine$double.eps))
{
stopifnot(is.matrix(xmat) || (is.list(xmat) && is.matrix(xmat[["xmat"]])),
          (is.list(xmat) && "strata"%in%names(xmat)) || !is.null(thematch),
          (is.null(thematch) && inherits(xmat$strata, "optmatch")) ||
          (!is.null(thematch) && inherits(thematch, "optmatch")),
          sum(matched(if (is.null(thematch)) xmat$strata else thematch))==nrow(xmat),
          is.list(xmat) || round(sum(xmat[,1]),log10(.Machine$double.eps)/2)==0)

if (is.null(thematch)) thematch <- xmat$strata 
if (is.list(xmat)) xmat <- xmat[[1]]

getN <- function(thematch)
{
stopifnot(inherits(thematch, "optmatch"))
nlevels(thematch[matched(thematch), drop=T])
}
numstrat <- getN(thematch)

dvec <- oneOverStratSize(thematch)
K <- ncol(xmat)


ssqx0 <- apply(xmat*xmat, 2, sum)/numstrat
sdxreciprocals <- ssqx0^(-.5)
x.times.sxm1 <- xmat*matrix(sdxreciprocals, nrow(xmat), K, byrow=TRUE)
invertme <- crossprod(x.times.sxm1, dvec*x.times.sxm1)
Xsvd <- svd(invertme)
Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
inverted <- if (all(Positive)) 
  Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
else if (!any(Positive)) 
  array(0, dim(X)[2L:1L])
else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))

list(V=sdxreciprocals*inverted*matrix(sdxreciprocals, K, K, byrow=TRUE),
     df=sum(Positive) )
}

###getXbetahat <- function(X,...) UseMethod("getXbetahat")
getXbetahat <- function(X, varCombdiffX)
  {
    stopifnot(is.list(X),
              all.equal(names(X), c("xmat", "combdiff", "strata")),
              is.list(varCombdiffX),
              all.equal(names(varCombdiffX),c("V", "df")))

    getNH.over.2 <- function(thematch)
{
  stopifnot(inherits(thematch, "optmatch"))
  ms <- table(thematch[matched(thematch), drop=T])
  sum((ms-1)/ms)
}

    nHover2 <- getNH.over.2(X$strata)
    betahat <- nHover2*varCombdiffX$V %*%X$combdiff
    X$xmat %*% betahat
  }
ReplaceXdeviationsWithXdiffs <- function(X)
  {
   stopifnot(is.list(X),
              all.equal(names(X), c("xmat", "combdiff", "strata"))
             )
mdm <- MatchedDiffMaker(X$strata[matched(X$strata), drop=TRUE])$diffMaker # MatchedDiffMaker() lives in matchedDifferencesFromOptmatch.R
mat <- mdm %*% X$xmat
mat <- as.matrix(mat)
colnames(mat) <- colnames(X$xmat)
X$xmat <- mat
X
}
###getXbetahatSDs <- function(X,...) UseMethod("getXbetahatSDs")
getXbetahatSDs <- function(X, varCombdiffX)
  {
    stopifnot(is.list(X),
              all.equal(names(X), c("xmat", "combdiff", "strata")),
              is.list(varCombdiffX),
              all.equal(names(varCombdiffX),c("V", "df")))
   
    sqrt( rowSums( (X$xmat %*% varCombdiffX$V)*X$xmat ) )
  }

imbalChisq <- function(X,...) UseMethod("imbalChisq")
imbalChisq.default <- function(X, varCombdiffX) # option to get this from a clogit fit?
  {
stopifnot(is.list(X),
          all.equal(names(X), c("xmat", "combdiff", "strata")),
          is.list(varCombdiffX),
          all.equal(names(varCombdiffX),c("V", "df")))

xbh <- getXbetahat(X, varCombdiffX)
dvec <- oneOverStratSize(X$strata)
chis <- crossprod(xbh*dvec, xbh)
df <- varCombdiffX$df
c(chisquare=chis, df=df, p.value=pchisq(chis, df, lower.tail=FALSE))
}
imbalChisq.clogit <- function(X,tol=.Machine$double.eps,...)
  { 
covmat <- survival:::vcov.coxph(X)
coefs <- coef(X)
stopifnot(setequal(names(coefs), colnames(covmat)))
if (!all.equal(names(coefs), colnames(covmat))) coefs <- coefs[colnames(covmat)]
keepme <- !is.na(coefs)
coefs <- coefs[keepme]
covmat <- covmat[keepme, keepme, drop=FALSE]
### Invert covariance:
Xsvd <- svd(covmat)
Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
inverted <- if (all(Positive)) 
  Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
else if (!any(Positive)) 
  array(0, dim(X)[2L:1L])
else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
### end inversion

chis <- crossprod(coefs, inverted%*%coefs)
df <- sum(Positive)
c(chisquare=chis, df=df, p.value=pchisq(chis, df, lower.tail=FALSE))
}

extrapZ <- function(X,...) UseMethod("extrapZ")
extrapZ.default <- function(X, varCombdiffX=NULL, verbose=FALSE, useMatchedDifferencesNotMatchedDeviations=FALSE,...) 
  {
stopifnot(is.list(X),
          setequal(names(X), c("xmat", "combdiff", "strata")) ||
          setequal(names(X), c("xbetahat", "xbetahat.sds")),
          is.list(varCombdiffX) || (is.null(varCombdiffX) && names(X)[1]=="xbetahat"),
          is.null(varCombdiffX) || all.equal(names(varCombdiffX),c("V", "df")))

if (names(X)[1]=="xbetahat")
   {
     xbh <- X$xbetahat
     xbh.sds <- X$xbetahat.sds
   } else {
     if (useMatchedDifferencesNotMatchedDeviations) X <- ReplaceXdeviationsWithXdiffs(X)
xbh <- as.vector(getXbetahat(X, varCombdiffX))
xbh.sds <- getXbetahatSDs(X, varCombdiffX)
}
mx.xbh <- max(abs(xbh))
mx.xbh.sds <- max(xbh.sds)
z.newRatio <- max(abs(xbh)/pmax(median(xbh.sds), xbh.sds))
z.star <- mx.xbh/mx.xbh.sds

nc <- length(xbh)
ans <- list(maxPSdiff=mx.xbh,
            n.comps=nc,
            dev.to.maxSD=z.star,
            p=min(1,pnorm(-z.star)*nc*2),
            dev.to.medianSD.or.SD=z.newRatio,
            `p `=min(1,pnorm(-z.newRatio)*nc)
            )

if (verbose)
  {
    ans$PShat <- summary(xbh)
    ans$PShatSDs <- summary(xbh.sds)
    return(ans)
  } else return(unlist(ans))
}
extrapZ.clogit <- function(X,thematch,useMatchedDifferencesNotMatchedDeviations=FALSE,intermediateResults=FALSE,...)
  {
xm <- model.matrix(X)
if (!all(row.names(xm) %in% names(thematch))) stop("row names of 'X', given fitted clogit, don't match names of 'thematch'.")
them <- thematch[match(row.names(xm), names(thematch)), drop=TRUE]

newx <-  if (useMatchedDifferencesNotMatchedDeviations)
  {
    DM <- MatchedDiffMaker(them)
    mat <- as.matrix(DM$diffMaker %*% xm)
    colnames(mat) <- colnames(xm)
    rownames(mat) <- DM$rownames
    mat
  } else resid(lm(xm ~ them))

  coef <- ifelse(is.na(X$coefficients), 0, X$coefficients)
  newx <- newx[,names(coef)]
  xbh <- newx %*% coef
  xbh.sds <- sqrt(rowSums((newx %*% X$var) * newx))

ans <- list(xbetahat=xbh,xbetahat.sds=xbh.sds)
if (intermediateResults) ans else extrapZ(ans,...)
}
