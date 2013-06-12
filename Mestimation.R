if (!("GetInvObservedInfo"%in%ls())) stop("Have to source 'sandwich_glm.R' first.")
IDmanyonematchedsets <- function(omatch)
  {
    manyone <- unsplit(tapply(attr(omatch, "contrast.group"), omatch,
                              function(x) rep(as.logical(pmax(0,sum(x)-1)), length(x)) ),
                       omatch, drop=F)
    manyone[is.na(manyone)] <- FALSE
    manyone
  }
getpisFromPSfit <- function(object, ...) UseMethod("getpisFromPSfit")
getpisFromPSfit.logical <- function(object, strata=NULL,...)
  {
stopifnot(!any(is.na(object)), is.null(strata) || length(strata)==length(object),
          is.null(strata) || inherits(strata, "optmatch"))
if (inherits(strata, "optmatch") && any(object!=attr(strata, "contrast.group")))
  warning("specified treatment assignments differ from those implicit in matching; using specified.")
if (is.null(strata))
  {
pis <- rep(mean(object), length(object))
  } else {
pis <- as.numeric(object)
mtd <- matched(strata)
object <- object[mtd]
strata <- strata[mtd]
pis[mtd] <- unsplit(tapply(object, strata, function(x) rep(mean(x), length(x))),
                    strata, drop=F)
}

names(pis) <- names(object)
pis
}
getpisFromPSfit.glm <- function(object,...)
  {
    stopifnot(object$family$family=="binomial")
getpisFromPSfit(object$linear.predictor,...)
  }
getpisFromPSfit.numeric <- function(object, strata=NULL, ...)
  {
stopifnot(!any(is.na(object)), !any(is.infinite(object)),
          is.null(strata) || inherits(strata, "optmatch"),
          is.null(strata) || length(object)==length(strata),
          is.null(names(object)) || is.null(strata) || is.null(names(strata)) || all.equal(names(object), names(strata))
          )

if (is.null(strata)) { return(plogis(object)) } else
{
manyone <- IDmanyonematchedsets(strata)
mfd <- optmatch::matchfailed(strata)
mtd <- !is.na(strata) & !mfd
manyone <- manyone[mtd]
pis <- as.numeric(attr(strata, "contrast.group"))
names(pis) <- names(object)

strata <- factor(strata[mtd])
centeredthetas <- object[mtd]
centeredthetas <- SparseM::residuals.slm(SparseM::slm(centeredthetas~strata))
centeredthetas <- centeredthetas*(1-2*manyone)
numers <- exp(centeredthetas)
denoms <- tapply(numers, strata, function(x) rep(sum(x), length(x)))
denoms <- unsplit(denoms, strata, drop=F)
num.over.denom <- numers/denoms

pis[mtd] <- ifelse(manyone, 1-num.over.denom, num.over.denom)
pis
    }
}

covZtx <- function(object, ...) UseMethod("covZtx")
covZtx.default <- function(object, fittedpsmodel, outcome, type="IPW",...)
{
  if (!is.null(object)) stop("covZtx only implemented for full matched or unstratified (object==NULL) designs.")
  if (type!="IPW") stop("Only inverse probability weighting (IPW) for mu_c implemented at present.")
  pis <- getpisFromPSfit(fittedpsmodel)
  stopifnot(!any(is.na(pis)),
          is.null(fittedpsmodel) || inherits(fittedpsmodel, "glm"),
          is.data.frame(outcome) || is.numeric(outcome),
          nrow(as.data.frame(outcome))==length(pis),
          !any(is.na(outcome))
            )

  contribs <- pis<1
  pis <- pis[contribs]

  fitwts <- pis*(1-pis)

  outcome <- as.matrix(as.data.frame(outcome)[contribs,])
  outcome.weights <- 1/(1-pis)

  wted.mean.outcomes <- apply(pis*outcome, 2, sum)/sum(fitwts)
  outcome.centeredandweighted <- sweep(outcome*outcome.weights, 2, wted.mean.outcomes, check.margin=FALSE)

  if (is.null(fittedpsmodel)) {thexes <- thexes.centeredandweighted <- matrix(numeric(0),sum(contribs),0)} else
  {
  thexes <- model.matrix(fittedpsmodel)[contribs,]
  wted.mean.xes <- apply(thexes*fitwts, 2, sum)/sum(fitwts)
  thexes.centeredandweighted <- fitwts*sweep(thexes, 2, wted.mean.xes, check.margin=FALSE)
}
  list(VZtx=crossprod(thexes.centeredandweighted, thexes),
     CZtxZty=crossprod(thexes.centeredandweighted, outcome.centeredandweighted),
     diag.VZty=colSums(pis*outcome*outcome.centeredandweighted)
     )

}
  covZtx.optmatch <- function(omatch, fittedpsmodel, outcome,
                            pis=getpisFromPSfit(fittedpsmodel,strata=omatch), type="IPW",...)
{ # figure pis from omatch and fittedpsmodel?
stopifnot(is.numeric(pis), !any(is.na(pis)),
          inherits(omatch, "optmatch"),
          is.null(fittedpsmodel) || inherits(fittedpsmodel, "glm"),
          length(pis)==length(omatch),
          is.data.frame(outcome)||is.numeric(outcome),
          nrow(as.data.frame(outcome))==length(pis),
          !any(is.na(outcome)),
          "SparseM" %in% loadedNamespaces()
          )

if (type!="IPW") stop("Only inverse probability weighting (IPW) for mu_c implemented at present.")

mtd <- matched(omatch)
pis <- pis[mtd]
manyone <- IDmanyonematchedsets(omatch)[mtd]
omatch <- omatch[mtd]
fitwts <- ifelse(manyone, 1-pis, pis)

outcome <- as.data.frame(outcome)[mtd,]
outcome.weights <- 1/(1-pis)
badow <- is.infinite(outcome.weights)
outcome.owted <- as.matrix(outcome*outcome.weights)
outcome.owted.and.centered <-
  SparseM::residuals.slm(SparseM::slm(outcome.owted~omatch, weights=fitwts))

if (is.null(fittedpsmodel)) {thexes <- thexes.centeredandweighted <- matrix(numeric(0),sum(mtd),0)} else
{
  thexes <- model.matrix(fittedpsmodel)[mtd,]
  if (!is.na(icpt <- match("(Intercept)", names(thexes)))) thexes <- thexes[,-icpt]
  thexes.centeredandweighted <-
    fitwts*SparseM::residuals.slm(SparseM::slm(thexes~omatch, weights=fitwts))
}

list(VZtx=crossprod(thexes.centeredandweighted, thexes),
     CZtxZty=crossprod(thexes.centeredandweighted, outcome.owted),
     diag.VZty=colSums(fitwts*outcome.owted.and.centered*outcome.owted)
     )
}

innerSandwichMeat <- # ie, 2*Iobs^-1 - Iobs^-1 %*%I%*%Iobs^-1
  function(fittedpsmodel, esteqncovariances=NULL) # null esteqncovariances indicates unmatched M-estimation
{
stopifnot(is.null(fittedpsmodel) || inherits(fittedpsmodel, "glm"), is.matrix(esteqncovariances))

if (is.null(fittedpsmodel)) return(matrix(numeric(0), 0,0))

inv.i.obs <- GetInvObservedInfo(fittedpsmodel)

if (is.null(esteqncovariances)) return(inv.i.obs)

activevars <- dimnames(inv.i.obs)[[1]]
i <- esteqncovariances[activevars, activevars]
2*inv.i.obs - inv.i.obs%*%i%*%inv.i.obs

}

Vhatmu_c <- function(yc, fittedpsmodel=NULL, omatch=NULL, Tx.group=NULL, type="IPW")
  {
if (!is.null(dim(yc))) {
  dmyc <- dim(yc)
  dnyc <- dimnames(yc)[2L:length(dmyc)]
  dmyc <- dmyc[-1L]
} else { dmyc <- NULL ; dnyc <- NULL}

  yc <- as.data.frame(yc)

stopifnot(!is.null(fittedpsmodel) || !is.null(omatch) || !is.null(Tx.group),
          is.null(fittedpsmodel) ||
          (inherits(fittedpsmodel, "glm") && fittedpsmodel$family$family=="binomial" &&
           nrow(yc)==length(fitted(fittedpsmodel))),
          is.null(omatch) || (inherits(omatch, "optmatch") && nrow(yc)==length(omatch)),
          is.null(Tx.group) || (is.logical(Tx.group) &&nrow(yc)==length(Tx.group)),
          type=="IPW"
          )

if (is.null(fittedpsmodel) && is.null(Tx.group)) Tx.group <- attr(omatch, "contrast.group")
pis <- if (is.null(fittedpsmodel))
  getpisFromPSfit(Tx.group, strata=omatch) else getpisFromPSfit(fittedpsmodel,strata=omatch)
N <- if (is.null(omatch)) length(pis) else sum(matched(omatch))

vs <- covZtx(omatch, fittedpsmodel=fittedpsmodel, outcome=yc,
             pis=pis, type=type)

sm <- innerSandwichMeat(fittedpsmodel, vs$VZtx)
activeXes <- dimnames(sm)[[1]]
###browser()
thevars <- with(vs, diag.VZty - colSums(CZtxZty[activeXes,] * (sm %*% CZtxZty[activeXes,])))
thevars <- thevars/N^2
dim(thevars) <- dmyc
dimnames(thevars) <- dnyc
###thevars

thevars.u <- vs$diag.VZty/N^2
dim(thevars.u) <- dmyc
dimnames(thevars.u) <- dnyc

list(uncorrected=thevars.u, corrected=thevars)
}

Ehatmu_c <- function(yc, fittedpsmodel=NULL, omatch=NULL, type="IPW")
  {
    stopifnot(type=="IPW")
    nobs <- if (is.null(dim(yc))) length(yc) else dim(yc)[1]
    contribs <- if (!is.null(omatch)) matched(omatch) else
    {
      if (!is.null(fittedpsmodel)) fitted(fittedpsmodel)<1 else rep(TRUE, nobs)
    }
if (is.null(dim(yc))) mean(yc[contribs]) else {
  dm <- dim(yc)[-1L]
  dn <- dimnames(yc)[-1L]
  if (is.data.frame(yc)) {yc <- unlist(yc)} else dim(yc) <- NULL
  yc <- yc[contribs]
  dim(yc) <- c(sum(contribs), dm)
  dimnames(yc) <- c(list(1L:sum(contribs)), dn)
  apply(yc, 2:length(dim(yc)), mean)
}
}

hatmu_c <-  function(outcomes, fittedpsmodel=NULL, omatch=NULL, Tx.group=NULL, type="IPW")
  {
if (is.null(dim(outcomes))) outcomes <- as.data.frame(outcomes)


stopifnot(!is.null(fittedpsmodel) || !is.null(omatch) || !is.null(Tx.group),
          is.null(fittedpsmodel) ||
          (inherits(fittedpsmodel, "glm") && fittedpsmodel$family$family=="binomial" &&
           nrow(outcomes)==length(fitted(fittedpsmodel))),
          is.null(omatch) || (inherits(omatch, "optmatch") && nrow(outcomes)==length(omatch)),
          is.null(Tx.group) || (is.logical(Tx.group) &&nrow(outcomes)==length(Tx.group)),
          type=="IPW"
          )

if (is.null(Tx.group))
  {
Tx.group <- if (is.null(fittedpsmodel)) attr(omatch, "contrast.group") else {
  as.logical(fittedpsmodel$y) }
}
ctl.group <- !Tx.group

pis <- if (is.null(fittedpsmodel))
  getpisFromPSfit(Tx.group, strata=omatch) else getpisFromPSfit(fittedpsmodel,strata=omatch)
contribs <- if (is.null(omatch)) (pis<1) else matched(omatch)

N <- sum(contribs)
ctl.group <- ctl.group[contribs]
pis <- pis[contribs]

dm <- dim(outcomes )[-1L]
dn <- c(list(dimnames(outcomes)[[1]][contribs]), dimnames(outcomes)[2L:length(dimnames(outcomes))])
outcomes <- if (length(dm)>1) as.array(outcomes) else as.matrix(outcomes)
dim(outcomes) <- NULL
outcomes <- outcomes[rep(contribs, prod(dm))]
dim(outcomes) <- c(sum(contribs), dm)
dimnames(outcomes) <- dn
apply(outcomes*ctl.group/(1-pis), 2L:length(dim(outcomes)), sum)/N
}

z.muc <- function(yc, fittedpsmodel=NULL, omatch=NULL, withuncorrectedSE=FALSE)
  {
hat.muc <- hatmu_c(outcomes=yc, fittedpsmodel=fittedpsmodel, omatch=omatch)
E.muc <- Ehatmu_c(yc=yc, fittedpsmodel=fittedpsmodel, omatch=omatch)
if (withuncorrectedSE)
  {
SE.muc <- sqrt(as.data.frame(Vhatmu_c(yc=yc, fittedpsmodel=fittedpsmodel, omatch=omatch)))
z.muc <- (hat.muc - E.muc)/SE.muc$c
} else {
SE.muc <- sqrt(Vhatmu_c(yc=yc, fittedpsmodel=fittedpsmodel, omatch=omatch)$c)
z.muc <- (hat.muc - E.muc)/SE.muc
}
p.muc <- pnorm(abs(z.muc), lower=F)
cbind(estimate=hat.muc, E=E.muc, SE=SE.muc, z=z.muc, p.value=p.muc)
  }
