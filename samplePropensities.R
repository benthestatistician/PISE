### DEPENDENCIES: sandwich_glm.R  
samplePropensities <- function(glmm, row.names=NULL,samples=1)
{
covmat <- summary(glmm)$cov.scaled
muhats <- glmm$coeff[!is.na(glmm$coeff)]
stopifnot(length(muhats)==nrow(covmat))
if (!is.null(dimnames(covmat)) & !all(names(muhats)%in%dimnames(covmat)[[1]])) 
  warning("names of covariance matrix don't seem to match coefficient names; proceeding")
if (!is.null(dimnames(covmat)) & all(names(muhats)%in%dimnames(covmat)[[1]])) 
  muhats <- muhats[dimnames(covmat)[[1]]]

glmm.mm <- if (is.null(row.names)) model.matrix(glmm) else {
  model.matrix(glmm)[row.names,,drop=FALSE]}
glmm.mm <- glmm.mm[,names(muhats)]
thetstar <- MASS::mvrnorm(n=samples, mu=muhats, Sigma=covmat)
dim(thetstar) <- c(samples, length(muhats))
nus <- glmm.mm %*% t(thetstar)
##nus <- plogis(nus)
nus <- as.data.frame(nus)
names(nus) <- paste("draw", 1:length(nus), sep="")
row.names(nus) <- row.names(glmm.mm)
nus
}
##' <description>
##' A function to reconstruct \sQuote{posterior draws} of a propensity score from the fitted PS model and from normal deviates.
##' <details>
##' @title Reconstruct posterior draws of a propensity score
##' @param normaldraws Standard normal draws, one for each non-NA coeff, or a data frame of such
##' @param glmm fitted propensity model
##' @param strata if non-null, a stratifying variable; used to calculate PS posterior SD. Defaults to NULL.
##' @param row.names Row names for glmm's model frame, for subsetting it; Defaults to NULL.
##' @param vcov.glmm.eigendecomposition What it sounds like. 
##' @return data frame, the columns of which are sampled propensity scores
##' @author Ben Hansen
reconstructPropensities <- function(normaldraws,glmm,strata=NULL,row.names=NULL, vcov.glmm.eigendecomposition=eigen(SandwichVcov(glmm, strata=strata)))
  {
stopifnot(is.numeric(normaldraws), inherits(glmm, "glm") && glmm$family$family=="binomial")
muhats <- glmm$coeff[!is.na(glmm$coeff)]
coefnames <- names(muhats)
if ((!is.null(dim(normaldraws)) && nrow(normaldraws)!=length(muhats)) ||
    (is.null(dim(normaldraws)) && length(normaldraws)!=length(muhats)))
  stop("size of normaldraws doesn't match number of non-NA coefficients")

normaldraws <- as.matrix(normaldraws)
betahats <- vcov.glmm.eigendecomposition$vectors %*%
  (normaldraws * sqrt(vcov.glmm.eigendecomposition$values))
###thecov <-
###  if (!is.null(strata)) {
###       SandwichVcov(glmm, strata=strata)
###      } else vcov(glmm)
###    thecov <- thecov[coefnames, coefnames]
###eig <- eigen(thecov, symmetric=TRUE)
betahats <- betahats + muhats
row.names(betahats) <- coefnames

glmm.mm <- if (is.null(row.names)) model.matrix(glmm) else {
  model.matrix(glmm)[row.names,,drop=FALSE]}
glmm.mm <- glmm.mm[,coefnames]
nus <- glmm.mm %*% betahats
nus <- as.data.frame(nus)
names(nus) <- paste("draw", 1:length(nus), sep="")
row.names(nus) <- row.names(glmm.mm)
nus
}
