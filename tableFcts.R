###summaryTable <- function(x,...) UseMethod("summaryTable")
summaryTable.rdObj <- function(rdo, getZstatsFct= function(x) x,
                               propensity=NULL,...)
{
tptiles <- qnorm(c(.05, .95))
tptile.labels <- c("\\underline{\\PP\\{t_Z < z_{.05}\\}}","\\underline{\\PP\\{z_{.95}<t_Z\\}}")
tconts <- list(Row0="\\begin{array}{c@{\\hspace{1.5\\tabcolsep}}c}") # Could add a table label here

themean <- estimateMean(rdo, FUN=getZstatsFct, propensity=propensity)
thesd <- estimateSD(rdo, FUN=getZstatsFct, propensity=propensity)
digits.mn <- pmax(1, pmin(floor(-log10(themean$se)), 3))
digits.sd <- pmax(1, pmin(floor(-log10(thesd$se)), 3))
tconts$Row1 <- paste(paste("\\underline{\\EE t_Z}", "\\underline{\\var t_Z}", sep="&"), "\\\\ ")
tconts$Row2 <- paste(paste(formatC(round(themean$muhat, digits.mn), digits=digits.mn, format='f'),
                           formatC(round(thesd$sigmahat, digits.sd), digits=digits.sd, format='f'),
                           sep="&"),
                 "\\\\")
               
#browser()
getZstatsFct2 <- function(x) {
                   c(x<tptiles[1], x>tptiles[2])+
                     .5*(x==tptiles[1:2])
                 }
thesizes <-
    estimateMean(rdo,FUN=getZstatsFct2, propensity=propensity)
digits.sz <- max(1, min(floor(-log10(thesizes$se)), 3))
tconts$Row3 <- paste(paste(tptile.labels, collapse="&"), "\\\\ ")
tconts$Row4 <- paste(paste(round(thesizes$muhat, digits.sz), collapse="&"), "\\\\")
tconts$Row5 <- "\\end{array}\n"
paste(tconts, collapse="\n ")
}


### Tests -- run after simulateMSshuffle-tests.R
### nullSEs <- (sapply(as.data.frame(R), function(x) sumStatVar(x, attr(fm3, "contrast.group"), fm3))/
### cat(summaryTable.rdObj(rdObj2, getZstatsFct=function(x) x['te.m100']/nullSEs['te.m100']))
