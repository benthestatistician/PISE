ett <- function(x,y) x
harmonic <-  function (x, y) 2*(1/x + 1/y)^-1
##'
##' ##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' - Build the model frame as in `lm`, except use `na.pass`
##' - Look for a column inheriting from class `optmatch`.  (If there are none, or two or more, then bail.)
##' - use `model.offset` to assemble vector of offsets (if applicable)
##' - Strip the optmatch variable from the provided formula and push it, data frame through `fill.NAs` (for now)
##' - Remove from resulting model frame rows that have NAs on either the response or the optmatch or the offset.
##' - Remove corresponding entries from the optmatch and from the offset
##' - Build model matrix, then linearly transform into matched differences model matrix, `x`
##' - Also build vector `y` of matched diffs on response, vectors of n.tx, n.ctl
##' - Strip out rows/entries for which n.tx=0 or n.ctl=0
##' - Look in x for a column of 1's, other than "(Intercept)".  If found, remove "(Intercept)" column.  If not found, relabel "(Intercept)" as "z"
##' - fitweights = weights(n.tx, n.ctl)
##' - wlm.fit(x, y, w=fitweights)

##' The matched differences model matrix has a col for each covariate but w/
##' one row for each matched set.  Entries in this row are per matched set differences of means between treatment and control groups, with the treatment group being identified by the `contrast.group` attribute of the optmatch object.  
##'
##' @title Ordinary least squares for matched differences
##' @param formula 
##' @param data 
##' @param weights Numeric vector-valued function of 2 vector args `n.t`, `n.c`
##' @param ... additional arguments passed to `model.frame`
##' @return object of class `lm`
##' @author Ben B Hansen
mlm <- function(formula, data, weights=c("ett", "harmonic"), ...)
    {

    }
