ett <- function(n.t,n.c) n.t
harmonic <-  function (n.t, n.c) 2*(1/n.t + 1/n.c)^-1
##'
##' ##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' - Build the model frame as in `lm`, except use `na.pass`
##' - Look for a column inheriting from class `optmatch`.  (If there are none, or two or more, then bail.)
##' - use `model.response`, `model.offset` (if applicable), `model.weights` (if applicable) to assemble corresponding vectors
##' - Strip the optmatch variable from the provided formula and push reduced formula, model frame through `fill.NAs` (for now)
##' - Remove from resulting model frame rows that have NAs on either the response or the optmatch or the offset.
##' - Remove corresponding entries from the optmatch, the response and (if applicable) the offset and the weights
##' - Build model matrix, then linearly transform into matched differences model matrix, `x`
##' - Also build vectors `y` of matched diffs on response, `offset` of matched diffs in offset
##' - build vectors of n.tx, n.ctl.  If weights provided, these are sums of weights.
##' - Strip out rows/entries for which n.tx=0 or n.ctl=0
##' - Look in x for a column of 1's, other than "(Intercept)".  If found, remove "(Intercept)" column.  If not found, relabel "(Intercept)" as "z"
##' - fitweights = ms.weights(n.tx, n.ctl)  If any of fitweights are NA, bail.  
##' - If fit.type="lm", then do and return wlm.fit(x, y, w=fitweights, offset=offset)
##' - If/when we enable fit.type="robust", call rlm on x, y, w, offset
##' The matched differences model matrix has a col for each covariate but w/
##' one row for each matched set.  Entries in this row are per matched set differences of means between treatment and control groups, with the treatment group being identified by the `contrast.group` attribute of the optmatch object.  
##'
##' @title Ordinary least squares for matched differences
##' @param formula 
##' @param data 
##' @param ms.weights Function of 2 vector args `n.t`, `n.c`, sums of weights from treatment and control group members by matched set, returning vector of matched-set specific weights. `harmonic` returns harmonic means of n.t and n.c; `ett` simply returns n.t.
##' @param fit.type character string indicating type of fit. For now, only "lm", but may expand to include rlm
##' @param fit.control optional list of additional arguments to the fitter
##' @param ... additional arguments passed to `model.frame`
##' @return object of class `lm`
##' @author Ben B Hansen
mlm <- function(formula, data, ms.weights=c(ett, harmonic), fit.type="lm", fit.type=list(NULL),...)
    {

    }
