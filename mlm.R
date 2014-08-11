library(optmatch)
library(SparseM)

##'
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' - Build the model frame as in `lm`, except use `na.pass`
##' - Look for a column inheriting from class `optmatch`.  (If there are none, or two or more, then bail.)
##' - use `model.response`, `model.offset` (if applicable), `model.weights` (if applicable) to assemble corresponding vectors
##' - Strip the optmatch variable from the provided formula and push reduced formula, model frame through `fill.NAs` (for now)
##' - Remove from resulting model frame rows that have NAs on either the response or the optmatch or the offset.
##' - Remove corresponding entries from the optmatch, the response and (if applicable) the offset and the weights
##' - Coerce the optmatch to a matrix.csr
##' - Build model matrix, then use matrix.csr version of the optmatch to linearly transform into matched differences model matrix, `x`
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
mlm <- function(formula, data, ms.weights = ett, fit.type = "lm", fit.control = list(NULL), na.action = na.pass, ...) {
  parsed <- parseMatchingProblem(formula, data, na.action, ...)

  outcome <- model.response(parsed$mf)
  weights <- model.weights(parsed$mf)
  offset <- model.offset(parsed$mf)
  
  # this will be a nearly filled in model matrix (ie. all factors expanded), but without an intercept
  noNAs <- fill.NAs(parsed$fmla, parsed$mf)

  theMatch <- parsed$match

  checkNA <- function(i) {
    if (is.null(i)) {
      return(FALSE)
    }
    return(is.na(i))
  }

  remove <- !with(parsed,
                  checkNA(weights) |
                  is.na(outcome) |  
                  is.na(theMatch))

  noNAs <- model.matrix(parsed$fmla, noNAs[remove, , drop = FALSE])
  if ("(Intercept)" %in% colnames(noNAs)) {
    noNAs <- noNAs[, -which("(Intercept)" %in% colnames(noNAs)), drop = FALSE]
  }

  theMatch <- theMatch[remove]
  outcome <- outcome[remove]
  weights <- weights[remove]

  z <- attr(theMatch, "contrast.group")
  nt <- sapply(levels(theMatch), function(l) { sum(theMatch == l & z) })
  nc <- sapply(levels(theMatch), function(l) { sum(theMatch == l & !z) })

  missingTorC <- nt == 0 | nc == 0 
  nt <- nt[!missingTorC]
  nc <- nc[!missingTorC]

  matchCsr <- as(theMatch, "matrix.csr")
  matchCsr <- matchCsr[!missingTorC, ]

  if (ncol(noNAs) == 0) {
    X <- matrix(1, ncol = 1, nrow = nlevels(theMatch), dimnames = list(levels(theMatch), "z"))
  } else {
    # make the design matrix for the matched sets
    # switching back to dense representation since the we don't expect many zero's in the design matrix
    X <- as.matrix(matchCsr %*% as.matrix(noNAs)) 
    colnames(X) <- colnames(noNAs)
    # the rows are the matched sets (but we don't need to include those)

    # add treatment indicator if there is not one
    # an intercept column will turn up as all zeros
    isConst <- apply(X, 2, function(col) { all(col == col[1]) })
    if (!any(isConst)) {
      X <- cbind(X, z = 1)
    }
  }

  Y <- matchCsr %*% outcome

  fit.weights <- ms.weights(nt, nc)


  if ((fit.type == "robust" || fit.type == "rlm")) {
    if(require(MASS)) {
      return(rlm(X, Y, weights = fit.weights))
    } else {
      warning("MASS package not found, cannot use robust regression fit.") # this should probably never happen. MASS in in base now.
    }
  }

  # can't fit with rlm package. fall back to good old lm
  lm.wfit(X, Y, w = fit.weights, offset = offset)
}


# This next line should get put in the makeOptmatch.R file. It provides S4 compatability.
setOldClass(c("optmatch", "factor"))

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Sparse matrices with which to assemble treatment minus control differences by matched set
##' @param from An optmatch object 
##' @return A matrix.csr object by which to left-multiply vectors
##' and model matrices in order to assemble matched differences.
##' @author Ben B Hansen
setAs("optmatch", "matrix.csr", function(from) {
  # treatment variable, a logical
  zz <- optmatch:::toZ(attr(from, "contrast.group")) # can remove the explicit namespace when this goes in the optmatch pkg

  # list of positions of treatment member(s), then
  # control group members; by matched set
  pos.tc <- lapply( levels(from), function(lev) c(which(from==lev & zz),
                                                  which(from==lev & !zz)))

  # starting positions for rows of the csr matrix
  rowstarts <- 1 + c(0, cumsum(sapply(pos.tc, length)))
  rowstarts <- as.integer(rowstarts)

  # each row has 1st tx and then ctl, but we need to know how many of each      
  n.t <- as.integer(table(from[zz, drop = FALSE]))
  n.c <- as.integer(table(from))-n.t

  # if either t or c unrepresented in a matched set, null out other group's contrib
  # (this can happen due to `from` having been subsetted, perhaps b/c of NAs elsewhere)
  tscale <- ifelse(n.c&n.t, 1/n.t, 0)
  cscale <- ifelse(n.c&n.t, -1/n.c, 0)

  # multipliers to go in positions pos.tc
  multipliers <- rep(as.vector(rbind(tscale, cscale)),
                     as.vector(rbind(n.t, n.c)) )

  new("matrix.csr",
      ra = multipliers,
      ja = unlist(pos.tc),
      ia = rowstarts,
      dimension = c(nlevels(from), length(from)))
})

ett <- function(n.t,n.c) n.t
harmonic <-  function (n.t, n.c) 2*(1/n.t + 1/n.c)^-1

#' Helper to parse a matched analysis from a formula and a data.frame.
#'
#' @param formula The formula.
#' @param data The data.frame containing terms in the formula.
#' @param ... Other arguments passed to `model.frame`.
#' @return A list with: `mf` a model frame stripped of the optmatch argument, `match` the matched factor
parseMatchingProblem <- function(formula, data, na.action = na.pass, ...) {
  mf <- model.frame(formula, data, na.action = na.action, ...)

  isMatch <- sapply(mf, function(i) { inherits(i, "optmatch") })

  if (sum(isMatch) != 1) {
    stop("You must include precisely one matching in the formula.")
  }

  match <- mf[, isMatch, drop = TRUE]
  names(match) <- rownames(data)

  mname <- colnames(mf)[isMatch]

  # now that we've peeled out the match, the match out of the formula
  # I dislike string hacking on formulas, but I can't think of a better way to do this.
  fparts <- as.character(formula)
  fstr <- strsplit(fparts[3], " \\+ ")[[1]]
  keep <- grep(pattern = mname, fstr, value = TRUE, invert = TRUE)

  if (length(keep) == 0) {
    keep <- "1"
  }
  newf <- as.formula(paste(fparts[2], "~", paste(keep, collapse = "+")))

  # now make a new model frame, using the reduce form of the formula
  mf <- model.frame(newf, data, na.action = na.action, ...)
  
  return(
      list(
          fmla = newf,
          mf = mf,
          match = match))
}
