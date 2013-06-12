`lm.ridge.safe` <-
function (formula, data, subset, na.action, lambda = 0, model = FALSE, 
    x = FALSE, y = FALSE, contrasts = NULL, ...) 
{
    m <- match.call(expand.dots = FALSE)
    m$model <- m$x <- m$y <- m$contrasts <- m$... <- m$lambda <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    Y <- model.response(m)
    X <- model.matrix(Terms, m, contrasts)
    n <- nrow(X)
    p <- ncol(X)
    offset <- model.offset(m)
    if (!is.null(offset)) 
        Y <- Y - offset
    if (Inter <- attr(Terms, "intercept")) {
        Xm <- colMeans(X[, -Inter])
        Ym <- mean(Y)
        p <- p - 1
        X <- X[, -Inter] - rep(Xm, rep(n, p))
        Y <- Y - Ym
    }
    else Ym <- Xm <- NA
    Xscale <- drop(rep(1/n, n) %*% X^2)^0.5
    if (any(badXes <- zapsmall(Xscale)==0))
      {
        orig.xnms <- names(Xscale)
        omitted.xnms <- names(Xscale)[badXes]
        X <- X[,!badXes]
        Xscale <- Xscale[!badXes]
        p <- sum(!badXes)
        warning(paste("lm.ridge.safe() dropped", sum(badXes), "of", length(badXes),
                      "indep vars due to lack of variation."))
      }
    X <- X/rep(Xscale, rep(n, p))
    Xs <- try(svd(X))
    if (inherits(Xs, 'try-error'))
      {
        Xs <- svd(X, LINPACK=T)
        Positive <- Xs$d > max(sqrt(.Machine$double.eps)*Xs$d[1], 0)
        p <- sum(Positive)
        Xs$v <- Xs$v[, Positive, drop=FALSE]
        Xs$d <- Xs$d[Positive]
        Xs$u <- Xs$u[,Positive, drop=FALSE]
        warning(paste("lm.ridge.safe() dropped", sum(!Positive), 
                      "degenerate singular values/vectors, of", length(Positive),"."))

      }
    rhs <- t(Xs$u) %*% Y
    d <- Xs$d
    lscoef <- Xs$v %*% (rhs/d)
    lsfit <- X %*% lscoef
    resid <- Y - lsfit
    s2 <- sum(resid^2)/(n - p - Inter)
    HKB <- (p - 2) * s2/sum(lscoef^2)
    LW <- (p - 2) * s2 * n/sum(lsfit^2)
    k <- length(lambda)
    dx <- length(d)
    div <- d^2 + rep(lambda, rep(dx, k))
    a <- drop(d * rhs)/div
    dim(a) <- c(dx, k)
    coef <- Xs$v %*% a
    dimnames(coef) <- list(names(Xscale), format(lambda))
    GCV <- colSums((Y - X %*% coef)^2)/(n - colSums(matrix(d^2/div, 
        dx)))^2
    if (any(badXes))
      {
        bxs <- rep(0, sum(badXes))
        names(bxs) <- omitted.xnms
        Xscale <- c(Xscale, bxs)[orig.xnms]
        coef <- rbind(coef, matrix(bxs, sum(badXes), dim(coef)[2], dimnames=list(omitted.xnms, format(lambda))))
        coef <- coef[orig.xnms,]
      }
    res <- list(coef = drop(coef), scales = Xscale, Inter = Inter, 
        lambda = lambda, ym = Ym, xm = Xm, GCV = GCV, kHKB = HKB, 
        kLW = LW)
    class(res) <- "ridgelm"
    res
}
