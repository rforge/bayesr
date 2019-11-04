exp2 <- function(x) {
  i <- x > 20
  x[i] <- 20
  x <- exp(x)
  x
}

sigmoid.fit <- function(X, y, weights = NULL) {
  if(any(is.na(y)) | any(is.na(X)))
    stop("NA values in data!")
  nr <- nrow(X)
  nc <- ncol(X)
  if(is.null(weights))
    weights <- rep(1, nr)
  objfun <- function(w) {
    sum((weights * (y - (w[1] + w[2] / (1 + exp2(-(X %*% w[-c(1:2)]))))))^2)
  }
  gradfun <- function(w) {
    ez <- exp2(-(X %*% w[-c(1:2)]))
    ez1 <- 1 + ez
    eta <- w[1] + w[2] / ez1
    s1 <- -(2 * (y - eta))
    g <- matrix(0, nrow = nr, ncol = nc + 2)
    g[, 1] <- s1
    g[, 2] <- s1 * 1 / ez1
    for(j in 1:nc) {
      g[, j + 2] <- s1 * w[2] * 1 / ez1 * (1 - ez1) * X[, j]
    }
    colSums(g * weights)
  }
  opt <- optim(rep(0.1, ncol(X) + 2), fn = objfun, method = "L-BFGS-B")
  w <- opt$par
  rval <- list(
    "fitted.values" = drop(w[1] + w[2] / (1 + exp(-(X %*% w[-c(1:2)])))),
    "fit.fun" = function(X) drop(w[1] + w[2] / (1 + exp(-(X %*% w[-c(1:2)])))),
    "coefficients" = w[-c(1:2)],
    "weights" = weights
  )
  rval$rank <- 2 + ncol(X)
  rval$residuals <- y - rval$fitted.values
  class(rval) <- "sigmoid"
  rval
}

logLik.sigmoid <- function(object, ...) {
  res <- object$residuals
  w <- object$weights
  N <- length(res)
  val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * res^2))))
  attr(val, "nall") <- N
  attr(val, "nobs") <- N
  attr(val, "df") <- object$rank + 1
  class(val) <- "logLik"
  val
}

sigmoid <- function(x, ...) {
  UseMethod("sigmoid")
}

sigmoid.formula <- function(formula, data, weights, ..., subset, na.action, contrasts = NULL) 
{
  class.ind <- function(cl) {
    n <- length(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1L:n) + n * (as.vector(unclass(cl)) - 1L)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
  }
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval.parent(m$data))) 
    m$data <- as.data.frame(data)
  m$... <- m$contrasts <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  x <- model.matrix(Terms, m, contrasts)
  cons <- attr(x, "contrast")
  w <- model.weights(m)
  if(length(w) == 0L) 
    w <- rep(1, nrow(x))
  y <- model.response(m)
  res <- sigmoid.default(x, y, w, ...)
  res$terms <- Terms
  res$coefnames <- colnames(x)
  res$call <- match.call()
  res$na.action <- attr(m, "na.action")
  res$contrasts <- cons
  res$xlevels <- .getXlevels(Terms, m)
  class(res) <- c("sigmoid.formula", "sigmoid")
  res
}

sigmoid.default <- function(x, y, weights, ...) {
  sigmoid.fit(x, y, weights)
}

estfun.sigmoid <- function(x, ...) {
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  if(any(alias <- is.na(coef(x)))) 
    xmat <- xmat[, !alias, drop = FALSE]
  wts <- weights(x)
  if(is.null(wts)) 
    wts <- 1
  res <- residuals(x)
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  return(rval)
}


## Trying to split.
split_fun <- function(x, y, X, weights = rep(TRUE, length(y)), plot = FALSE, smax = 5)
{
  warn <- getOption("warn")
  options("warn" = -1)
  on.exit(options("warn" = warn))
  np <- ncol(X) + 2
  splits <- sort(unique(x[weights]))
  if(length(splits) > smax) {
    splits <- quantile(x[weights], prob = seq(0, 1, length = smax))
  }
  sse <- rep(Inf, length(splits))
  for(i in 1:length(splits)) {
    if(i > 1 & i < length(splits)) {
      j1 <- (x <= splits[i]) & weights
      j2 <- (x > splits[i]) & weights
      if((sum(j1) > np) & (sum(j2) > np)) {
        b1 <- sigmoid.fit(X[j1,,drop=FALSE], y[j1])
        b2 <- sigmoid.fit(X[j2,,drop=FALSE], y[j2])
        r1 <- diff(range(fitted(b1)))
        r2 <- diff(range(fitted(b2)))
        sse[i] <- -1 * (r1 + r2)
      }
    }
  }
  split_at <- splits[i <- which.min(sse)]
  j1 <- (x <= split_at) & weights
  j2 <- (x > split_at) & weights
  b1 <- sigmoid.fit(X[j1,,drop=FALSE], y[j1])
  b2 <- sigmoid.fit(X[j2,,drop=FALSE], y[j2])

  if(plot) {
    plot(y ~ x)
    points(x[j1], y[j1], col = 4, pch = 16)
    points(x[j2], y[j2], col = 2, pch = 16)
    plot2d(fitted(b1) ~ x[j1], add = TRUE, col.lines = 4, lwd = 2)
    plot2d(fitted(b2) ~ x[j2], add = TRUE, col.lines = 2, lwd = 2)
  }

  cb1 <- coef(b1)
  cb2 <- coef(b2)
  return(list("sse" = sse[i], "split" = split_at,
    "coefficients" = list("L" = cb1, "R" = cb2),
    "weights" = list("L" = j1, "R" = j2), "b1" = b1, "b2" = b2))
}

ntree <- function(x, y, k = 20, smax = 5, verbose = TRUE, plot = TRUE, ...) {
  if(is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
    X <- cbind(1, x)
  } else {
    x <- as.matrix(x)
    has_itcpt <- apply(x, 2, function(x) all(x == 1) )
    if(!any(has_itcpt)) {
      X <- cbind(1, x)
    } else {
      x <- x[, -which(has_itcpt), drop = FALSE]
      X <- cbind(1, x)
    }
  }
  y <- scale2(y, 0.01, 0.99)
  if(plot)
    plot(y, ylim = c(-1, 1), main = paste(0, "node"))
  np <- ncol(X) + 2
  splits <- list()
  for(l in 1:ncol(x)) {
    splits[[l]] <- split_fun(x[, l], y, X, smax = smax)
  }
  sse <- sapply(splits, function(x) { x$sse })
  split_take <- splits[[which.min(sse)]]
  fit <- rep(0, length(y))
  y[split_take$weights$L] <- y[split_take$weights$L] - fitted(split_take$b1)
  y[split_take$weights$R] <- y[split_take$weights$R] - fitted(split_take$b2)
  fit[split_take$weights$L] <- fit[split_take$weights$L] + fitted(split_take$b1)
  fit[split_take$weights$R] <- fit[split_take$weights$R] + fitted(split_take$b2)
  i <- 1
  if(plot)
    plot(y, ylim = c(-1, 1), main = paste(1, "node"))
  coef <- do.call("cbind", split_take$coefficients)
  lw <- as.data.frame(do.call("cbind", split_take$weights))
  while(i < k) {
    if(verbose)
      cat(".. node", i, "of", k, "\n")
    splits <- list()
    for(j in 1:ncol(lw)) {
      if(sum(lw[, j]) > np) {
        splitsj <- list()
        for(l in 1:ncol(x)) {
          splitsj[[l]] <- split_fun(x[, l], y, X, weights = lw[, j], smax = smax)
        }
        l <- which.min(sapply(splitsj, function(x) { x$sse }))
        splits[[j]] <- splitsj[[l]]
      }
    }
    if(length(splits) > 0) {
      ssej <- sapply(splits, function(x) { x$sse })
      if(!all(!is.finite(ssej))) {
        j <- which.min(ssej)
        split_take <- splits[[j]]
        lw[, j] <- do.call("cbind", split_take$weights)
        lw <- as.data.frame(as.matrix(lw))
        coef <- cbind(coef, do.call("cbind", split_take$coefficients))
        y[split_take$weights$L] <- y[split_take$weights$L] - fitted(split_take$b1)
        y[split_take$weights$R] <- y[split_take$weights$R] - fitted(split_take$b2)
        fit[split_take$weights$L] <- fit[split_take$weights$L] + fitted(split_take$b1)
        fit[split_take$weights$R] <- fit[split_take$weights$R] + fitted(split_take$b2)
        if(plot)
          plot(y, ylim = c(-1, 1), main = paste(i, "nodes"))
        i <- ncol(coef)
      } else {
        i <- k
      }
    } else {
      i <- k
    }
  }
  i <- ncol(coef)
  if(i < k) {
    coef <- cbind(coef, build_net_w(X, y, k = k - i, ..., I = i, plot = plot))
  }
  i <- ncol(coef)
  if(verbose)
    cat(".. node", i, "of", k, "\n")
  return(coef)
}

build_net_w <- function(X, y, k = 10, n = 10, plot = FALSE, eps = 0.3, ...) {
  I <- list(...)$I
  if(is.null(I))
    I <- 1
  ind <- 1:nrow(X)
  tX <- t(X)
  err0 <- 1e+20
  w <- NULL
  i <- 0
  iter <- 1
  while(i < k) {
    j <- sample(ind, size = 1)
    tx <- as.numeric(X[j, , drop = FALSE])
    cs <- colSums((tX - tx)^2)
    if(length(n) < 2)
      n2 <- n
    else
      n2 <- sample(n, size = 1)
    take <- order(cs)[1:n2]
    yn <- y[take]
    xn <- X[take, , drop = FALSE]
    m <- sigmoid.fit(xn, yn)
    err <- sum((y[take] - fitted(m))^2)
    erry <- sum(y[take]^2)
    epsr <- (err - erry) / erry
    if(epsr < -1 * eps) {
      y[take] <- y[take] - fitted(m)
      if(plot)
        plot(y, main = paste(I, "nodes"), ylim = c(-1, 1))
      w <- cbind(w, coef(m))
      i <- ncol(w)
      I <- I + 1
    }
    eps <- eps * 0.99
    iter <- iter + 1
    if(iter > 100 * k) {
      i <- k
      stop("could not compute all weights, set argument eps!")
    }
  }
  return(w)
}

