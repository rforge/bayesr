## Some useful links:
## http://adv-r.had.co.nz/C-interface.html
## http://stackoverflow.com/questions/7457635/calling-r-function-from-c
## http://gallery.rcpp.org/articles/r-function-from-c++/
## BayesG propose function in C.
propose_iwls <- function(x, family,
  response, eta, id, rho, ...)
{
  .Call("do_propose", x, family, response, eta, id, rho)
}

propose_iwls0 <- function(x, family, response, eta, id, ...)
{
  require("mvtnorm")

  args <- list(...)

  if(!is.null(args$no.mcmc)) {
    if(!x$fixed & is.null(x$sp)) {
      for(j in seq_along(x$S)) {
        a <- x$rank[j] / 2 + x$a
        b <- 0.5 * crossprod(x$state$g, x$S[[j]]) %*% x$state$g + x$b
        x$state$tau2[j] <- 1 / rgamma(1, a, b)
      }
    }
  }

  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute weights.
  weights <- family$weights[[id]](response, peta)

  ## Score.
  score <- family$score[[id]](response, peta)

  ## Compute working observations.
  z <- eta[[id]] + 1 / weights * score

  ## Compute old log likelihood and old log coefficients prior.
  pibeta <- family$loglik(response, peta)
  p1 <- x$prior(x$state$g, x$state$tau2)

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit

  ## Compute mean and precision.
  S <- 0
  XW <- t(x$X * weights)
  P <- if(x$fixed) {
    if(k <- ncol(x$X) < 2) {
      1 / (XW %*% x$X)
    } else chol2inv(chol(XW %*% x$X))
  } else {
    for(j in seq_along(x$S))
      S <- S + 1 / x$state$tau2[j] * x$S[[j]]
    chol2inv(chol(XW %*% x$X + S))
  }
  P[P == Inf] <- 0
  M <- P %*% (XW %*% (z - eta[[id]]))

  ## Save old coefficients
  g0 <- drop(x$state$g)

  ## Sample new parameters.
  x$state$g <- drop(rmvnorm(n = 1, mean = M, sigma = P))

  ## Compute log priors.
  p2 <- x$prior(x$state$g, x$state$tau2)
  qbetaprop <- dmvnorm(x$state$g, mean = M, sigma = P, log = TRUE)

  ## Compute fitted values.        
  x$state$fit <- drop(x$X %*% x$state$g)

  ## Set up new predictor.
  eta[[id]] <- eta[[id]] + x$state$fit

  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute new log likelihood.
  pibetaprop <- family$loglik(response, peta)

  ## Compute new weights
  weights <- family$weights[[id]](response, peta)

  ## New score.
  score <- family$score[[id]](response, peta)

  ## New working observations.
  z <- eta[[id]] + 1 / weights * score

  ## Compute mean and precision.
  XW <- t(x$X * weights)
  P2 <- if(x$fixed) {
    if(k < 2) {
      1 / (XW %*% x$X)
    } else chol2inv(chol(XW %*% x$X))
  } else {
    chol2inv(L <- chol(P0 <- XW %*% x$X + S))
  }
  P2[P2 == Inf] <- 0
  M2 <- P2 %*% (XW %*% (z - (eta[[id]] - x$state$fit)))

  ## Get the log prior.
  qbeta <- dmvnorm(g0, mean = M2, sigma = P2, log = TRUE)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & is.null(args$no.mcmc)) {
    if(!x$fixed & is.null(x$sp)) {
      for(j in seq_along(x$S)) {
        a <- x$rank[j] / 2 + x$a
        b <- 0.5 * crossprod(x$state$g, x$S[[j]]) %*% x$state$g + x$b
        x$state$tau2[j] <- 1 / rgamma(1, a, b)
      }
    }
  }

  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

  return(x$state)
}


## Sampling with multivariate normal proposals.
propose_mvn <- function(x, family, response, eta, id, ...)
{
  require("mvtnorm")
  args <- list(...)

  if(is.null(x$state$XX)) x$state$XX <- crossprod(x$X)

  ll1 <- family$loglik(response, family$map2par(eta))
  p1 <- x$prior(x$state$g, x$state$tau2)
  eta[[id]] <- eta[[id]] - x$state$fit
  P <- matrix_inv(x$state$XX + 1 / x$state$tau2 * x$S[[1]])

  if((x$state$iter < x$adapt) & x$xt$adaptive) {
    eta2 <- eta
    do <- TRUE
    while(do) {
      g <- drop(rmvnorm(n = 1, mean = x$state$g, sigma = x$state$scale[1] * P))
      f <- x$get.mu(x$X, g)
      if(!is.null(x$xbin.ind))
        f <- f[x$xbin.ind]
      eta2[[id]] <- eta[[id]] + f
      ll2 <- family$loglik(response, family$map2par(eta2))
      p2 <- x$prior(g, x$state$tau2)
      alpha <- drop((ll2 + p2) - (ll1 + p1))
      accepted <- if(is.na(alpha)) FALSE else log(runif(1)) <= alpha
      x$state$scale[1] <- if(alpha < log(0.23)) {
        x$state$scale[1] - 0.1 * x$state$scale[1]
      } else {
        x$state$scale[1] + 0.1 * x$state$scale[1]
      }
      if(accepted) do <- FALSE
    }
  }

  if(x$state$iter %% x$xt$step == 0) {
    if(!is.null(x$state$mode)) {
      x$state$g <- x$state$mode$g
    } else {
      opt <- update_optim2(x, family, response, eta, id, ...)
      x$state$g <- opt$g
    }
  }

  x$state$g <- drop(rmvnorm(n = 1, mean = x$state$g, sigma = x$state$scale[1] * P))
  x$state$fit <- x$get.mu(x$X, x$state$g)
  eta[[id]] <- eta[[id]] + x$state$fit
  ll2 <- family$loglik(response, family$map2par(eta))
  p2 <- x$prior(x$state$g, x$state$tau2)

  if(!x$fixed & is.null(x$sp) & is.null(args$no.mcmc)) {
    a <- x$rank / 2 + x$a
    b <- 0.5 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g + x$b
    x$state$tau2 <- 1 / rgamma(1, a, b)
  }

  x$state$alpha <- drop((ll2 + p2) - (ll1 + p1))

  return(x$state)
}


logPost <- function(g, x, family, response, eta, id)
{
  ## Set up new predictor.
  f <- x$get.mu(x$X, g)
  if(!is.null(x$xbin.ind))
    f <- f[x$xbin.ind]
  eta[[id]] <- eta[[id]] + f

  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute log likelihood and log coefficients prior.
  ll <- family$loglik(response, peta)
  lp <- x$prior(x$state$g, x$state$tau2)

  -1 * (ll + lp)
}


## Numerical derivatives.
num_deriv2 <- function(y, eta, family, id = NULL,
  d = 1, method = "simple", eps = 1e-04)
{
  require("numDeriv")

  if(is.null(id)) id <- family$names[1L]

  d1fun <- function(x, y, eta) {
    eta[[id]] <- x
    family$d(y, eta, log = TRUE)
  }

  d2fun <- function(x, y, eta) {
    grad(d1fun, x, y = y, eta = eta, method = method, method.args = list("eps" = eps))
  }

  dp <- grad(if(d < 2) d1fun else d2fun, eta[[id]],
    y = y, eta = eta, method = method,
    method.args = list("eps" = eps))

  return(dp)
}

num_deriv <- function(y, eta, family, id = NULL, d = 1, eps = 1e-04)
{
  d1 <- function(fn, x) {
    fn1 <- fn(x + eps)
    fn2 <- fn(x)
    d <- (fn1 - fn2) / eps
    d
  }

  d2 <- function(fn, x){
    fn1 <- d1(fn, x + eps)
    fn2 <- d1(fn, x)
    d <- (fn1 - fn2) / eps
    d
  }

  fun <- function(x) {
    eta[[id]] <- x
    family$d(y, eta, log = TRUE)
  }
 
  d <- if(d < 2) {
    d1(fun, eta[[id]])
  } else {
    optimHess(eta[[id]], fun)
  }

  return(d) 
}


## Slice sampling
## See: http://www.cs.toronto.edu/~radford/ftp/slice-R-prog
## Log-posterior used by propose_slice().
logPost2 <- function(g, x, family, response, eta, id)
{
  f <- x$get.mu(x$X, g)
  if(!is.null(x$xbin.ind))
    f <- f[x$xbin.ind]
  eta[[id]] <- eta[[id]] + f
  ll <- family$loglik(response, family$map2par(eta))
  lp <- x$prior(g, x$state$tau2)
  return(ll + lp)
}

## For rational splines.
logPost3 <- function(g, x, family, response, eta, id)
{
  f <- x$get.mu(x$X, g)
  if(!is.null(x$xbin.ind))
    f <- f[x$xbin.ind]
  eta[[id]] <- eta[[id]] + f
  ll <- family$loglik(response, family$map2par(eta))
  lp <- x$prior(g, x$state$tau2)
  return(ll + lp)
}

## For variance parameter sampling.
logPost4 <- function(tau2, x, family, response, eta, id)
{
  ll <- family$loglik(response, family$map2par(eta))
  lp <- x$prior(x$state$g, tau2)
  return(ll + lp)
}

## Rational splines variance parameter.
logPost5 <- function(tau2, x, family, response, eta, id)
{
  ll <- family$loglik(response, family$map2par(eta))
  lp <- x$prior(x$state$g, tau2)
  return(ll + lp)
}

## Univariate slice sampling.
uni.slice2 <- function(g, x, family, response, eta, id, j,
  w = 1, m = 100, lower = -Inf, upper = +Inf, logPost, rho)
{
  .Call("uni_slice", g, x, family, response, eta, id,
     as.integer(j), w, as.integer(m), lower, upper, logPost, rho)
}

uni.slice <- function(g, x, family, response, eta, id, j, ...,
  w = 1, m = 30, lower = -Inf, upper = +Inf, logPost)
{
  x0 <- g[j]
  gL <- gR <- g

  gx0 <- logPost(g, x, family, response, eta, id)

  ## Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)

  ## Find the initial interval to sample from.
  ## w <- w * abs(x0) ## FIXME???
  u <- runif(1, 0, w)
  gL[j] <- g[j] - u
  gR[j] <- g[j] + (w - u)  ## should guarantee that g[j] is in [L, R], even with roundoff

  ## Expand the interval until its ends are outside the slice, or until
  ## the limit on steps is reached.
  if(is.infinite(m)) {
    repeat {
      if(gL[j] <= lower) break
      if(logPost(gL, x, family, response, eta, id) <= logy) break
      gL[j] <- gL[j] - w
    }
    repeat {
      if(gR[j] >= upper) break
      if(logPost(gR, x, family, response, eta, id) <= logy) break
      gR[j] <- gR[j] + w
    }
  } else {
    if(m > 1) {
      J <- floor(runif(1, 0, m))
      K <- (m - 1) - J
      while(J > 0) {
        if(gL[j] <= lower) break
        if(logPost(gL, x, family, response, eta, id) <= logy) break
        gL[j] <- gL[j] - w
        J <- J - 1
      }
      while(K > 0) {
        if(gR[j] >= upper) break
        if(logPost(gR, x, family, response, eta, id) <= logy) break
        gR[j] <- gR[j] + w
        K <- K - 1
      }
    }
  }

  ## Shrink interval to lower and upper bounds.
  if(gL[j] < lower) {
    gL[j] <- lower
  }
  if(gR[j] > upper) {
    gR[j] <- upper
  }

  ## Sample from the interval, shrinking it on each rejection.
  repeat {
    g[j] <- runif(1, gL[j], gR[j])

    gx1 <- logPost(g, x, family, response, eta, id)

    if(gx1 >= logy) break

    if(g[j] > x0) {
      gR[j] <- g[j]
    } else {
      gL[j] <- g[j]
    }
  }

  ## Return the point sampled
  return(g)
}


## Actual univariate slice sampling propose() function.
propose_slice <- function(x, family,
  response, eta, id, rho, ...)
{
  args <- list(...)

  if(!is.null(args$no.mcmc)) {
    for(j in seq_along(x$state$tau2)) {
      x$state$tau2 <- uni.slice(x$state$tau2, x, family, response, eta, id, j,
        logPost = logPost4, rho = rho, lower = 0)
    }
  }

  ## Remove fitted values.
  eta[[id]] <- eta[[id]] - x$state$fit

  for(j in seq_along(x$state$g)) {
    x$state$g <- uni.slice(x$state$g, x, family, response, eta, id, j,
      logPost = logPost2, rho = rho)
  }

  ## Setup return state.
  x$state$alpha <- log(1)
  x$state$fit <- x$get.mu(x$X, x$state$g)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & is.null(args$no.mcmc)) {
    for(j in seq_along(x$state$tau2)) {
      x$state$tau2 <- uni.slice(x$state$tau2, x, family, response, eta, id, j,
        logPost = logPost4, rho = rho, lower = 0)
    }
  }

  return(x$state)
}


propose_slice2 <- function(x, family,
  response, eta, id, rho, ...)
{
  args <- list(...)

  if(!is.null(args$no.mcmc)) {
    for(j in seq_along(x$state$tau2)) {
      x$state$tau2 <- uni.slice(x$state$tau2, x, family, response, eta, id, j,
        logPost = logPost4, rho = rho, lower = 0)
    }
  }

  ## Remove fitted values.
  eta[[id]] <- eta[[id]] - x$state$fit

  x$state$g <- stepout.slice.sample(logPost2, x$state$g, sample.size = 2, tuning = 1,
    step.out = TRUE, limit = length(x$state$g) * 100,
    response = response, eta = eta, id = id, x = x, family = family)

  ## Setup return state.
  x$state$alpha <- log(1)
  f <- x$get.mu(x$X, x$state$g)
  if(!is.null(x$xbin.ind))
    f <- f[x$xbin.ind]
  x$state$fit <- f

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & is.null(args$no.mcmc)) {
    for(j in seq_along(x$state$tau2)) {
      x$state$tau2 <- uni.slice(x$state$tau2, x, family, response, eta, id, j,
        logPost = logPost4, rho = rho, lower = 0)
    }
  }

  return(x$state)
}


stepout.slice.sample <- function(logPost, g0, sample.size, tuning = 1,
  step.out = TRUE, limit = length(g0) * 100, ...) {
  p <- length(g0)
  X <- array(NA, c(sample.size, p))
  X[1, ] <- g0
  nevals <- 0

  for(obs in 2:sample.size) {
    ## Set g0 to current state and y.slice to a new slice level.
    g0 <- X[obs-1,,drop=TRUE]
    y.slice <- logPost(g0, ...) - rexp(1)
    nevals <- nevals + 1

    ## Transition each coordinate in turn.
    for(i in 1:p) {
      ## If step.out is set, find the boundaries of the slice in
      ## this coordinate using stepout.one.coord.  Otherwise, choose
      ## a random interval.
      if(step.out) {
	      lr <- stepout.one.coord(logPost, g0, i, y.slice,
          tuning, limit * sample.size - nevals, ...)
        nevals <- nevals + lr$nevals
        if(is.null(lr$left))
          return(list(X = X[1:(obs - 1), , drop = FALSE], evals = nevals))
      } else {
        lr <- list(left = g0[i] - runif(1) * tuning)
        lr$right <- lr$left + tuning
      }

      ## Draw proposals, shrinking the slice estimate along the way,
      ## until a proposal is accepted.
      repeat {
	      ## Draw a proposal by modifying the current coordinate of
	      ## g0 to be uniformly chosen from the current slice estimate..
        x1 <- g0
        x1[i] <- lr$left + (lr$right-lr$left) * runif(1)

	      ## If the proposal is in the slice, accept it and move on
	      ## to the next coordinate.
        y1 <- logPost(x1, ...)
        nevals <- nevals + 1
        if(y1 >= y.slice) {
          g0 <- x1
          break
        }

        ## Otherwise, shrink the slice towards g0.
        if(x1[i] < g0[i])
          lr$left <- x1[i]
        else
          lr$right <- x1[i]
      }
    }

    ## Having updated each coordinate, save the accepted proposal.
    X[obs, ] <- x1

    ## Make sure we haven't run too long.
    if(nevals > limit * sample.size) {
      X <- X[1:obs,,drop=FALSE]
      break
    }
  }

  return(X[nrow(X),])
}

## stepout.one.coord performs the stepping-out procedure on a single
## coordinate, returning an estimate of the slice endpoints.  The
## arguments are:
##
##   L           log.density function
##   x           current point
##   i           index of coordinate to estimate interval in
##   y.slice     slice level
##   w           initial slice size estimate
##   max.evals   maximum number of times to call L before aborting
stepout.one.coord <- function(L, g, i, y.slice, w, max.evals = NULL, ...) {
  ## Choose left and right to be the same as x, but with one coordinate
  ## changed in each of them so that together they bound a randomly
  ## positioned interval around g in coordinate i.
  left <- g
  left[i] <- left[i] - runif(1) * w
  right <- g
  right[i] <- left[i] + w

  ## Compute the log density at each end of the proposed interval.
  left.y <- L(left, ...)
  right.y <- L(right, ...)
  nevals <- 2

  ## Keep expanding as long as either end has log density smaller
  ## than the slice level.
  while(left.y > y.slice || right.y > y.slice) {
    ## Expand each side of the proposed slice with equal probability.
    if(runif(1) > 0.5) {
      right <- right + w
      right.y <- L(right, ...)
    } else {
      left <- left - w
      left.y <- L(left, ...)
    }
    nevals <- nevals + 1

    ## Make sure we haven't run too long.  If we have, return without
    ## an interval.
    if(!is.null(max.evals) && nevals >= max.evals)
      return(list(nevals = nevals))
  }

  ## We found an acceptable interval; return it.
  return(list(left = left[i], right = right[i], nevals = nevals))
}


propose_wslice <- function(x, family,
  response, eta, id, rho, ...)
{
  args <- list(...)
  iter <- args$iter  

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit

  if(iter %% x$xt$step == 0) {
    ## Map predictor to parameter scale.
    peta <- family$map2par(eta)

    ## Compute weights.
    weights <- family$weights[[id]](response, peta)

    ## Which obs to take.
    ok <- !(weights %in% c(NA, -Inf, Inf, 0))
    weights <- weights[ok]

    ## Score.
    score <- family$score[[id]](response, peta)

    ## Compute working observations.
    z <- eta[[id]][ok] + 1 / weights[ok] * score[ok]

    ## Compute mean and precision.
    XW <- t(x$X[ok, , drop = FALSE] * weights[ok])
    P <- if(x$fixed) {
      if(k <- ncol(x$X) < 2) {
        1 / (XW %*% x$X[ok, , drop = FALSE])
      } else matrix_inv(XW %*% x$X[ok, , drop = FALSE])
    } else {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / x$state$tau2[j] * x$S[[j]]
      matrix_inv(XW %*% x$X[ok, , drop = FALSE] + S)
    }
    P[P == Inf] <- 0
    x$state$g <- drop(P %*% (XW %*% (z - eta[[id]][ok])))
  }

  for(j in seq_along(x$state$g)) {
    x$state$g <- uni.slice(x$state$g, x, family, response, eta, id, j,
      logPost = logPost2, rho = rho)
  }

  ## Setup return state.
  x$state$alpha <- log(2)
  x$state$fit <- drop(x$X %*% x$state$g)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & is.null(args$no.mcmc)) {
    for(j in seq_along(x$state$tau2)) {
      x$state$tau2 <- uni.slice(x$state$tau2, x, family, response, eta, id, j,
        logPost = logPost4, rho = rho, lower = 0)
    }
  }

  return(x$state)
}


propose_oslice <- function(x, family,
  response, eta, id, rho, ...)
{
  args <- list(...)
  iter <- args$iter  

  if(iter %% x$xt$step == 0) {
    if(!is.null(x$state$mode)) {
      x$state$g <- x$state$mode$g
    } else {
      opt <- update_optim2(x, family, response, eta, id, ...)
      x$state$g <- opt$g
    }
  }

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit

  for(j in seq_along(x$state$g)) {
    x$state$g <- uni.slice(x$state$g, x, family, response, eta, id, j,
      logPost = logPost2, rho = rho)
  }

  ## Setup return state.
  x$state$alpha <- log(2)
  x$state$fit <- x$get.mu(x$X, x$state$g)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & is.null(args$no.mcmc)) {
    for(j in seq_along(x$state$tau2)) {
      x$state$tau2 <- uni.slice(x$state$tau2, x, family, response, eta, id, j,
        logPost = logPost4, rho = rho, lower = 0)
    }
  }

  return(x$state)
}


## Backfitting updating functions.
update_iwls <- function(x, family, response, eta, id, ...)
{
  args <- list(...)

  peta <- family$map2par(eta)

  if(is.null(args$weights)) {
    ## Compute weights.
    weights <- family$weights[[id]](response, peta)
  } else weights <- args$weights

  ## Which obs to take.
  ok <- is.finite(weights) & !is.na(weights) & (weights != 0) ##!(weights %in% c(NA, -Inf, Inf, 0))
  weights <- weights[ok]

  if(length(weights) == 0) return(x$state)

  if(is.null(args$z)) {
    ## Score.
    score <- family$score[[id]](response, peta)

    ## Compute working observations.
    z <- eta[[id]][ok] + 1 / weights * score[ok]
  } else z <- args$z[ok]

  ## Compute partial predictor.
  eta[[id]][ok] <- eta[[id]][ok] - x$state$fit[ok]

  ## Compute mean and precision.
  XW <- t(x$X[ok, ] * weights)
  XWX <- XW %*% x$X[ok, ]
  if(!x$optimize | x$fixed | !is.null(x$sp)) {
    if(x$fixed) {
      P <- matrix_inv(XWX)
    } else {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / x$state$tau2[j] * x$S[[j]]
      P <- matrix_inv(XWX + S)
    }
    x$state$g <- drop(P %*% (XW %*% (z - eta[[id]][ok])))
  } else {
    args <- list(...)
    edf0 <- args$edf - x$state$edf
    eta2 <- eta
    e <- z - eta[[id]][ok]

    objfun <- function(tau2, ...) {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      P <- matrix_inv(XWX + S)
      if(inherits(P, "try-error")) return(NA)
      g <- drop(P %*% (XW %*% e))
      if(any(is.na(g)) | any(g %in% c(-Inf, Inf))) g <- rep(0, length(g))
      fit <- drop(x$X %*% g)
      edf <- sum(diag(P %*% XWX))
      if(!is.null(x$xt$center)) {
        if(x$xt$center) edf <- edf - 1
      }
      eta2[[id]] <- eta2[[id]] + fit
      IC <- get.ic(family, response, family$map2par(eta2), edf0 + edf, length(e), x$criterion)
      return(IC)
    }

    if(length(x$state$tau2) < 2) {
      x$state$tau2 <- try(optimize(objfun, interval = x$interval, grid = x$grid)$minimum, silent = TRUE)
      if(inherits(x$state$tau2, "try-error"))
        x$state$tau2 <- optimize2(objfun, interval = x$interval, grid = x$grid)$minimum
      if(!length(x$state$tau2)) x$state$tau2 <- x$interval[1]
    } else {
      i <- grep("tau2", names(x$lower))
      opt <- try(optim(x$state$tau2, fn = objfun, method = "L-BFGS-B",
        lower = x$lower[i], upper = x$upper[i]), silent = TRUE)
      if(!inherits(opt, "try-error"))
        x$state$tau2 <- opt$par  
    }
    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / x$state$tau2[j] * x$S[[j]]
    P <- matrix_inv(XWX + S)
    x$state$g <- drop(P %*% (XW %*% e))
  }

  ## Compute fitted values.
  if(any(is.na(x$state$g)) | any(x$state$g %in% c(-Inf, Inf))) x$state$g <- rep(0, length(x$state$g))
  x$state$fit <- drop(x$X %*% x$state$g)
  x$state$edf <- sum(diag(P %*% XWX))
  if(!is.null(x$xt$center)) {
    if(x$xt$center) x$state$edf <- x$state$edf - 1
  }

  return(x$state)
}


## Update function using optim() including smoothing parameter selection.
update_optim <- function(x, family, response, eta, id, ...)
{
  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit
  eta2 <- eta

  if(is.null(x$state$XX)) x$state$XX <- crossprod(x$X)
  if(!x$fixed) {
    args <- list(...)
    edf0 <- args$edf - x$state$edf
    if(is.null(x$state$XX))
      x$state$XX <- crossprod(x$X)
  }

  if(!x$fixed & is.null(x[["sp"]])) {
    ## Objective function for variance parameter.
    objfun2 <- function(tau2) {
      ## Objective for regression coefficients.
      objfun <- function(gamma) {
        eta2[[id]] <- eta[[id]] + x$get.mu(x$X, gamma)
        ll <- family$loglik(response, family$map2par(eta2))
        lp <- x$prior(gamma, tau2)
        -1 * (ll + lp)
      }

      ## Gradient function.
      grad <- if(!is.null(family$score[[id]])) {
        function(gamma) {
          eta2[[id]] <- eta[[id]] + x$get.mu(x$X, gamma)
          peta <- family$map2par(eta2)
          score <- drop(family$score[[id]](response, peta))
          grad <- x$grad(score, gamma, tau2, full = FALSE)
          return(drop(-1 * grad))
        }
      } else NULL

      suppressWarnings(opt <- try(optim(x$state$g, fn = objfun, gr = grad,
        method = "BFGS", control = list()), silent = TRUE))

      if(!inherits(opt, "try-error")) {
        x$state$g <- opt$par
        x$state$fit <- x$get.mu(x$X, x$state$g)
      }

      edf <- x$edf(x, tau2)

      eta2[[id]] <- eta[[id]] + x$state$fit

      IC <- get.ic(family, response, family$map2par(eta2), edf0 + edf, length(eta2[[id]]), x$criterion)
      IC
    }

    if(length(x$state$tau2) < 2) {
      x$state$tau2 <- try(optimize(objfun2, interval = x$interval, grid = x$grid)$minimum, silent = TRUE)
      if(inherits(x$state$tau2, "try-error"))
        x$state$tau2 <- optimize2(objfun2, interval = x$interval, grid = x$grid)$minimum
      if(!length(x$state$tau2)) x$state$tau2 <- x$interval[1]
    } else {
      i <- grep("tau2", names(x$lower))
      opt <- try(optim(x$state$tau2, fn = objfun2, method = "L-BFGS-B",
        lower = x$lower[i], upper = x$upper[i]), silent = TRUE)
      if(!inherits(opt, "try-error"))
        x$state$tau2 <- opt$par  
    }
  }

  ## Objective for regression coefficients.
  objfun <- function(gamma) {
    eta2[[id]] <- eta[[id]] + x$get.mu(x$X, gamma)
    ll <- family$loglik(response, family$map2par(eta2))
    lp <- x$prior(gamma, x$state$tau2)
    -1 * (ll + lp)
  }

  ## Gradient function.
  grad <- if(!is.null(family$score[[id]])) {
    function(gamma) {
      eta2[[id]] <- eta[[id]] + x$get.mu(x$X, gamma)
      peta <- family$map2par(eta2)
      score <- drop(family$score[[id]](response, peta))
      grad <- x$grad(score, gamma, x$state$tau2, full = FALSE)
      return(drop(-1 * grad))
    }
  } else NULL

  a <- suppressWarnings(opt <- try(optim(x$state$g, fn = objfun, gr = grad,
    method = "BFGS", control = list()), silent = TRUE))

  if(!inherits(opt, "try-error")) {
    x$state$g <- opt$par
    x$state$fit <- x$get.mu(x$X, x$state$g)
  }

  if(!x$fixed)
    x$state$edf <- x$edf(x, x$state$tau2)

  return(x$state)
}


## Simple update function for regression coefficients.
update_optim2 <- function(x, family, response, eta, id, ...)
{
  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit
  eta2 <- eta

  ## Objective function.
  objfun <- function(gamma) {
    eta2[[id]] <- eta[[id]] + x$get.mu(x$X, gamma)
    peta <- family$map2par(eta2)
    ll <- family$loglik(response, peta)
    lp <- x$prior(gamma, x$state$tau2)
    -1 * (ll + lp)
  }

  grad <- if(!is.null(family$score[[id]]) & !is.null(x$grad)) {
    function(gamma) {
      eta2[[id]] <- eta[[id]] + x$get.mu(x$X, gamma)
      peta <- family$map2par(eta2)
      score <- drop(family$score[[id]](response, peta))
      grad <- x$grad(score, gamma, x$state$tau2, full = FALSE)
      return(drop(-1 * grad))
    }
  } else NULL

  suppressWarnings(opt <- try(optim(x$state$g, fn = objfun, method = "BFGS",
    control = list()), silent = TRUE))

  if(!inherits(opt, "try-error")) {
    x$state$g <- opt$par
    x$state$fit <- x$get.mu(x$X, x$state$g)
  }

  if(!x$fixed) {
    if(is.null(x$state$XX))
      x$state$XX <- crossprod(x$X)
    x$state$edf <- x$edf(x, x$state$tau2)
  }

  return(x$state)
}

## 3rd optimzer updater.
update_optim3 <- function(x, family, response, eta, id, ...)
{
  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit
  eta2 <- eta

  if(is.null(x$state$XX)) x$state$XX <- crossprod(x$X)
  args <- list(...)
  edf0 <- args$edf - x$state$edf
  if(!x$fixed) {
    if(is.null(x$state$XX))
      x$state$XX <- crossprod(x$X)
  }

  par <- c(x$state$g, if(!x$fixed) x$state$tau2 else NULL)
  names(par) <- x$s.colnames
  n <- length(eta[[id]])

  ## Objective function.
  objfun <- function(par) {
    gamma <- par[grep("c", x$s.colnames)]
    tau2 <- if(x$fixed) {
      NULL
    } else {
      if(!is.null(x$sp)) x$sp else par[grep("tau2", x$s.colnames)]
    }
    eta2[[id]] <- eta[[id]] + x$get.mu(x$X, gamma)
    ll <- family$loglik(response, family$map2par(eta2))
    edf <- edf0 + x$edf(x, tau2)
    val <- switch(x$criterion,
      "AIC" = -1 * (2 * ll + 2 * edf),
      "BIC" = -1 * (2 * ll + edf * log(n)),
      "AICc" = -1 * (2 * ll + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1))
    )
    return(val)
  }

  a <- suppressWarnings(opt <- try(optim(par, fn = objfun,
    method = "L-BFGS-B", lower = x$lower, upper = x$upper), silent = TRUE))

  if(!inherits(opt, "try-error")) {
    par <- opt$par
    x$state$g <- par[grep("c", x$s.colnames)]
    x$state$fit <- x$get.mu(x$X, x$state$g)
    x$state$tau2 <- if(!x$fixed) {
      if(!is.null(x$sp)) x$sp else par[grep("tau2", x$s.colnames)]
    }
  } else print(opt)

  x$state$edf <- x$edf(x, x$state$tau2)

  return(x$state)
}


## Setup for BayesG engine, handling
## sampling functions.
transformBayesG <- function(x, ...)
{
  call <- x$call; x$call <- NULL
  x <- assign.weights(x)

  family <- attr(x, "family")
  cat <- if(is.null(family$cat)) FALSE else family$cat

  if(cat) {
    if(length(x) != length(family$names)) {
      family$names <- paste(family$names[1], 1:length(x), sep = "")
      names(x) <- family$names
      family$links <- rep(family$links, length.out = length(x))
      names(family$links) <- family$names
      linkinv <- vector(mode = "list", length = length(family$names))
      for(j in family$names)
        linkinv[[j]] <- make.link2(family$links[j])$linkinv
      family$map2par <- function(eta) {
        for(j in names(eta)) {
          eta[[j]] <- linkinv[[j]](eta[[j]])
          eta[[j]][is.na(eta[[j]])] <- 0
          if(any(jj <- eta[[j]] == Inf))
            eta[[j]][jj] <- 10
          if(any(jj <- eta[[j]] == -Inf))
            eta[[j]][jj] <- -10
        }
        return(eta)
      }
      attr(x, "family") <- family
    }
  }

  tBayesG <- function(obj, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(obj)
      if(length(unique(nx)) < length(obj)) nx <- 1:length(obj)
      for(j in nx)
        obj[[j]] <- tBayesG(obj[[j]], ...)
    } else {
      if(!is.null(dim(obj$X))) {
        if(nrow(obj$X) > 0 & !is.na(mean(unlist(obj$X), na.rm = TRUE))) {
          if(is.null(obj$smooth)) obj$smooth <- list()
          obj$smooth[["parametric"]] <- list(
            "X" = obj$X,
            "S" = list(diag(0, ncol(obj$X))),
            "rank" = ncol(obj$X),
            "term" = "linear",
            "label" = "linear",
            "bs.dim" = ncol(obj$X),
            "fixed" = TRUE,
            "by" = "NA",
            "is.parametric" = TRUE
          )
          obj$sterms <- c(obj$strems, "parametric")
          obj$X <- NULL
        }
      }
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          obj$smooth[[j]] <- smooth.BayesG(obj$smooth[[j]])
          if(!is.null(obj$smooth[[j]]$rank))
            obj$smooth[[j]]$rank <- as.numeric(obj$smooth[[j]]$rank)
          if(!is.null(obj$smooth[[j]]$Xf)) {
            obj$smooth[[j]]$Xfcn <- paste(paste(paste(obj$smooth[[j]]$term, collapse = "."),
              "Xf", sep = "."), 1:ncol(obj$smooth[[j]]$Xf), sep = ".")
            colnames(obj$smooth[[j]]$Xf) <- obj$smooth[[j]]$Xfcn
            if(is.null(obj$smooth[["parametric"]])) {
              obj$smooth[["parametric"]] <- list(
                "X" = obj$smooth[[j]]$Xf,
                "S" = list(diag(0, ncol(obj$X))),
                "rank" = ncol(obj$smooth[[j]]$Xf),
                "term" = "linear",
                "label" = "linear",
                "bs.dim" = ncol(obj$smooth[[j]]$Xf),
                "fixed" = TRUE,
                "by" = "NA",
                "is.parametric" = TRUE
              )
              obj$sterms <- c(obj$strems, "parametric")
            } else {
              cn <- colnames(obj$smooth[["parametric"]]$X)
              obj$smooth[["parametric"]]$X <- cbind(obj$smooth[["parametric"]]$X, obj$smooth[[j]]$Xf)
              obj$smooth[["parametric"]]$S <- list(diag(0, ncol(obj$smooth[["parametric"]]$X)))
              obj$smooth[["parametric"]]$bs.dim <- list(diag(0, ncol(obj$smooth[["parametric"]]$X)))
              cn <- gsub("Intercept.", "Intercept", gsub("X.", "", c(cn , obj$smooth[[j]]$Xfcn), fixed = TRUE))
              obj$smooth[["parametric"]]$s.colnames <- colnames(obj$smooth[["parametric"]]$X) <- cn 
            }
          }
        }
      }
    }
    obj
  }

  x <- tBayesG(x, ...)

  attr(x, "call") <- call
  attr(x, "response.vec") <- attr(x, "model.frame")[, attr(attr(x, "model.frame"), "response.name")]

  if(cat) {
    response <- attr(x, "response.vec")
    if(!is.factor(response)) {
      if(!is.null(dim(response))) {
        ref <- apply(response, 1, function(x) { all(x == 0) * 1 })
        if(all(ref < 1)) stop("too many categories specified in formula!")
        response <- cbind(response, ref)
        response <- t(t(response) * 1:ncol(response))
        response <- drop(apply(response, 1, sum))
      }
    } else response <- as.integer(response)
    attr(x, "response.vec") <- response
  }

  x
}


optimize2 <- function(f, interval, ...,
  lower = min(interval), upper = max(interval), 
  maximum = FALSE, grid = 1000)
{
  f2 <- function(arg) f(arg, ...)
  gridvals <- seq(lower, upper, length = grid)
  fvals <- rep(NA, length = grid)
  fvals <- unlist(sapply(gridvals, f2))
  val <- gridvals[i <- if(maximum) which.max(fvals) else which.min(fvals)]
  list(minimum = val, objective = fvals[i])
}


## Function to setup BayesG smooths.
smooth.BayesG <- function(x, ...) {
  UseMethod("smooth.BayesG")
}

smooth.BayesG.default <- function(x, ...)
{
  if(!is.null(x$xbin.ind))
    x$X <- unique(x$X)
  x$a <- if(is.null(x$xt$a)) 1e-04 else x$xt$a
  x$b <- if(is.null(x$xt$b)) 1e-04 else x$xt$b
  if(is.null(x$fixed)) {
    x$fixed <- if(!is.null(x$fx)) x$fx[1] else FALSE
  }
  if(!x$fixed & is.null(x$interval)) {
    x$interval <- if(is.null(x$xt$interval)) tau2interval(x) else x$xt$interval
  }
  x$grid <- if(is.null(x$xt$grid)) 40 else x$xt$grid
  ntau2 <- length(x$S)
  if(length(ntau2) < 1)
    x$sp <- NULL
  if(!is.null(x$sp)) x$sp <- rep(x$sp, length.out = ntau2)
  if(is.null(x$state)) {
    x$p.save <- c("g", "tau2")
    x$state <- list()
    if(is.logical(x$np)) x$np <- NULL
    x$state$g <- rep(0, if(is.null(x$np)) ncol(x$X) else x$np)
    x$state$tau2 <- if(is.null(x$sp)) {
      if(x$fixed) 1e-20 else rep(if(!is.null(x$xt$tau2)) x$xt$tau2 else 10, length.out = ntau2)
    } else rep(x$sp, length.out = ntau2)
    x$s.colnames <- if(is.null(x$s.colnames)) {
      c(paste("c", 1:length(x$state$g), sep = ""),
        if(!x$fixed) paste("tau2", 1:length(x$state$tau2), sep = "") else NULL)
    } else x$s.colnames
    x$np <- c(length(x$state$g), length(x$state$tau2))
    x$state$XX <- crossprod(x$X)
  }
  if(is.null(x$xt$adaptive))
    x$xt$adaptive <- TRUE
  if(is.null(x$xt$step))
    x$xt$step <- 40
  if(is.null(x$get.mu) | !is.function(x$get.mu)) {
    x$get.mu <- function(X, b) {
      drop(as.matrix(X) %*% as.numeric(b))
    }
  }
  if(!is.null(x$xt$prior))
    x$prior <- x$xt$prior
  if(is.null(x$prior) | !is.function(x$prior)) {
    x$prior <- function(gamma, tau2 = NULL) {
      if(x$fixed | is.null(tau2)) {
        lp <- sum(dnorm(gamma, sd = 10, log = TRUE))
      } else {
        if(!is.null(x$sp)) tau2 <- x$sp
        lp <- 0
        for(j in seq_along(tau2)) {
          lp <- lp + -log(tau2[j]) * x$rank[j] / 2 + drop(-0.5 / tau2[j] * crossprod(gamma, x$S[[j]]) %*% gamma) +
            log((x$b^x$a)) - log(gamma(x$a)) + (-x$a - 1) * log(tau2[j]) - x$b / tau2[j]
        }
      }
      return(lp)
    }
  }
  if(is.null(x$edf) | !is.function(x$edf)) {
    x$edf <- function(x, tau2 = 0) {
      if(x$fixed) return(ncol(x$X))
      if(is.null(x$state$XX))
        x$state$XX <- crossprod(x$X)
      if(!is.null(x$sp)) tau2 <- x$sp
      S <- 0
      for(j in seq_along(tau2))
        S <- S + 1 / tau2[j] * x$S[[j]]
      P <- matrix_inv(x$state$XX + S)
      edf <- sum(diag(x$state$XX %*% P))
      if(!is.null(x$xt$center)) {
        if(x$xt$center) edf <- edf - 1
      }
      return(edf)
    }
  }
  if(is.null(x$grad) | !is.function(x$grad)) {
    x$grad <- function(score, gamma, tau2 = NULL, full = TRUE) {
      grad2 <- NULL
      if(x$fixed) {
        grad <- 0
      } else {
        grad <- 0; grad2 <- NULL
        for(j in seq_along(tau2)) {
          gS <- crossprod(gamma, x$S[[j]])
          grad <- grad + drop(-0.5 / tau2[j] * gS)
          if(full & !is.null(tau2[j])) {
            grad2 <- c(grad2, drop(-x$rank[j] / (2 * tau2[j]) - 1 / (2 * tau2[j]^2) * gS %*% gamma + (-x$a - 1) / tau2[j] + x$b / (tau2[j]^2)))
            x$X <- cbind(x$X, 0)
          }
        }
      }
      if(!is.null(x$xbin.ind))
        x$X <- x$X[x$xbin.ind, , drop = FALSE]
      grad <- drop(crossprod(x$X, score)) + c(grad, grad2)
      return(grad)
    }
  } else {
    if(!is.function(x$grad))
      x$grad <- NULL
  }

  x$state$edf <- x$edf(x, x$state$tau2)

  x$lower <- c(rep(-Inf, length(x$state$g)),
    if(is.list(x$interval)) unlist(sapply(x$interval, function(x) { x[1] })) else x$interval[1])
  x$upper <- c(rep(Inf, length(x$state$g)),
    if(is.list(x$interval)) unlist(sapply(x$interval, function(x) { x[2] })) else x$interval[2])
  names(x$lower) <- names(x$upper) <- x$s.colnames
  if(!is.null(x$sp)) {
    if(length(x$sp) < 1)
      x$sp <- NULL
    if(is.logical(x$sp))
      x[["sp"]] <- NULL
  }

  x
}


get.ic <- function(family, response, eta, edf, n, type = c("AIC", "BIC", "AICc", "MP"))
{
  type <- match.arg(type)
  ll <- family$loglik(response, eta)
  pen <- switch(type,
    "AIC" = -2 * ll + 2 * edf,
    "BIC" = -2 * ll + edf * log(n),
    "AICc" = -2 * ll + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
    "MP" = -1 * (ll + edf)
  )
  return(pen)
}


## The generic Bayes engine.
BayesG <- function(x, n.iter = 12000, thin = 10, burnin = 2000, accept.only = TRUE,
  verbose = TRUE, step = 20, svalues = TRUE, eps = .Machine$double.eps^0.25, maxit = 400,
  n.adapt = floor(0.1 * burnin), tdir = NULL, method = "MP", outer = FALSE, inner = FALSE, n.samples = 200,
  criterion = c("AICc", "BIC", "AIC"), lower = 1e-09, upper = 1e+04,
  optim.control = NULL, digits = 3,  ## list(pgtol = 1e-04, maxit = 5)
  update = c("optim", "iwls", "optim2", "optim3"),
  propose = c("mvn", "iwls", "slice", "slice2", "oslice", "wslice", "iwls0"),
  sample = c("slice", "slice2", "iwls", "iwls0"), optimize = TRUE, ...)
{
  known_methods <- c("backfitting", "MCMC", "backfitting2",
    "backfitting3", "backfitting4", "mcmc", "MP", "mp", "LD", "ld", "mp2", "MP2")
  if(is.integer(method))
    method <- known_methods[method]
  tm <- NULL
  for(m in method)
    tm <- c(tm, match.arg(m, known_methods))
  method <- tm

  family <- attr(x, "family")
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  response <- attr(x, "response.vec")
  criterion <- match.arg(criterion)
#  scipen <- getOption("scipen")
#  options("scipen" = 20)
#  on.exit(options("scipen" = scipen))

  ## Actual number of samples to save.
  if(!any(grepl("MCMC", method)) & !any(grepl("LD", method))) {
    if(nosamp <- n.samples < 1)
      n.samples <- 1
    n.iter <- n.samples
    thin = 1
    burnin = 0
  }
  if(n.iter < burnin) stop("argument burnin exceeds n.iter!")
  if(thin > (n.iter - burnin)) stop("argument thin is set too large!")
  iterthin <- as.integer(seq(burnin, n.iter, by = thin))
  n.save <- length(iterthin)
  nstep <- step
  step <- floor(n.iter / step)

  getDigits <- function(x) {
    nchar(gsub("(.*\\.)|([0]*$)", "", as.character(x)))
  }
  epsdigits <- getDigits(eps)

  ## The update function to be used within backfitting.
  if(!is.function(update)) {
    update <- match.arg(update)
    update <- switch(update,
      "optim" = update_optim,
      "iwls" = update_iwls,
      "optim2" = update_optim2,
      "optim3" = update_optim3
    )
  }

  ## The proposal function to be used for smooth terms.
  if(!is.function(propose)) {
    propose <- match.arg(propose)
    propose <- switch(propose,
      "iwls" = propose_iwls,
      "iwls0" = propose_iwls0,
      "slice" = propose_slice,
      "slice2" = propose_slice2,
      "oslice" = propose_oslice,
      "wslice" = propose_wslice,
      "mvn" = propose_mvn
    )
  }

  ## Function for creating samples when using backfitting.
  if(!is.function(sample)) {
    sample <- match.arg(sample)
    sample <- switch(sample,
      "iwls" = propose_iwls,
      "iwls0" = propose_iwls0,
      "slice" = propose_slice,
      "slice2" = propose_slice2,
      "oslice" = propose_oslice,
      "wslice" = propose_wslice,
      "mvn" = propose_mvn
    )
  }
  
  ## Add acceptance rate and fitted values vectors.
  smBayesG <- function(obj, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(obj)
      if(length(unique(nx)) < length(obj)) nx <- 1:length(obj)
      for(j in nx)
        obj[[j]] <- smBayesG(obj[[j]], ...)
    } else {
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          if(!is.null(obj$smooth[[j]]$is.parametric)) {
            obj$smooth[[j]]$np <- ncol(obj$smooth[[j]]$X)
            obj$smooth[[j]]$p.save = "g"
            obj$smooth[[j]]$s.colnames <- paste("c", 1:ncol(obj$smooth[[j]]$X), sep = "")
          }
          obj$smooth[[j]]$s.alpha <- rep(0, nrow = n.save)
          obj$smooth[[j]]$s.accepted <- rep(0, nrow = n.save)
          obj$smooth[[j]]$s.samples <- matrix(0, nrow = n.save, ncol = sum(obj$smooth[[j]]$np))
          obj$smooth[[j]]$state$fit <- rep(0, nrow(obj$smooth[[j]]$X))
          obj$smooth[[j]]$state$g <- rep(0, obj$smooth[[j]]$np[1])
          obj$smooth[[j]]$fxsp <- if(!is.null(obj$smooth[[j]]$sp)) TRUE else FALSE
          if(!is.null(update) & is.null(obj$smooth[[j]]$update))
            obj$smooth[[j]]$update <- update
          if(!is.null(propose) & is.null(obj$smooth[[j]]$propose))
            obj$smooth[[j]]$propose <- propose
          if(!is.null(sample) & is.null(obj$smooth[[j]]$sample))
            obj$smooth[[j]]$sample <- sample
          obj$smooth[[j]]$state$accepted <- FALSE
          obj$smooth[[j]]$state$scale <- rep(1, length = obj$smooth[[j]]$np[1])
          obj$smooth[[j]]$state$iter <- 1
          obj$smooth[[j]]$state$maxit <- 50
          obj$smooth[[j]]$adapt <- n.adapt
          obj$smooth[[j]]$optimize <- optimize
          obj$smooth[[j]]$criterion <- criterion
        }
      }
    }
    obj
  }

  x <- smBayesG(x, ...)

  ## Formatting for printing.
  fmt <- function(x, width = 8, digits = 2) {
    txt <- formatC(round(x, digits), format = "f", digits = digits , width = width)
    if(nchar(txt) > width) {
      txt <- strsplit(txt, "")[[1]]
      txt <- paste(txt[1:width], collapse = "", sep = "")
    }
    txt
  }

  ## Number of parameters
  np <- length(nx)

  ## Set up predictors.
  eta <- vector(mode = "list", length = np)
  names(eta) <- nx
  for(j in 1:np)
    eta[[j]] <- rep(0, length(response))

  ## Extract edf or logpriors.
  get_edf_lp <- function(x, logprior = FALSE) {
    rval <- 0
    for(j in 1:np) {
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        if(logprior) {
          selp <- x[[nx[j]]]$smooth[[sj]]$prior(
            x[[nx[j]]]$smooth[[sj]]$state$g,
            x[[nx[j]]]$smooth[[sj]]$state$tau2)
        } else {
          selp <- x[[nx[j]]]$smooth[[sj]]$edf(x[[nx[j]]]$smooth[[sj]],
            x[[nx[j]]]$smooth[[sj]]$state$tau2)
        }
        rval <- rval + selp
      }
    }
    rval
  }

  fxMP <- "mp2" %in% tolower(method)
  nobs <- if(is.null(dim(response))) length(response) else nrow(response)

  ## Function to create full parameter vector.
  make_par <- function(x) {
    par <- npar <- lower <- upper <- NULL
    grad <- TRUE
    for(j in 1:np) {
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        par <- c(par,
          x[[nx[j]]]$smooth[[sj]]$state$g,
          if(!x[[nx[j]]]$smooth[[sj]]$fixed & !fxMP) x[[nx[j]]]$smooth[[sj]]$state$tau2 else NULL
        )
        if(is.null(x[[nx[j]]]$smooth[[sj]]$grad)) grad <- FALSE
        ng <- length(x[[nx[j]]]$smooth[[sj]]$state$g)
        nv <- if(fxMP) NULL else length(x[[nx[j]]]$smooth[[sj]]$state$tau2)
        lower <- c(lower, x[[nx[j]]]$smooth[[sj]]$lower)
        upper <- c(upper, x[[nx[j]]]$smooth[[sj]]$upper)
        npar <- c(npar, paste("p", j, "t", sj,
          c(paste("c", 1:ng, sep = ""),
          if(!x[[nx[j]]]$smooth[[sj]]$fixed & !fxMP) paste("v", 1:nv, sep = "") else NULL),
          sep = ""))
      }
    }
    names(par) <- npar
    return(list("par" = par, "lower" = lower, "upper" = upper, "grad" = grad))
  }

  log_posterior_BayesG <- function(par, par2 = NULL, rx = FALSE, type = 2) {
    lprior <- edf2 <- 0
    for(j in 1:np) {
      eta[[nx[j]]] <- 0
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        gamma <- par[grep(paste("p", j, "t", sj, "c", sep = ""), names(par))]
        x[[nx[j]]]$smooth[[sj]]$state$g <- gamma
        if(!x[[nx[j]]]$smooth[[sj]]$fixed) {
          if(!fxMP) {
            tau2 <- par[grep(paste("p", j, "t", sj, "v", sep = ""), names(par))]
            x[[nx[j]]]$smooth[[sj]]$state$tau2 <- tau2
          } else {
            if(!is.null(par2)) {
              id <- paste("id", j, sj, sep = "")
              tau2 <- par2[grep(id, names(par2), fixed = TRUE)]
            } else {
              tau2 <- x[[nx[j]]]$smooth[[sj]]$state$tau2
            }
          }
        } else tau2 <- NULL
        if(rx) {
          x[[nx[j]]]$smooth[[sj]]$state$mode <- list(
            "g" = gamma,
            "tau2" = x[[nx[j]]]$smooth[[sj]]$state$tau2
          )
        }
        x[[nx[j]]]$smooth[[sj]]$state$fit <- x[[nx[j]]]$smooth[[sj]]$get.mu(x[[nx[j]]]$smooth[[sj]]$X, gamma)
        if(!is.null(x[[nx[j]]]$smooth[[sj]]$xbin.ind))
          x[[nx[j]]]$smooth[[sj]]$state$fit <- x[[nx[j]]]$smooth[[sj]]$state$fit[x[[nx[j]]]$smooth[[sj]]$xbin.ind]
        eta[[nx[j]]] <- eta[[nx[j]]] + x[[nx[j]]]$smooth[[sj]]$state$fit
        x[[nx[j]]]$smooth[[sj]]$state$edf <- x[[nx[j]]]$smooth[[sj]]$edf(x[[nx[j]]]$smooth[[sj]], tau2)
        lprior <- lprior + x[[nx[j]]]$smooth[[sj]]$prior(gamma, tau2)
        edf2 <- edf2 + x[[nx[j]]]$smooth[[sj]]$state$edf
      }
    }

    ll <- family$loglik(response, family$map2par(eta))
    if(rx) {
      rval <- list("x" = x, "eta" = eta, "lp" = ll + lprior, "ll" = ll)
    } else {
      edf <- get_edf_lp(x, logprior = FALSE)
      ic <- get.ic(family, response, family$map2par(eta), edf, nobs, criterion)
      rval <- ll + lprior
      if(type > 1)
        rval <- -1 * rval
      if(!is.finite(rval))
        rval <- NA
    }

    if(verbose & !rx & type > 1) {
      cat("\r")
      cat(criterion, fmt(ic, width = 8, digits = digits),
        "logPost", fmt(ll + lprior, width = 8, digits = digits),
        "logLik", fmt(ll, width = 8, digits = digits),
        "edf", fmt(edf2, width = 6, digits = digits))
    }

    return(rval)
  }

  hessian <- NULL
  ## Posterior mode estimation.
  if(any(grepl("mp", method, ignore.case = TRUE))) {
    tpar <- make_par(x = x)
    par <- tpar$par; lower2 <- tpar$lower; upper2 <- tpar$upper
    if(!is.null(family$score) & tpar$grad) {
      grad_posterior_BayesG <- function(par, par2 = NULL, type = 1, ...) {
        grad <- NULL
        for(j in 1:np) {
          eta[[nx[j]]] <- 0
          for(sj in seq_along(x[[nx[j]]]$smooth)) {
            gamma <- par[grep(paste("p", j, "t", sj, "c", sep = ""), names(par))]
            f <- x[[nx[j]]]$smooth[[sj]]$get.mu(x[[nx[j]]]$smooth[[sj]]$X, gamma)
            if(!is.null(x[[nx[j]]]$smooth[[sj]]$xbin.ind))
              f <- f[x[[nx[j]]]$smooth[[sj]]$xbin.ind]
            eta[[nx[j]]] <- eta[[nx[j]]] + f
          }
        }
        for(j in 1:np) {
          score <- family$score[[nx[j]]](response, family$map2par(eta))
          for(sj in seq_along(x[[nx[j]]]$smooth)) {
            gamma <- par[grep(paste("p", j, "t", sj, "c", sep = ""), names(par))]
            if(!x[[nx[j]]]$smooth[[sj]]$fixed) {
              if(!fxMP) {
                tau2 <- par[grep(paste("p", j, "t", sj, "v", sep = ""), names(par))]
                x[[nx[j]]]$smooth[[sj]]$state$tau2 <- tau2
              } else {
                if(!is.null(par2)) {
                  id <- paste("id", j, sj, sep = "")
                  tau2 <- par2[grep(id, names(par2), fixed = TRUE)]
                } else {
                  tau2 <- x[[nx[j]]]$smooth[[sj]]$state$tau2
                }
              }
            } else tau2 <- NULL
            tgrad <- x[[nx[j]]]$smooth[[sj]]$grad(score, gamma, tau2)
            if(fxMP) tgrad <- tgrad[1:length(gamma)]
            grad <- c(grad, tgrad)
          }
        }
        if(type < 2) grad <- grad * -1
        return(grad)
      }
    } else grad_posterior_BayesG <- NULL

    opt <- optim(par, fn = log_posterior_BayesG, gr = NULL, ##grad_posterior_BayesG,
      method = "L-BFGS-B", lower = lower2, upper = upper2,
      control = optim.control)
    par <- opt$par
    hessian <- opt$hessian
    rm(opt)
    opt <- log_posterior_BayesG(par, rx = TRUE)

    x <- opt$x; eta <- opt$eta
    rm(opt)

    if(verbose) cat("\n")
    svalues <- FALSE
  }

  ## Find starting values with backfitting.
  if(svalues | any(grepl("backfitting", method))) {
    inner_bf <- function(x, response, eta, family, edf, id, ...) {
      eps0 <- eps + 1; iter <- 1
      while(eps0 > eps & iter < maxit) {
        eta0 <- eta
        for(sj in seq_along(x)) {
          ## Get updated parameters.
          p.state <- x[[sj]]$update(x[[sj]], family, response, eta, id, edf = edf, ...)

          ## Compute equivalent degrees of freedom.
          edf <- edf - x[[sj]]$state$edf + p.state$edf

          ## Update predictor and smooth fit.
          eta[[id]] <- eta[[id]] - x[[sj]]$state$fit + p.state$fit
          x[[sj]]$state <- p.state
          x[[sj]]$state$mode <- list("g" = p.state$g, "tau2" = p.state$tau2)
        }
        eps0 <- do.call("cbind", eta)
        eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
        if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1
        iter <- iter + 1
      }
      return(list("x" = x, "eta" = eta, "edf" = edf))
    }

    ## 1st backfitting main function.
    backfit <- function(x, eta, verbose = TRUE) {
      edf <- get_edf_lp(x, logprior = FALSE)

      eps0 <- eps + 1; iter <- 1
      while(eps0 > eps & iter < maxit) {
        eta0 <- eta
        ## Cycle through all parameters
        for(j in 1:np) {
          if(outer) {
            peta <- family$map2par(eta)

            ## Compute weights.
            weights <- family$weights[[nx[j]]](response, peta)

            ## Score.
            score <- family$score[[nx[j]]](response, peta)

            ## Compute working observations.
            z <- eta[[nx[j]]] + 1 / weights * score
          } else z <- weights <- NULL

          ## And all terms.
          if(inner) {
            tbf <- inner_bf(x[[nx[j]]]$smooth, response, eta, family,
              edf = edf, id = nx[j], z = z, weights = weights)
            x[[nx[j]]]$smooth <- tbf$x
            edf <- tbf$edf
            eta <- tbf$eta
            rm(tbf)
          } else {
            for(sj in seq_along(x[[nx[j]]]$smooth)) {
              ## Get updated parameters.
              p.state <- x[[nx[j]]]$smooth[[sj]]$update(x[[nx[j]]]$smooth[[sj]],
                family, response, eta, nx[j], edf = edf, z = z, weights = weights)

              ## Compute equivalent degrees of freedom.
              edf <- edf - x[[nx[j]]]$smooth[[sj]]$state$edf + p.state$edf

              ## Update predictor and smooth fit.
              eta[[nx[j]]] <- eta[[nx[j]]] - x[[nx[j]]]$smooth[[sj]]$state$fit + p.state$fit

              x[[nx[j]]]$smooth[[sj]]$state <- p.state
              x[[nx[j]]]$smooth[[sj]]$state$mode <- list("g" = p.state$g, "tau2" = p.state$tau2)
            }
          }
        }

        eps0 <- do.call("cbind", eta)
        eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
        if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

        peta <- family$map2par(eta)

        if(any(method == "backfitting") | verbose) {
          IC <- get.ic(family, response, peta, edf, nobs, criterion)
          edf2 <- get_edf_lp(x, logprior = FALSE)
          cat("\r")
          vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
            " logLik ", fmt(family$loglik(response, peta), width = 8, digits = digits),
            " edf ", fmt(edf2, width = 6, digits = digits),
            " eps ", fmt(eps0, width = 6, digits = digits + 2),
            " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
          cat(vtxt)

          if(.Platform$OS.type != "unix") flush.console()
        }

        iter <- iter + 1
      }

      if(criterion == "MP")
        edf <- get_edf_lp(x, logprior = TRUE)
      IC <- get.ic(family, response, peta, edf, nobs, criterion)
      edf2 <- get_edf_lp(x, logprior = FALSE)

      if(any(method %in% c("backfitting", "backfitting2", "backfitting4")) | verbose) {
        cat("\r")
        vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
          " logLik ", fmt(family$loglik(response, peta), width = 8, digits = digits),
          " edf ", fmt(edf2, width = 6, digits = digits),
          " eps ", fmt(eps0, width = 6, digits = digits + 2),
          " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
        cat(vtxt)
        if(.Platform$OS.type != "unix") flush.console()
        if(any(method %in% "backfitting"))cat("\n")
      }

      if(iter == maxit)
        warning("the backfitting algorithm did not converge, please check argument eps and maxit!")

      return(list("x" = x, "eta" = eta, "ic" = IC))
    }

    bf <- backfit(x, eta, verbose = verbose)
    x <- bf$x; eta <- bf$eta
    rm(bf)

    if(any("backfitting2" %in% method)) {
      tau2 <- NULL
      for(j in 1:np) {
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          tmp <- x[[nx[j]]]$smooth[[sj]]$state$tau2
          names(tmp) <- paste("id", j, sj, ".", 1:length(tmp), sep = "")
          tau2 <- c(tau2, tmp)
        }
      }

      objfun <- function(tau2, retbf = FALSE) {
        for(j in 1:np) {
          for(sj in seq_along(x[[nx[j]]]$smooth)) {
            id <- paste("id", j, sj, sep = "")
            tmp <- tau2[grep(id, names(tau2), fixed = TRUE)]
            names(tmp) <- names(x[[nx[j]]]$smooth[[sj]]$state$tau2)
            x[[nx[j]]]$smooth[[sj]]$state$tau2 <- tmp
          }
        }

        bf <- backfit(x, eta, verbose = verbose)
        return(if(retbf) bf else bf$ic)
      }

      lower <- rep(lower, length.out = length(tau2))
      upper <- rep(upper, length.out = length(tau2))

      ## Use optim() to find variance parameters.
      tau2 <- optim(tau2, fn = objfun,
        method = "L-BFGS-B", lower = lower, upper = upper,
        control = optim.control)$par

      bf <- objfun(tau2, retbf = TRUE)
      eta <- bf$eta
      x <- bf$x
      rm(bf)
    }

    if(any("backfitting4" %in% method)) {
      objfun <- function(tau2, j, sj, retbf = FALSE) {
        if(length(tau2) > 1) stop('only single variances can be updated using "backftting4"!')
        x[[nx[j]]]$smooth[[sj]]$state$tau2 <- tau2
        bf <- backfit(x, eta, verbose = verbose)
        return(if(retbf) bf else bf$ic)
      }
      eps0 <- eps + 1; iter <- 1
      while(eps0 > eps & iter < 4) {
        eta0 <- eta
        for(j in 1:np) {
          for(sj in seq_along(x[[nx[j]]]$smooth)) {
            if(!x[[nx[j]]]$smooth[[sj]]$fixed) {
              tau2 <- optimize(objfun, x[[nx[j]]]$smooth[[sj]]$interval, j = j, sj = sj)$minimum
              bf <- objfun(tau2, j, sj, retbf = TRUE)
            } else {
              bf <- backfit(x, eta, verbose = verbose)
            }
            eta <- bf$eta
            x <- bf$x
            rm(bf)
          }
        }
        eps0 <- do.call("cbind", eta)
        eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
        iter <- iter + 1
      }
    }
  }

  if(verbose & any(method %in% c("backfitting2", "backfitting4")))
    cat("\n")

  deviance <- rep(0, length(iterthin))
  rho <- new.env()

  barfun <- function(ptm, n.iter, i, step, nstep, start = TRUE) {
    if(i == 10 & start) {
      elapsed <- c(proc.time() - ptm)[3]
      rt <- elapsed / i * (n.iter - i)
      rt <- if(rt > 60) {
        paste(formatC(format(round(rt / 60, 2), nsmall = 2), width = 5), "min", sep = "")
      } else paste(formatC(format(round(rt, 2), nsmall = 2), width = 5), "sec", sep = "")
      cat("|", rep(" ", nstep), "|   0% ", rt, sep = "")
      if(.Platform$OS.type != "unix") flush.console()
    }
    if(i %% step == 0) {
      cat("\r")
      p <- i / n.iter
      p <- paste("|", paste(rep("*", round(nstep * p)), collapse = ""),
        paste(rep(" ", round(nstep * (1 - p))), collapse = ""), "| ",
        formatC(round(p, 2) * 100, width = 3), "%", sep = "")
      elapsed <- c(proc.time() - ptm)[3]
      rt <- elapsed / i * (n.iter - i)
      rt <- if(rt > 60) {
        paste(formatC(format(round(rt / 60, 2), nsmall = 2), width = 5), "min", sep = "")
      } else paste(formatC(format(round(rt, 2), nsmall = 2), width = 5), "sec", sep = "")
      cat(p, rt, sep = " ")
      if(.Platform$OS.type != "unix") flush.console()
    }
  }

  if(any(grep("MCMC", method))) {
    ## Start sampling
    ptm <- proc.time()
    for(i in 1:n.iter) {
      if(save <- i %in% iterthin)
        js <- which(iterthin == i)
    
      ## Cycle through all parameters
      for(j in 1:np) {
        ## And all terms.
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          ## Get proposed states.
          p.state <- x[[nx[j]]]$smooth[[sj]]$propose(x[[nx[j]]]$smooth[[sj]], family, response, eta, nx[j], rho = rho, iter = i)

          ## If accepted, set current state to proposed state.
          accepted <- if(is.na(p.state$alpha)) FALSE else log(runif(1)) <= p.state$alpha

          if(accepted) {
            eta[[nx[j]]] <- eta[[nx[j]]] - x[[nx[j]]]$smooth[[sj]]$state$fit + p.state$fit
            x[[nx[j]]]$smooth[[sj]]$state <- p.state 
          }
          x[[nx[j]]]$smooth[[sj]]$state$accepted <- accepted
          x[[nx[j]]]$smooth[[sj]]$state$iter <- i

          ## Save the samples and acceptance.
          if(save) {
            x[[nx[j]]]$smooth[[sj]]$s.alpha[js] <- min(c(exp(p.state$alpha), 1), na.rm = TRUE)
            x[[nx[j]]]$smooth[[sj]]$s.accepted[js] <- accepted
            x[[nx[j]]]$smooth[[sj]]$s.samples[js, ] <- unlist(x[[nx[j]]]$smooth[[sj]]$state[x[[nx[j]]]$smooth[[sj]]$p.save])
          }

          ## Check.
#          if(any(abs(eta[[nx[j]]]) > 10)) {
#            eta[[nx[j]]][eta[[nx[j]]] > 10] <- 10
#            eta[[nx[j]]][eta[[nx[j]]] < -10] <- -10
#          }
        }
      }

      if(save) deviance[js] <- -2 * family$loglik(response, family$map2par(eta))

      if(verbose) barfun(ptm, n.iter, i, step, nstep)
    }
    if(verbose) cat("\n")
  }

  save.edf <- save.loglik <- NULL

  if(any(grep("LD", toupper(method)))) {
    stopifnot(require("LaplacesDemon"))

    method <- c(method, "MCMC")
    
    Model <- function(parm, Data) {
      lprior <- 0
      for(j in 1:np) {
        eta[[nx[j]]] <- 0
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          gamma <- parm[grep(paste("p", j, "t", sj, "c", sep = ""), Data$parm.names)]
          x[[nx[j]]]$smooth[[sj]]$state$g <- gamma
          if(!x[[nx[j]]]$smooth[[sj]]$fixed) {
            tau2 <- parm[grep(paste("p", j, "t", sj, "v", sep = ""), Data$parm.names)]
            x[[nx[j]]]$smooth[[sj]]$state$tau2 <- tau2
          } else tau2 <- NULL
          x[[nx[j]]]$smooth[[sj]]$state$fit <- x[[nx[j]]]$smooth[[sj]]$get.mu(x[[nx[j]]]$smooth[[sj]]$X, gamma)
          eta[[nx[j]]] <- eta[[nx[j]]] + x[[nx[j]]]$smooth[[sj]]$state$fit
          lprior <- lprior + x[[nx[j]]]$smooth[[sj]]$prior(gamma, tau2)
        }
      }

      ll <- drop(family$loglik(response, family$map2par(eta)))
      lp <- drop(ll + lprior)

      rval <- list("LP" = lp, "Dev" = -2 * ll, "Monitor" = lp, "yhat" = NA, "parm" = parm)

      return(rval)
    }

    par <- make_par(x = x)$par
    MyData <- list("mon.names" = "LP", "parm.names" = names(par), "N" = length(eta[[1]]))
    Initial.Values <- par

    Fit <- LaplacesDemon(Model, Data = MyData, Initial.Values,
      Iterations = n.iter, Thinning = 1, ...)
    samples <- Fit$Posterior1[iterthin, , drop = FALSE]
    deviance <- Fit$Deviance[iterthin]

    rm(Fit)

    for(j in 1:np) {
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        gamma <- samples[, grep(paste("p", j, "t", sj, "c", sep = ""), colnames(samples)), drop = FALSE]
        if(!x[[nx[j]]]$smooth[[sj]]$fixed) {
          tau2 <- samples[, grep(paste("p", j, "t", sj, "v", sep = ""), colnames(samples))]
        } else tau2 <- NULL
        x[[nx[j]]]$smooth[[sj]]$s.samples <- cbind(gamma, tau2)
      }
    }

    save.edf <- get_edf_lp(x)
    save.loglik <- family$loglik(response, family$map2par(eta))

    rm(samples)
  }


  ## Remove some settings.
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      x[[nx[j]]]$smooth[[sj]]$state$iter <- NULL
      x[[nx[j]]]$smooth[[sj]]$state$maxit <- 1
    }
  }

  if(!any(grepl("MCMC", method)) & TRUE) {   
    save.edf <- get_edf_lp(x)
    save.loglik <- family$loglik(response, family$map2par(eta))
    llim <- -Inf
    if(verbose) cat("generating samples\n")
    ptm <- proc.time()
    for(js in seq_along(iterthin)) {
      for(j in 1:np) {
        ## And all terms.
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          ## Get proposed states.
          p.state <- if(nosamp) {
            x[[nx[j]]]$smooth[[sj]]$state$alpha <- 1
            x[[nx[j]]]$smooth[[sj]]$state
          } else {
            x[[nx[j]]]$smooth[[sj]]$sample(x[[nx[j]]]$smooth[[sj]], family, response, eta, nx[j], rho = rho, no.mcmc = TRUE, llim = llim)
          }

          accepted <- if(is.na(p.state$alpha)) FALSE else log(runif(1)) <= p.state$alpha

          ## Save the samples.
          x[[nx[j]]]$smooth[[sj]]$s.alpha[js] <- min(c(exp(p.state$alpha), 1), na.rm = TRUE)
          x[[nx[j]]]$smooth[[sj]]$s.accepted[js] <- accepted
          x[[nx[j]]]$smooth[[sj]]$s.samples[js, ] <- unlist(p.state[x[[nx[j]]]$smooth[[sj]]$p.save])
          llim <- p.state$llim
        }
      }
      if(verbose) barfun(ptm, length(iterthin), js, 1, 20, start = FALSE)
    }
    if(verbose) cat("\n")
    deviance <- rep(-2 * save.loglik, length = length(iterthin))
  }

  if(accept.only) {
    for(j in 1:np) {
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        accepted <- x[[nx[j]]]$smooth[[sj]]$s.accepted
        accepted <- if(!any(accepted > 0)) 1 else accepted > 0
        x[[nx[j]]]$smooth[[sj]]$s.samples[!accepted, ]  <- NA
      }
    }
  }

  ## Return all samples as mcmc matrix.
  ## (1) Write out all samples to tdir.
  if(is.null(tdir)) {
    dir.create(tdir <- tempfile())
    on.exit(unlink(tdir))
  }
  tdir <- path.expand(tdir)

  samplesBayesG <- function(obj, id = NULL, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- paste("p", 1:length(obj), sep = "")
      if(length(unique(nx)) < length(obj)) nx <- paste("p", 1:length(obj), sep = "")
      for(j in nx)
        samplesBayesG(obj[[j]], id = j, ...)
    } else {
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          slab <- gsub("/", "RSdivRS", obj$smooth[[j]]$label, fixed = TRUE)
          fn <- file.path(tdir, paste(id, if(!is.null(id)) ":", "h",
            obj$hlevel, ":", slab, ".raw", sep = ""))
          obj$smooth[[j]]$s.samples <- cbind(obj$smooth[[j]]$s.samples, obj$smooth[[j]]$s.alpha)
          colnames(obj$smooth[[j]]$s.samples) <- c(obj$smooth[[j]]$s.colnames, "alpha")
          write.table(obj$smooth[[j]]$s.samples, file = fn, row.names = FALSE, quote = FALSE)
        }
      }
    }
    NULL
  }

  samplesBayesG(x)

  ## (2) Remove old object.
  rm(x)

  ## (3) Read back in samples
  sf <- grep(".raw", dir(tdir), value = TRUE)
  samples <- NULL
  for(j in sf) {
    st <- as.matrix(read.table(file.path(tdir, j), header = TRUE, colClasses = "numeric"))
    colnames(st) <- paste(gsub(".raw", "", j), gsub("c", "", colnames(st)), sep = ".")
    colnames(st) <- gsub("RSdivRS", "/", colnames(st))
    samples <- cbind(samples, st)
  }
  samples <- if(!any(grepl("MCMC", method))) {
    cbind(samples, "log Lik." = save.loglik, "save.edf" = save.edf)
  } else cbind(samples, "deviance" = deviance)

  return(as.mcmc(samples))
}


resultsBayesG <- function(x, samples)
{
  family <- attr(x, "family")
  grid <- attr(x, "grid")
  if(is.null(grid)) grid <- 100
  if(is.function(family))
    family <- family()

  createBayesGresults <- function(obj, samples, id = NULL)
  {
    if(inherits(samples[[1]], "mcmc.list")) {
      samples <- do.call("c", samples)
    } else {
      if(!inherits(samples, "mcmc.list"))
        samples <- as.mcmc.list(if(!inherits(samples, "list")) list(samples) else samples)
    }
    chains <- length(samples)
    rval <- vector(mode = "list", length = chains)
    snames <- colnames(samples[[1]])

    for(j in 1:chains) {
      DIC <- pd <- NA
      if(any(grepl("deviance", snames))) {
        DIC <- as.numeric(samples[[j]][, grepl("deviance", snames)])
        pd <- var(DIC, na.rm = TRUE) / 2
        DIC <- mean(DIC, na.rm = TRUE)
      }
      if(any(grepl("logLik", snames))) {
        DIC <- -2 * as.numeric(samples[[j]][, grepl("logLik", snames)])
        pd <- var(DIC, na.rm = TRUE) / 2
        DIC <- mean(DIC, na.rm = TRUE)
      }

      if(any(grepl("save.edf", snames))) {
        ic <- grep("save.edf", snames)
        nic <- snames[(ic - 1):ic]
        IC <- samples[[j]][, grepl(nic[1], snames)][1]
        edf <- samples[[j]][, grepl(nic[2], snames)][1]
        names(IC) <- nic[1]
        DIC <- pd <- NULL
      } else {
        IC <- edf <- NULL
      }

      ## Compute model term effects.
      param.effects <- effects <- effects.hyp <- NULL
      fitted.values <- 0

      ## Smooth terms.
      if(length(obj$smooth)) {
        for(i in 1:length(obj$smooth)) {
          ## Get coefficient samples of smooth term.
          if(any(grepl("h1", snames))) {
            pn <- grep(paste(id, "h1", obj$smooth[[i]]$label, sep = ":"), snames,
              value = TRUE, fixed = TRUE) ## FIXME: hlevels!
          } else {
            pn <- grep(paste(id, obj$smooth[[i]]$label, sep = "."), snames,
              value = TRUE, fixed = TRUE)
          }
          pn <- pn[!grepl("tau2", pn) & !grepl("alpha", pn) & !grepl("accepted", pn)]
          k <- sum(snames %in% pn)
          psamples <- matrix(samples[[j]][, snames %in% pn], ncol = k)
          nas <- apply(psamples, 1, function(x) { any(is.na(x)) } )
          psamples <- psamples[!nas, , drop = FALSE]

          if(!is.null(obj$smooth[[i]]$Xf)) {
            kx <- if(is.null(obj$smooth[[i]]$Xf)) 0 else ncol(obj$smooth[[i]]$Xf)
            if(kx) {
              pn <- paste(paste(id, ":h1:linear.",
                paste(paste(obj$smooth[[i]]$term, collapse = "."), "Xf", sep = "."), sep = ""),
                1:kx, sep = ".")
              xsamples <- matrix(samples[[j]][, snames %in% pn], ncol = kx)
              psamples <- cbind("ra" = psamples, "fx" = xsamples)
              re_trans <- function(g) {
                g <- obj$smooth[[i]]$trans.D * g
                if(!is.null(obj$smooth[[i]]$trans.U))
                  g <- obj$smooth[[i]]$trans.U %*% g
                g
              }
              psamples <- t(apply(psamples, 1, re_trans))
            }
          }

          ## Possible variance parameter samples.
          vsamples <- NULL
          tau2 <- paste(id, "h1", paste(obj$smooth[[i]]$label, "tau2", sep = "."), sep = ":")

          if(length(tau2 <- grep(tau2, snames, fixed = TRUE))) {
            vsamples <- samples[[j]][, tau2, drop = FALSE]
            vsamples <- vsamples[!nas, , drop = FALSE]
          }

          ## Acceptance probalities.
          asamples <- NULL
          alpha <- paste(id, "h1", paste(obj$smooth[[i]]$label, "alpha", sep = "."), sep = ":")
          if(length(alpha <- grep(alpha, snames, fixed = TRUE))) {
            asamples <- as.numeric(samples[[j]][, alpha])
            asamples <- asamples[!nas]
          }

          ## Prediction matrix.
          get.X <- function(x) { ## FIXME: time(x)
            acons <- obj$smooth[[i]]$xt$center
            for(char in c("(", ")", "[", "]"))
              obj$smooth[[i]]$term <- gsub(char, ".", obj$smooth[[i]]$term, fixed = TRUE)
            X <- PredictMat(obj$smooth[[i]], x)
            X
          }

          ## Compute final smooth term object.
          tn <- c(obj$smooth[[i]]$term, if(obj$smooth[[i]]$by != "NA") obj$smooth[[i]]$by else NULL)

          if(is.null(obj$smooth[[i]]$is.parametric)) {
            if(!is.list(effects))
              effects <- list()
            if(length(effects)) {
              if(obj$smooth[[i]]$label %in% names(effects)) {
                ct <- gsub(".smooth.spec", "", class(obj$smooth[[i]]))[1]
                if(ct == "random.effect") ct <- "re"
                obj$smooth[[i]]$label <- paste(obj$smooth[[i]]$label, ct, sep = ":")
              }
            }
            if(is.null(obj$smooth[[i]]$get.mu)) {
              obj$smooth[[i]]$get.mu <- function(X, b) {
                drop(X %*% b)
              }
            }

            fst <- compute_term(obj$smooth[[i]], get.X = get.X, get.mu = obj$smooth[[i]]$get.mu,
              psamples = psamples, vsamples = vsamples, asamples = asamples,
              FUN = NULL, snames = snames, effects.hyp = effects.hyp,
              fitted.values = fitted.values, data = attr(x, "model.frame")[, tn, drop = FALSE],
              grid = grid)

            attr(fst$term, "specs")$get.mu <- obj$smooth[[i]]$get.mu

            ## Add term to effects list.
            effects[[obj$smooth[[i]]$label]] <- fst$term
            effects.hyp <- fst$effects.hyp

            fitted.values <- fst$fitted.values
            rm(fst)
          } else {
            nx <- colnames(obj$smooth[[i]]$X)
            qu <- t(apply(psamples, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
            sd <- drop(apply(psamples, 2, sd, na.rm = TRUE))
            me <- drop(apply(psamples, 2, mean, na.rm = TRUE))
            param.effects <- cbind(me, sd, qu)
            rownames(param.effects) <- nx
            colnames(param.effects) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
            if(!is.null(asamples))
              param.effects <- cbind(param.effects, "alpha" = mean(asamples))
            fitted.values <- as.vector(fitted.values + obj$smooth[[i]]$X %*% param.effects[, 1])
            attr(param.effects, "samples") <- as.mcmc(psamples)
            colnames(attr(param.effects, "samples")) <- nx
          }
        }
      }

      ## Compute partial residuals.
      if(!is.null(effects)) {
        if(length(obj$response)) {
          if(obj$response %in% names(attr(x, "model.frame"))) {
            effects <- partial.residuals(effects, attr(x, "model.frame")[[obj$response]],
              fitted.values, family)
          }
        }
      }

      ## Stuff everything together.
      rval[[j]] <- list(
        "model" = list("DIC" = DIC, "pd" = pd,
          "N" = nrow(attr(x, "model.frame")), "formula" = obj$formula,
          "IC" = IC, "edf" = edf),
        "param.effects" = param.effects, "effects" = effects,
        "effects.hyp" = effects.hyp, "fitted.values" = fitted.values
      )
      
#      ## Clean.
#      rval[[j]] <- delete.NULLs(rval[[j]])

      class(rval[[j]]) <- "bamlss"
    }
    names(rval) <- paste("Chain", 1:chains, sep = "_")
    if(length(rval) < 2) {
      rval <- rval[[1]]
    }
    class(rval) <- "bamlss"
    return(rval)
  }

  if(inherits(x, "bamlss.input") & !all(c("formula", "fake.formula", "response") %in% names(x))) {
    nx <- names(x)
    nx <- nx[nx != "call"]
    rval <- list()
    fn <- family$names
    cat <- if(!is.null(family$cat)) family$cat else FALSE
    if(cat) {
      if(length(attr(attr(x, "model.frame"), "response.name")) < 2)
        fn <- gsub(attr(attr(x, "model.frame"), "response.name"), "", names(x))
    }
    if(length(fn) != length(nx))
      fn <- paste(fn, 1:length(nx), sep = "")
    for(j in seq_along(nx)) {
      rval[[nx[j]]] <- createBayesGresults(x[[nx[j]]], samples, id = fn[j])
      if(!is.null(rval[[nx[j]]]$effects)) {
        for(i in seq_along(rval[[nx[j]]]$effects)) {
          specs <- attr(rval[[nx[j]]]$effects[[i]], "specs")
          specs$label <- paste(specs$label, fn[j], sep = ":")
          attr(rval[[nx[j]]]$effects[[i]], "specs") <- specs
        }
        names(rval[[nx[j]]]$effects) <- paste(names(rval[[nx[j]]]$effects), fn[j], sep = ":")
      }
    }
    names(rval) <- fn
    if(cat) {
      reference <- attr(x, "reference")
      if(!is.null(reference)) {
        if(any(reference %in% names(rval)))
          rval <- rval[!grepl(reference, fn)]
      }
    }
    attr(rval, "family") <- family
    class(rval) <- "bamlss"
    return(rval)
  } else {
    return(createBayesGresults(x, samples))
  }
}

