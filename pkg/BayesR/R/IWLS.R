## Some useful links:
## http://adv-r.had.co.nz/C-interface.html
## http://stackoverflow.com/questions/7457635/calling-r-function-from-c
## http://gallery.rcpp.org/articles/r-function-from-c++/
## IWLS propose function in C.
propose_iwls <- function(x, family,
  response, eta, id, rho, ...)
{
  .Call("do_propose", x, family, response, eta, id, rho)
}

propose_iwls0 <- function(x, family, response, eta, id, ...)
{
  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute weights.
  weights <- family$iwls$weights[[id]](response, peta)

  ## Score.
  score <- family$iwls$score[[id]](response, peta)

  ## Compute working observations.
  z <- eta[[id]] + 1 / weights * score

  ## Compute old log likelihood and old log coefficients prior.
  pibeta <- family$loglik(response, peta)
  p1 <- if(x$fixed) {
    0
  } else drop(-0.5 / x$state$tau2 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g)

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit

  ## Compute mean and precision.
  XW <- t(x$X * weights)
  P <- if(x$fixed) {
    if(k <- ncol(x$X) < 2) {
      1 / (XW %*% x$X)
    } else chol2inv(chol(XW %*% x$X))
  } else chol2inv(chol(XW %*% x$X + 1 / x$state$tau2 * x$S[[1]]))
  P[P == Inf] <- 0
  M <- P %*% (XW %*% (z - eta[[id]]))

  ## Save old coefficients
  g0 <- drop(x$state$g)

  ## Sample new parameters.
  x$state$g <- drop(rmvnorm(n = 1, mean = M, sigma = P))

  ## Compute log priors.
  p2 <- if(x$fixed) {
    0
  } else drop(-0.5 / x$state$tau2 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g)
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
  weights <- family$iwls$weights[[id]](response, peta)

  ## New score.
  score <- family$iwls$score[[id]](response, peta)

  ## New working observations.
  z <- eta[[id]] + 1 / weights * score

  ## Compute mean and precision.
  XW <- t(x$X * weights)
  P2 <- if(x$fixed) {
    if(k < 2) {
      1 / (XW %*% x$X)
    } else chol2inv(chol(XW %*% x$X))
  } else chol2inv(L <- chol(P0 <- XW %*% x$X + 1 / x$state$tau2 * x$S[[1]]))
  P2[P2 == Inf] <- 0
  M2 <- P2 %*% (XW %*% (z - (eta[[id]] - x$state$fit)))

  ## Get the log prior.
  qbeta <- dmvnorm(g0, mean = M2, sigma = P2, log = TRUE)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp)) {
    a <- x$rank / 2 + x$a
    b <- 0.5 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g + x$b
    x$state$tau2 <- 1 / rgamma(1, a, b)
  }

  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

  return(x$state)
}


## Random walk propose function.
propose_rw <- function(x, family,
  response, eta, id, ...)
{
  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute old log likelihood and old log coefficients prior.
  pibeta <- family$loglik(response, peta)
  p1 <- if(x$fixed) {
    sum(dnorm(x$state$g, sd = 10, log = TRUE), na.rm = TRUE)
  } else drop(-0.5 / x$state$tau2 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g)

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit

  ## Number of parameters.
  k <- length(x$state$g)

  ## Sample new parameters.
  do <- TRUE; accepted <- x$state$accepted
  eta2 <- eta; j <- 0
  while(do & j < x$state$maxit) {
    ## Adaptive scales.
    if(!is.null(x$state$iter)) {
      aprop <- 0.44
      if(x$state$iter < x$adapt) {
        for(i in 1:k) {
          x$state$scale[i] <- if(accepted) {
            x$state$scale[i] + (x$state$scale[i] / (aprop * (1 - aprop))) * (1 - aprop) / x$state$iter
          } else {
            abs(x$state$scale[i] - (x$state$scale[i] / (aprop * (1 - aprop))) * aprop / x$state$iter)
          }
        }
      }
    }

    g <- if(TRUE) {
      theta <- rnorm(k)
      d <- theta / sqrt(sum(theta * theta))
      x$state$g + runif(k, 0, x$state$scale) * d
    } else x$state$g + rnorm(k, mean = 0, sd = x$state$scale)

    ## Compute log priors.
    p2 <- if(x$fixed) {
      sum(dnorm(g, sd = 10, log = TRUE), na.rm = TRUE)
    } else drop(-0.5 / x$state$tau2 * crossprod(g, x$S[[1]]) %*% g)
  
    ## Compute fitted values.        
    fit <- drop(x$X %*% g)

    ## Set up new predictor.
    eta2[[id]] <- eta[[id]] + fit

    ## Map predictor to parameter scale.
    peta <- family$map2par(eta2)

    ## Compute new log likelihood.
    pibetaprop <- family$loglik(response, peta)

    ## Compute log acceptance probablity.
    x$state$alpha <- drop((pibetaprop + p2) - (pibeta + p1))

    if((accepted <- if(is.na(x$state$alpha)) FALSE else log(runif(1)) <= x$state$alpha) | x$state$maxit < 2) {
      x$state$g <- g
      x$state$fit <- fit
      do <- FALSE
    }

    j <- j + 1
  }

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp)) {
    a <- x$rank / 2 + x$a
    b <- 0.5 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g + x$b
    x$state$tau2 <- 1 / rgamma(1, a, b)
  }

  return(x$state)
}


## t-walk propose function.
## Functions for the OneMove().
## see http://www.cimat.mx/~jac/twalk/
## Author J. Andres Christen
IntProd <- function(x) { sum(x * x)  }
DotProd <- function(x, y) { sum(x * y)  }

Simh1 <- function(dim, pphi, x, xp, beta)
{	
	phi <- (runif(dim) < pphi) 
	rt <- NULL 
	for(i in 1:dim)
		if(phi[i])
			rt <- append(rt, xp[i] + beta*(xp[i] - x[i]))
		else
			rt <- append(rt, x[i]) 
	list(rt = rt, nphi = sum(phi)) 
}

Simfbeta <- function(at)
{	
	if(runif(1) < (at - 1) / (2 * at))
		exp(1 / (at + 1) * log(runif(1)))
	else
		exp(1 / (1 - at) * log(runif(1))) 
}

Simh2 <- function(dim, pphi, aw, x, xp)
{
  u <- runif(dim) 
  phi <- (runif(dim) < pphi)
  z <- (aw / (1 + aw)) * (aw * u^2 + 2*u - 1) 
  z <- z * phi
  list(rt = x + (x - xp) * z, nphi = sum(phi))
}

Simh3 <- function(dim, pphi, x, xp)
{
  phi <- (runif(dim) < pphi)
  sigma <- max(phi * abs(xp - x))
  x + sigma * rnorm(dim) * phi
  list(rt = xp * phi + sigma * rnorm(dim) * phi + x * (1 - phi), nphi = sum(phi), phi = phi)
}

G3U <- function(nphi, phi, h, x, xp)
{
  sigma <- max(phi * abs(xp - x))
  if(nphi > 0)
    (nphi / 2) * log(2 * pi) + nphi * log(sigma) + 0.5 * IntProd(h - xp) / (sigma^2)
  else
    0 
}

Simh4 <- function(dim, pphi, x, xp)
{
  phi <- (runif(dim) < pphi)
	sigma <- max(phi * abs(xp - x)) / 3
	rt <- NULL
  for(i in 1:dim)
    if(phi[i])
      rt <- append(rt, x[i] + sigma * rnorm(1))
    else
      rt <- append(rt, x[i])
  list(rt = rt, nphi = sum(phi), phi = phi) 
}

log2pi <- log(2 * pi); log3 <- log(3)

G4U <- function(nphi, phi, h, x, xp)
{
  sigma <- max(phi * abs(xp - x)) / 3
  if(nphi > 0)
    (nphi / 2) * log2pi - nphi * log3 + nphi * log(sigma) + 0.5 * 9 * IntProd((h - x)) / (sigma^2)
  else
    0 
}

OneMove <- function(dim, Supp = function(x) { TRUE },
  x, U, xp, Up, at = 6, aw = 1.5, pphi = min(dim, 4) / dim,
  F1 = 0.4918, F2 = 0.9836, F3 = 0.9918,
  .X, .FAMILY, .RESPONSE, .ETA, .ID)
{
  ker <- runif(1)

  if(ker < F1) {	
    dir <- runif(1) 
    funh <- 1
    if((0 <= dir) && (dir < 0.5)) {	
      beta <- Simfbeta(at) 
      tmp <- Simh1(dim, pphi, xp, x, beta) 
      yp <- tmp$rt
      nphi <- tmp$nphi
      y  <- x 
      propU <- U 
      if(Supp(yp)) {
        propUp <- logPost(yp, .X, .FAMILY, .RESPONSE, .ETA, .ID)
        if(nphi == 0)
          A <- 1
        else
          A <- exp((U - propU) + (Up - propUp) +  (nphi-2)*log(beta)) 
      } else {
        propUp <- NULL
        A <- 0
      }
    }
    if((0.5 <= dir) && (dir < 1.0)) {
      beta <- Simfbeta(at)
      tmp <- Simh1(dim, pphi, x, xp, beta) 
      y <- tmp$rt
      nphi <- tmp$nphi
      yp  <- xp 
      propUp <- Up 
      if(Supp(y)) {
        propU <- logPost(y, .X, .FAMILY, .RESPONSE, .ETA, .ID)
        if(nphi == 0)
          A <- 1
        else
          A <- exp((U - propU) + (Up - propUp) +  (nphi - 2) * log(beta))
      } else {
        propU <- NULL
        A <- 0
      }
    }
  }
	
  if((F1 <= ker) && (ker < F2)) {	
    dir <- runif(1)
    funh <- 2
    if((0 <= dir) && (dir < 0.5)) {
      tmp <- Simh2(dim, pphi, aw, xp, x) 
      yp <- tmp$rt
      nphi <- tmp$nphi
      y  <- x 
      propU <- U
      if((Supp(yp)) && (all(abs(yp - y) > 0))) {
        propUp <- logPost(yp, .X, .FAMILY, .RESPONSE, .ETA, .ID) 
        A <- exp((U - propU) + (Up - propUp))  
      } else {
        propUp <- NULL
        A <- 0
      }
    }
    if((0.5 <= dir) && (dir < 1.0)) {
      tmp <- Simh2(dim, pphi, aw, x, xp)
      y <- tmp$rt
      nphi <- tmp$nphi
      yp  <- xp 
      propUp <- Up
      if((Supp(y)) && (all(abs(yp - y) > 0))) {
        propU <- logPost(y, .X, .FAMILY, .RESPONSE, .ETA, .ID)
        A <- exp((U - propU) + (Up - propUp))
      } else {
        propU <- NULL
        A <- 0
      }
    }
  }

  if((F2 <= ker) && (ker < F3)) {	
    dir <- runif(1) 
    funh <- 3
    if((0 <= dir) && (dir < 0.5)) {
      tmp <- Simh3( dim, pphi, xp, x) 
      yp <- tmp$rt
      nphi <- tmp$nphi
      phi <- tmp$phi
      y  <- x 
      propU <- U 
      if((Supp(yp)) && all(yp != x)) {
        propUp <- logPost(yp, .X, .FAMILY, .RESPONSE, .ETA, .ID) 
        W1 <- G3U(nphi, phi, yp, xp,  x) 
        W2 <- G3U(nphi, phi, xp, yp,  x)  
        A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))
      } else {
        propUp <- NULL
        A <- 0
      }
    }
    if((0.5 <= dir) && (dir < 1.0)) {
      tmp <- Simh3(dim, pphi, x, xp) 
      y <- tmp$rt
      nphi <- tmp$nphi
      phi <- tmp$phi
      yp  <- xp 
      propUp <- Up
      if((Supp(y)) && all(y != xp)) {
        propU <- logPost(y, .X, .FAMILY, .RESPONSE, .ETA, .ID) 
        W1 <- G3U(nphi, phi, y,  x, xp) 
        W2 <- G3U(nphi, phi, x,  y, xp) 
        A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))
      } else {
        propU <- NULL
        A <- 0
      }
    }
  }
		
  if(F3 <= ker) {	
    dir <- runif(1) 
    funh <- 4
    if((0 <= dir) && (dir < 0.5)) {
      tmp <- Simh4(dim, pphi, xp, x) 
      yp <- tmp$rt
      nphi <- tmp$nphi
      phi <- tmp$phi
      y  <- x 
      propU <- U
      if((Supp(yp)) && all(yp != x)) {
        propUp <- logPost(yp, .X, .FAMILY, .RESPONSE, .ETA, .ID) 
        W1 <- G4U(nphi, phi, yp, xp,  x) 
        W2 <- G4U(nphi, phi, xp, yp,  x) 
        A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))
      } else {
        propUp <- NULL
        A <- 0
      }
    }
    if((0.5 <= dir) && (dir < 1.0)) {
      tmp <- Simh4(dim, pphi, x, xp) 
      y <- tmp$rt
      nphi <- tmp$nphi
      phi <- tmp$phi
      yp  <- xp 
      propUp <- Up
      if((Supp(y)) && all(y != xp)) {
        propU <- logPost(y, .X, .FAMILY, .RESPONSE, .ETA, .ID) 
        W1 <- G4U(nphi, phi, y,  x, xp) 
        W2 <- G4U( nphi, phi, x,  y, xp) 
        A <- exp((U - propU) + (Up - propUp) +  (W1 - W2))
      } else {
        propU <- NULL
        A <- 0
      }
    }
  }
	
  return(list(y = y, propU = propU, yp = yp,
    propUp = propUp, A = A, funh = funh, nphi = nphi))
}

logPost <- function(g, x, family, response, eta, id)
{
  ## Set up new predictor.
  eta[[id]] <- eta[[id]] + drop(x$X %*% g)

  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute log likelihood and log coefficients prior.
  ll <- family$loglik(response, peta)
  lp <- if(x$fixed) {
    sum(dnorm(g, sd = 10, log = TRUE), na.rm = TRUE)
  } else drop(-0.5 / x$state$tau2 * crossprod(g, x$S[[1]]) %*% g)

  -1 * (ll + lp)
}

propose_twalk <- function(x, family,
  response, eta, id, ...)
{
  ## Remove fitted values.
  eta[[id]] <- eta[[id]] - x$state$fit

  ## Number of parameters.
  k <- length(x$state$g)

  ## Get log posteriors if not available.
  if(is.null(x$state$U))
    x$state$U <- logPost(x$state$g, x, family, response, eta, id)
  if(is.null(x$state$sg))
    x$state$sg <- rep(0, length = k)
  if(is.null(x$state$Up))
    x$state$Up <- logPost(x$state$sg, x, family, response, eta, id)

  ## Do one t-walk step
  p <- OneMove(dim = k, x = x$state$g, U = x$state$U, xp = x$state$sg, Up = x$state$Up,
    .X = x, .FAMILY = family, .RESPONSE = response, .ETA = eta, .ID = id)

  ## Setup return state.
  x$state$alpha <- log(p$A)
  x$state$g <- p$y
  x$state$sg <- p$yp
  x$state$U <- p$propU
  x$state$Up <- p$propUp
  x$state$fit <- drop(x$X %*% x$state$g)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp)) {
    a <- x$rank / 2 + x$a
    b <- 0.5 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g + x$b
    x$state$tau2 <- 1 / rgamma(1, a, b)
  }

  return(x$state)
}


## Numerical derivatives.
num_deriv <- function(y, eta, family, id = NULL,
  d = 1, method = "simple", eps = 1e-4)
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


## Slice sampling.
## See: http://www.cs.toronto.edu/~radford/ftp/slice-R-prog
## Log-posterior used by propose_slice().
logPost2 <- function(g, x, family, response, eta, id)
{
  ## Set up new predictor.
  eta[[id]] <- eta[[id]] + drop(x$X %*% g)

  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute log likelihood and log coefficients prior.
  ll <- family$loglik(response, peta)
  lp <- if(x$fixed) {
    sum(dnorm(g, sd = 10, log = TRUE), na.rm = TRUE)
  } else drop(-0.5 / x$state$tau2 * crossprod(g, x$S[[1]]) %*% g)

  return(ll + lp)
}

logPost3 <- function(g, x, family, response, eta, id)
{
  ## Set up new predictor.
  fit <- x$get.mu(x$X, g)
  fit <- fit - mean(fit, na.rm = TRUE)
  eta[[id]] <- eta[[id]] + fit

  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute log likelihood and log coefficients prior.
  ll <- family$loglik(response, peta)
  lp <- sum(dnorm(g, sd = 10, log = TRUE), na.rm = TRUE)

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
  w = 1, m = 100, lower = -Inf, upper = +Inf, logPost)
{
  x0 <- g[j]
  gL <- gR <- g

  gx0 <- logPost(g, x, family, response, eta, id)

  ## Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)

  ## Find the initial interval to sample from.
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
  ## Remove fitted values.
  eta[[id]] <- eta[[id]] - x$state$fit

  for(j in seq_along(x$state$g)) {
    x$state$g <- uni.slice(x$state$g, x, family, response, eta, id, j,
      logPost = logPost2, rho = rho)
  }

  ## Setup return state.
  x$state$alpha <- log(1)
  x$state$fit <- drop(x$X %*% x$state$g)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp)) {
    a <- x$rank / 2 + x$a
    b <- 0.5 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g + x$b
    x$state$tau2 <- 1 / rgamma(1, a, b)
  }

  return(x$state)
}


propose_wslice <- function(x, family,
  response, eta, id, rho, ...)
{
  args <- list(...)
  iter <- args$iter  

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit

  if(iter %% 20 == 0) {
    ## Map predictor to parameter scale.
    peta <- family$map2par(eta)

    ## Compute weights.
    weights <- family$iwls$weights[[id]](response, peta)

    ## Score.
    score <- family$iwls$score[[id]](response, peta)

    ## Compute working observations.
    z <- eta[[id]] + 1 / weights * score

    ## Compute mean and precision.
    XW <- t(x$X * weights)
    P <- if(x$fixed) {
      if(k <- ncol(x$X) < 2) {
        1 / (XW %*% x$X)
      } else chol2inv(chol(XW %*% x$X))
    } else chol2inv(chol(XW %*% x$X + 1 / x$state$tau2 * x$S[[1]]))
    P[P == Inf] <- 0
    x$state$g <- drop(P %*% (XW %*% (z - eta[[id]])))
  }

  for(j in seq_along(x$state$g)) {
    x$state$g <- uni.slice(x$state$g, x, family, response, eta, id, j,
      logPost = logPost2, rho = rho)
  }

  ## Setup return state.
  x$state$alpha <- log(2)
  x$state$fit <- drop(x$X %*% x$state$g)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp)) {
    a <- x$rank / 2 + x$a
    b <- 0.5 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g + x$b
    x$state$tau2 <- 1 / rgamma(1, a, b)
  }

  return(x$state)
}


propose_oslice <- function(x, family,
  response, eta, id, rho, ...)
{
  args <- list(...)
  iter <- args$iter  

  if(iter %% 20 == 0) {
    ## Compute mean.
    opt <- update_optim2(x, family, response, eta, id, ...)
    x$state$g <- opt$g
  }

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit

  for(j in seq_along(x$state$g)) {
    x$state$g <- uni.slice(x$state$g, x, family, response, eta, id, j,
      logPost = logPost2, rho = rho)
  }

  ## Setup return state.
  x$state$alpha <- log(2)
  x$state$fit <- drop(x$X %*% x$state$g)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp)) {
    a <- x$rank / 2 + x$a
    b <- 0.5 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g + x$b
    x$state$tau2 <- 1 / rgamma(1, a, b)
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
    weights <- family$iwls$weights[[id]](response, peta)
  } else weights <- args$weights

  ## Which obs to take.
  ok <- !(weights %in% c(NA, -Inf, Inf, 0))
  weights <- weights[ok]

  if(is.null(args$z)) {
    ## Score.
    score <- family$iwls$score[[id]](response, peta)

    ## Compute working observations.
    z <- eta[[id]][ok] + 1 / weights * score[ok]
  } else z <- args$z[ok]

  ## Compute partial predictor.
  eta[[id]][ok] <- eta[[id]][ok] - x$state$fit[ok]

  ## Compute mean and precision.
  XW <- t(x$X[ok, ] * weights)
  XWX <- XW %*% x$X[ok, ]
  if(is.null(x$optimize) | x$fixed | !is.null(x$sp)) {
    P <- if(x$fixed) {
      matrix_inv(XWX)
    } else matrix_inv(XWX + 1 / x$state$tau2 * x$S[[1]])
    x$state$g <- drop(P %*% (XW %*% (z - eta[[id]][ok])))
  } else {
    args <- list(...)
    edf0 <- args$edf - x$state$edf
    eta2 <- eta
    e <- z - eta[[id]][ok]

    objfun <- function(tau2, ...) {
      P <- matrix_inv(XWX + 1 / tau2 * x$S[[1]])
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

    x$state$tau2 <- try(optimize(objfun, interval = x$interval, grid = x$grid)$minimum, silent = TRUE)
    if(inherits(x$state$tau2, "try-error"))
      x$state$tau2 <- optimize2(objfun, interval = x$interval, grid = x$grid)$minimum
    if(!length(x$state$tau2)) x$state$tau2 <- x$interval[1]
    P <- matrix_inv(XWX + 1 / x$state$tau2 * x$S[[1]])
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

  if(!x$fixed) {
    ## Objective function for variance parameter.
    objfun2 <- function(tau2) {

      ## Objective for regression coefficients.
      objfun <- function(gamma) {
        eta2[[id]] <- eta[[id]] + drop(x$X %*% gamma)
        ll <- family$loglik(response, family$map2par(eta2))
        lp <- drop(-0.5 / tau2 * crossprod(gamma, x$S[[1]]) %*% gamma)
        -1 * (ll + lp)
      }

      suppressWarnings(opt <- try(optim(x$state$g, fn = objfun, method = "BFGS",
        control = list()), silent = TRUE))

      if(!inherits(opt, "try-error")) {
        x$state$g <- opt$par
        x$state$fit <- drop(x$X %*% x$state$g)
      }

      P <- matrix_inv(x$state$XX + 1 / tau2 * x$S[[1]])
      if(inherits(P, "try-error")) return(NA)
      edf <- sum(diag(x$state$XX %*% P))
      if(!is.null(x$xt$center)) {
        if(x$xt$center) edf <- edf - 1
      }

      eta2[[id]] <- eta[[id]] + x$state$fit

      IC <- get.ic(family, response, family$map2par(eta2), edf0 + edf, length(eta2[[id]]), x$criterion)
      IC
    }

    x$state$tau2 <- try(optimize(objfun2, interval = x$interval, grid = x$grid)$minimum, silent = TRUE)
    if(inherits(x$state$tau2, "try-error"))
      x$state$tau2 <- optimize2(objfun2, interval = x$interval, grid = x$grid)$minimum
    if(!length(x$state$tau2)) x$state$tau2 <- x$interval[1]
  }

  ## Objective for regression coefficients.
  objfun <- function(gamma) {
    eta2[[id]] <- eta[[id]] + drop(x$X %*% gamma)
    ll <- family$loglik(response, family$map2par(eta2))
    lp <- if(x$fixed) {
      sum(dnorm(gamma, sd = 10, log = TRUE))
    } else drop(-0.5 / x$state$tau2 * crossprod(gamma, x$S[[1]]) %*% gamma)
    -1 * (ll + lp)
  }

  suppressWarnings(opt <- try(optim(x$state$g, fn = objfun, method = "BFGS",
    control = list()), silent = TRUE))

  if(!inherits(opt, "try-error")) {
    x$state$g <- opt$par
    x$state$fit <- drop(x$X %*% x$state$g)
  }

  if(!x$fixed) {
    P <- matrix_inv(x$state$XX + 1 / x$state$tau2 * x$S[[1]])
    x$state$edf <- sum(diag(x$state$XX %*% P))
    if(!is.null(x$xt$center)) {
      if(x$xt$center) x$state$edf <- x$state$edf - 1
    }
  }

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
    eta2[[id]] <- eta[[id]] + drop(x$X %*% gamma)
    peta <- family$map2par(eta2)
    ll <- family$loglik(response, peta)
    lp <- if(x$fixed) {
      sum(dnorm(gamma, sd = 10, log = TRUE))
    } else drop(-0.5 / x$state$tau2 * crossprod(gamma, x$S[[1]]) %*% gamma)
    -1 * (ll + lp)
  }

  suppressWarnings(opt <- try(optim(x$state$g, fn = objfun, method = "BFGS",
    control = list()), silent = TRUE))

  if(!inherits(opt, "try-error")) {
    x$state$g <- opt$par
    x$state$fit <- drop(x$X %*% x$state$g)
  }

  if(!x$fixed) {
    if(is.null(x$state$XX))
      x$state$XX <- crossprod(x$X)
    P <- matrix_inv(x$state$XX + 1 / x$state$tau2 * x$S[[1]])
    if(!inherits(P, "try-error")) {
      x$state$edf <- sum(diag(x$state$XX %*% P))
      if(!is.null(x$xt$center)) {
        if(x$xt$center) x$state$edf <- x$state$edf - 1
      }
    }
  }

  return(x$state)
}


#############################
## Nadja's testing section ##
#############################
propose_nadja <- function(x, family,
  response, eta, id, ...)
{
  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute weights.
  weights <- family$iwls$weights[[id]](response, peta)

  ## Score.
  score <- family$iwls$score[[id]](response, peta)

  ## Compute working observations.
  z <- eta[[id]] + 1 / weights * score

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - x$state$fit

  ## Compute mean and precision.
  XW <- t(x$X * weights)
  P <- if(x$fixed) {
    if(k <- ncol(x$X) < 2) {
      1 / (XW %*% x$X)
    } else chol2inv(chol(XW %*% x$X))
  } else chol2inv(chol(XW %*% x$X + 1 / x$state$tau2 * x$S[[1]]))
  P[P == Inf] <- 0
  x$state$g <- drop(P %*% (XW %*% (z - eta[[id]])))

  for(j in seq_along(x$state$g)) {
    x$state$g <- uni.slice(x$state$g, x, family, response, eta, id, j, logPost = logPost2)
  }

  ## Setup return state.
  x$state$alpha <- log(2)
  x$state$fit <- drop(x$X %*% x$state$g)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp)) {
    a <- x$rank / 2 + x$a
    b <- 0.5 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g + x$b
    x$state$tau2 <- 1 / rgamma(1, a, b)
  }

  return(x$state)
}


if(FALSE) {
  require("mgcv")
  set.seed(111)
  n <- 300
  z <- runif(n, -3, 3)
  response <- 1.2 + sin(z) + rnorm(n, sd = 0.6)
  x <- smooth.construct(s(z), list("z" = z), NULL)
  x$state <- list("g" = runif(ncol(x$X)), "tau2" = 2.33, "fit" = runif(n))
  x$a <- x$b <- 0.00001
  family <- bayesr.family(gaussian.BayesR())
  eta <- list("mu" = rep(0, n), "sigma" = rep(0, n))
  id <- "mu"

  a <- propose_rw(x, family, response, eta, id, new.env())
}


## Setup for IWLS sampler, handling
## sampling functions.
transformIWLS <- function(x, ...)
{
  call <- x$call; x$call <- NULL
  x <- assign.weights(x)

  tIWLS <- function(obj, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(obj)
      if(length(unique(nx)) < length(obj)) nx <- 1:length(obj)
      for(j in nx)
        obj[[j]] <- tIWLS(obj[[j]], ...)
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
            "is.linear" = TRUE
          )
          obj$sterms <- c(obj$strems, "parametric")
          obj$X <- NULL
        }
      }
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          obj$smooth[[j]] <- smooth.IWLS(obj$smooth[[j]])
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
                "is.linear" = TRUE
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

  x <- tIWLS(x, ...)

  attr(x, "call") <- call
  attr(x, "response.vec") <- attr(x, "model.frame")[, attr(attr(x, "model.frame"), "response.name")]

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


## Function to setup IWLS smooths.
smooth.IWLS <- function(x, ...) {
  UseMethod("smooth.IWLS")
}

smooth.IWLS.default <- function(x, ...)
{
  if(is.null(x$get.mu)) {
    x$get.mu <- if(is.null(x$get.mu)) {
      function(X, g) {
        X %*% as.numeric(g)
      }
    } else x$get.mu
  }

  x$a <- if(is.null(x$xt$a)) 1e-04 else x$xt$a
  x$b <- if(is.null(x$xt$b)) 1e-04 else x$xt$b

  tau2interval <- function(x, lower = .Machine$double.eps^0.25, upper = 4000) {
    XX <- crossprod(x$X)
    objfun <- function(tau2, value) {
      df <- sum(diag(matrix_inv(XX + if(x$fixed) 0 else 1 / tau2 * x$S[[1]]) %*% XX))
      return((value - df)^2)
    }
    le <- try(optimize(objfun, c(lower, upper), value = 1)$minimum, silent = TRUE)
    ri <- try(optimize(objfun, c(lower, upper), value = ncol(x$X))$minimum, silent = TRUE)
    if(inherits(le, "try-error")) le <- 0.1
    if(inherits(ri, "try-error")) ri <- 1000
    return(c(le, ri))
  }

  x$interval <- if(is.null(x$xt$interval)) tau2interval(x) else x$xt$interval
  x$grid <- if(is.null(x$xt$grid)) 40 else x$xt$grid
  if(is.null(x$fixed)) {
    x$fixed <- if(!is.null(x$fx))  x$fx[1] else FALSE
  }

  if(is.null(x$state)) {
    x$p.save <- c("g", if(!x$fixed) "tau2" else NULL)
    x$state <- list()
    x$state$g <- rep(0, ncol(x$X))
    x$state$tau2 <- if(is.null(x$sp)) {
      if(x$fixed) 1e-20 else 10
    } else x$sp
    if(x$fixed) x$state$tau2 <- 1e-20
    x$s.colnames <- if(is.null(x$s.colnames)) {
      c(paste("c", 1:length(x$state$g), sep = ""),
        if(x$fixed) NULL else rep("tau2", length = length(x$state$tau2)))
    } else x$s.colnames
    x$np <- c(length(x$state$g), length(x$state$tau2))
    XX <- crossprod(x$X)
    x$state$edf <- sum(diag(matrix_inv(XX + if(x$fixed) 0 else 1 / x$state$tau2 * x$S[[1]]) %*% XX))
  }
  if(is.null(x$xt$adaptive))
    x$xt$adaptive <- TRUE

  x
}


get.ic <- function(family, response, eta, edf, n, type = c("AIC", "BIC", "AICc"))
{
  type <- match.arg(type)
  ll <- family$loglik(response, eta)
  pen <- switch(type,
    "AIC" = -2 * ll + 2 * edf,
    "BIC" = -2 * ll + edf * log(n),
    "AICc" = -2 * ll + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1)
  )
  return(pen)
}


## Sampler based on IWLS proposals.
samplerIWLS <- function(x, n.iter = 12000, thin = 10, burnin = 2000, accept.only = FALSE,
  verbose = TRUE, step = 20, svalues = TRUE, eps = .Machine$double.eps^0.25, maxit = 400,
  tdir = NULL, method = "backfitting", outer = FALSE, inner = FALSE, n.samples = 200,
  criterion = c("AICc", "BIC", "AIC"), lower = 1e-09, upper = 1e+04,
  optim.control = list(pgtol = 1e-04, maxit = 5), digits = 3,
  update = c("optim2", "iwls", "optim"),
  propose = c("oslice", "iwls", "rw", "twalk", "slice", "wslice", "nadja"),
  sample = c("slice", "iwls"), ...)
{
  known_methods <- c("backfitting", "MCMC", "backfitting2", "backfitting3", "backfitting4")
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
  scipen <- getOption("scipen")
  options("scipen" = 20)
  on.exit(options("scipen" = scipen))

  ## Actual number of samples to save.
  if(!any(grepl("MCMC", method))) {
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
      "optim2" = update_optim2
    )
  }

  ## The proposal function to be used for smooth terms.
  if(!is.function(propose)) {
    propose <- match.arg(propose)
    propose <- switch(propose,
      "iwls" = propose_iwls,
      "rw" = propose_rw,
      "twalk" = propose_twalk,
      "slice" = propose_slice,
      "oslice" = propose_oslice,
      "wslice" = propose_wslice,
      "nadja" = propose_nadja
    )
  }

  ## Function for creating samples when using backfitting.
  if(!is.function(sample)) {
    sample <- match.arg(sample)
    sample <- switch(sample,
      "iwls" = propose_iwls,
      "rw" = propose_rw,
      "twalk" = propose_twalk,
      "slice" = propose_slice,
      "oslice" = propose_oslice,
      "wslice" = propose_wslice,
      "nadja" = propose_nadja
    )
  }
  
  ## Add acceptance rate and fitted values vectors.
  smIWLS <- function(obj, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(obj)
      if(length(unique(nx)) < length(obj)) nx <- 1:length(obj)
      for(j in nx)
        obj[[j]] <- smIWLS(obj[[j]], ...)
    } else {
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          if(!is.null(obj$smooth[[j]]$is.linear)) obj$smooth[[j]]$np <- ncol(obj$smooth[[j]]$X)
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
          obj$smooth[[j]]$adapt <- burnin
          if("backfitting" %in% method) {
            obj$smooth[[j]]$optimize <- TRUE
            obj$smooth[[j]]$criterion <- criterion
          }
        }
      }
    }
    obj
  }

  x <- smIWLS(x, ...)

  ## Number of parameters
  np <- length(nx)

  ## Set up predictors.
  eta <- vector(mode = "list", length = np)
  names(eta) <- nx
  for(j in 1:np)
    eta[[j]] <- rep(0, length(response))

  ## Find starting values with backfitting.
  if(svalues | any(grepl("backfitting", method))) {
    nobs <- if(is.null(dim(response))) length(response) else nrow(response)

    get_edf <- function(x) {
      edf <- 0
      for(j in 1:np) {
        ## And all terms.
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          sedf <- if(is.null(x[[nx[j]]]$smooth[[sj]]$state$edf)) {
            ncol(x[[nx[j]]]$smooth[[sj]]$X)
          } else x[[nx[j]]]$smooth[[sj]]$state$edf
          edf <- edf + sedf
        }
      }
      edf
    }

    fmt <- function(x, width = 8, digits = 2) {
      formatC(round(x, digits), format = "f", digits = digits , width = width)
    }

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
        }
        eps0 <- do.call("cbind", eta)
        eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)

        iter <- iter + 1
      }
      return(list("x" = x, "eta" = eta, "edf" = edf))
    }

    ## 1st backfitting main function.
    backfit <- function(x, eta, verbose = TRUE) {
      edf <- get_edf(x)

      eps0 <- eps + 1; iter <- 1
      while(eps0 > eps & iter < maxit) {
        eta0 <- eta
        ## Cycle through all parameters
        for(j in 1:np) {
          if(outer) {
            peta <- family$map2par(eta)

            ## Compute weights.
            weights <- family$iwls$weights[[nx[j]]](response, peta)

            ## Score.
            score <- family$iwls$score[[nx[j]]](response, peta)

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
            }
          }
        }

        eps0 <- do.call("cbind", eta)
        eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
        if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

        peta <- family$map2par(eta)

        if(any(method == "backfitting") & verbose) {
          IC <- get.ic(family, response, peta, edf, nobs, criterion)

          cat("\r")
          vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
            " logLik ", fmt(family$loglik(response, peta), width = 8, digits = digits),
            " edf ", fmt(edf, width = 6, digits = digits + 2),
            " eps ", fmt(eps0, width = 6, digits = digits + 2),
            " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
          cat(vtxt)

          if(.Platform$OS.type != "unix") flush.console()
        }

        iter <- iter + 1
      }

      IC <- get.ic(family, response, peta, edf, nobs, criterion)

      if(any(method %in% c("backfitting", "backfitting2", "backfitting4")) & verbose) {
        cat("\r")
        vtxt <- paste(criterion, " ", fmt(IC, width = -8, digits = digits),
          " logLik ", fmt(family$loglik(response, peta), width = -9, digits = digits),
          " edf ", fmt(edf, width = -6, digits = digits + 2),
          " eps ", fmt(eps0, width = -6, digits = digits + 2),
          " iteration ", formatC(iter, width = -1 * nchar(maxit)), sep = "")
        cat(vtxt)
        if(.Platform$OS.type != "unix") flush.console()
        cat("\n")
      }

      if(iter == maxit)
        warning("the backfitting algorithm did not converge, please check argument eps and maxit!")

      return(list("x" = x, "eta" = eta, "ic" = IC))
    }

    verbose2 <- if(length(method) < 2) {
      if(method == "MCMC")  FALSE else verbose
    } else verbose

    bf <- backfit(x, eta, verbose = verbose2)
    x <- bf$x; eta <- bf$eta
    rm(bf)

    if(any("backfitting2" %in% method)) {
      tau2 <- NULL
      k <- 1
      for(j in 1:np) {
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          tau2 <- c(tau2, x[[nx[j]]]$smooth[[sj]]$state$tau2)
          k <- k + 1
        }
      }

      objfun <- function(tau2, retbf = FALSE) {
        k <- 1
        for(j in 1:np) {
          for(sj in seq_along(x[[nx[j]]]$smooth)) {
            x[[nx[j]]]$smooth[[sj]]$state$tau2 <- tau2[k]
            k <- k + 1
          }
        }

        bf <- backfit(x, eta, verbose = verbose2)
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
        x[[nx[j]]]$smooth[[sj]]$state$tau2 <- tau2
        bf <- backfit(x, eta, verbose = verbose2)
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
              bf <- backfit(x, eta, verbose = verbose2)
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


  ## Remove some settings.
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      x[[nx[j]]]$smooth[[sj]]$state$iter <- NULL
      x[[nx[j]]]$smooth[[sj]]$state$maxit <- 1
    }
  }

  save.edf <- save.loglik <- NULL
  if(!any(grepl("MCMC", method))) {
    save.edf <- get_edf(x)
    save.loglik <- family$loglik(response, family$map2par(eta))
    if(verbose) cat("generating samples\n")
    ptm <- proc.time()
    for(js in seq_along(iterthin)) {
      for(j in 1:np) {
        ## And all terms.
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          ## Get proposed states.
          p.state <- x[[nx[j]]]$smooth[[sj]]$sample(x[[nx[j]]]$smooth[[sj]], family, response, eta, nx[j], rho = rho)

          accepted <- if(is.na(p.state$alpha)) FALSE else log(runif(1)) <= p.state$alpha

          ## Save the samples.
          x[[nx[j]]]$smooth[[sj]]$s.alpha[js] <- min(c(exp(p.state$alpha), 1), na.rm = TRUE)
          x[[nx[j]]]$smooth[[sj]]$s.accepted[js] <- accepted
          x[[nx[j]]]$smooth[[sj]]$s.samples[js, ] <- unlist(p.state[x[[nx[j]]]$smooth[[sj]]$p.save])
        }
      }
      if(verbose) barfun(ptm, length(iterthin), js, 1, 20, start = FALSE)
    }
    if(verbose) cat("\n")
    deviance <- rep(-2 * family$loglik(response, family$map2par(eta)), length = length(iterthin))
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

  samplesIWLS <- function(obj, id = NULL, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- paste("p", 1:length(obj), sep = "")
      if(length(unique(nx)) < length(obj)) nx <- paste("p", 1:length(obj), sep = "")
      for(j in nx)
        samplesIWLS(obj[[j]], id = j, ...)
    } else {
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          fn <- file.path(tdir, paste(id, if(!is.null(id)) ":", "h",
            obj$hlevel, ":", obj$smooth[[j]]$label, ".raw", sep = ""))
          obj$smooth[[j]]$s.samples <- cbind(obj$smooth[[j]]$s.samples, obj$smooth[[j]]$s.alpha)
          colnames(obj$smooth[[j]]$s.samples) <- c(obj$smooth[[j]]$s.colnames, "alpha")
          write.table(obj$smooth[[j]]$s.samples, file = fn, row.names = FALSE, quote = FALSE)
        }
      }
    }
    NULL
  }

  samplesIWLS(x)

  ## (2) Remove old object.
  rm(x)

  ## (3) Read back in samples
  sf <- grep(".raw", dir(tdir), value = TRUE)
  samples <- NULL
  for(j in sf) {
    st <- as.matrix(read.table(file.path(tdir, j), header = TRUE, colClasses = "numeric"))
    colnames(st) <- paste(gsub(".raw", "", j), gsub("c", "", colnames(st)), sep = ".")
    samples <- cbind(samples, st)
  }
  samples <- if(!any(grepl("MCMC", method))) {
    cbind(samples, "log Lik." = save.loglik, "save.edf" = save.edf)
  } else cbind(samples, "deviance" = deviance)

  return(as.mcmc(samples))
}


resultsIWLS <- function(x, samples)
{
  family <- attr(x, "family")
  grid <- attr(x, "grid")
  if(is.null(grid)) grid <- 100
  if(is.function(family))
    family <- family()

  createIWLSresults <- function(obj, samples, id = NULL)
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
      if(any(grepl("deviance", snames))) {
        DIC <- as.numeric(samples[[j]][, grepl("deviance", snames)])
        pd <- var(DIC, na.rm = TRUE) / 2
        DIC <- mean(DIC, na.rm = TRUE)
      } else {
        DIC <- pd <- NA
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
          k <- obj$smooth[[i]]$np[1]
          pn <- grep(paste(id, "h1", obj$smooth[[i]]$label, sep = ":"), snames,
            value = TRUE, fixed = TRUE) ## FIXME: hlevels!
          pn <- pn[!grepl("tau2", pn) & !grepl("alpha", pn)]
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
            vsamples <- as.numeric(samples[[j]][, tau2])
            vsamples <- vsamples[!nas]
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

          if(is.null(obj$smooth[[i]]$is.linear)) {
            if(!is.list(effects))
              effects <- list()
            if(length(effects)) {
              if(obj$smooth[[i]]$label %in% names(effects)) {
                ct <- gsub(".smooth.spec", "", class(obj$smooth[[i]]))[1]
                if(ct == "random.effect") ct <- "re"
                obj$smooth[[i]]$label <- paste(obj$smooth[[i]]$label, ct, sep = ":")
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

      class(rval[[j]]) <- "bayesr"
    }
    names(rval) <- paste("Chain", 1:chains, sep = "_")
    if(length(rval) < 2) {
      rval <- rval[[1]]
    }
    class(rval) <- "bayesr"
    return(rval)
  }

  if(inherits(x, "bayesr.input") & !all(c("formula", "fake.formula", "response") %in% names(x))) {
    nx <- names(x)
    nx <- nx[nx != "call"]
    rval <- list()
    fn <- family$names
    cat <- if(!is.null(family$cat)) family$cat else FALSE
    if(cat)
      fn <- gsub(attr(attr(x, "model.frame"), "response.name"), "", names(x))
    if(length(fn) != length(nx))
      fn <- paste(fn, 1:length(nx), sep = "")
    for(j in seq_along(nx)) {
      rval[[nx[j]]] <- createIWLSresults(x[[nx[j]]], samples, id = fn[j])
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
      rval <- rval[!grepl(reference, fn)]
    }
    attr(rval, "family") <- family
    class(rval) <- "bayesr"
    return(rval)
  } else {
    return(createIWLSresults(x, samples))
  }
}

