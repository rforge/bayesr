MCMCpack <- function(x, n.iter = 1200, burnin = 200, thin = 1, verbose = 100, ...)
{
  require("MCMCpack")

  par <- make_par(x, type = 2)

  post.samp <- MCMCmetrop1R(log_posterior, theta.init = par$par,
    x = x, logfun = TRUE, V = attr(x, "hessian"),
    mcmc = n.iter, burnin = burnin, thin = thin, verbose = verbose)

  colnames(post.samp) <- names(par$par)

  post.samp
}


GMCMC <- function(x, n.iter = 1200, burnin = 200, thin = 1, verbose = 100,
  propose = gmcmc_iwls, ...)
{
  family <- attr(x, "family")
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  np <- length(nx)
  response <- attr(x, "response.vec")

  theta <- smooths <- propose2 <- fitfun <- list()
  for(i in names(x)) {
    theta[[i]] <- smooths[[i]] <- propose2[[i]] <- fitfun[[i]] <- list()
    nt <- NULL
    for(j in seq_along(x[[i]]$smooth)) {
      theta[[i]][[j]] <- x[[i]]$smooth[[j]]$state$parameters
      nt <- c(nt, x[[i]]$smooth[[j]]$label)
      smooths[[i]][[j]] <- x[[i]]$smooth[[j]]
      propose2[[i]][[j]] <- if(!is.null(propose)) propose else x[[i]]$smooth[[j]]$propose
      fitfun[[i]][[j]] <- function(x, p) {
        x$get.mu(x$X, p)
      }
    }
    names(theta[[i]]) <- names(smooths[[i]]) <- names(propose2[[i]]) <- names(fitfun[[i]]) <- nt
  }
  rm(x)

  logLik <- function(eta) {
    family$loglik(response, family$map2par(eta))
  }

  samps <- gmcmc(fun = family, theta = theta, fitfun = fitfun, data = smooths,
    propose = propose2, logLik = logLik, n.iter = n.iter, burnin = burnin,
    response = response, simplify = FALSE, ...)

  samps
}


gmcmc <- function(fun, theta, priors = NULL, propose = NULL,
  fitfun = NULL, logLik = NULL, data = NULL,
  n.iter = 12000, burnin = 2000, thin = 10, verbose = TRUE, step = 20,
  simplify = TRUE,  ...)
{
  require("coda")

  if(!is.list(theta)) {
    theta <- list(theta)
    names(theta) <- names(formals(fun))[1]
  }
  if(is.null(names(theta)))
    names(theta) <- paste("theta[", 1:length(theta), "]", sep = "")
  ntheta <- names(theta)
  k <- length(theta)
  for(i in 1:k) {
    if(!is.list(theta[[i]]))
      theta[[i]] <- list(theta[[i]])
    if(is.null(names(theta[[i]])))
      names(theta[[i]]) <- paste("term[", 1:length(theta[[i]]), "]", sep = "")
    for(j in seq_along(theta[[i]])) {
      if(is.null(names(theta[[i]][[j]])))
        names(theta[[i]][[j]]) <- paste("[", 1:length(theta[[i]][[j]]), "]", sep = "")
    }
  }

  if(is.null(propose))
    propose <- gmcmc_mvnorm
  if(!is.list(propose))
    propose <- rep(list(propose), length.out = k)
  if(is.null(names(propose)))
    names(propose) <- names(theta)
  if(!all(names(propose) %in% names(theta)))
    stop("the 'propose' list() names are different from theta!")
  for(i in ntheta) {
    if(!is.list(propose[[i]]))
      propose[[i]] <- rep(list(propose[[i]]), length.out = length(theta[[i]]))
    if(is.null(names(propose[[i]])))
      names(propose[[i]]) <- names(theta[[i]])
    if(!all(names(propose[[i]]) %in% names(theta[[i]])))
      stop(paste("propose list() names for parameter", i, "different from the theta list() names!"))
  }

  if(is.null(fitfun)) {
    fitfun <- function(x, p) {
      if(is.null(x)) {
        return(sum(p))
      } else return(drop(x %*% p))
    }
  }
  if(!is.list(fitfun))
    fitfun <- rep(list(fitfun), length.out = k)
  if(is.null(names(fitfun)))
    names(fitfun) <- names(theta)
  if(!all(names(fitfun) %in% names(theta)))
    stop("the 'fitfun' list() names are different from theta!")
  for(i in ntheta) {
    if(!is.list(fitfun[[i]]))
      fitfun[[i]] <- rep(list(fitfun[[i]]), length.out = length(theta[[i]]))
    if(is.null(names(fitfun[[i]])))
      names(fitfun[[i]]) <- names(theta[[i]])
    if(!all(names(fitfun[[i]]) %in% names(theta[[i]])))
      stop(paste("fitfun list() names for parameter", i, "different from the theta list() names!"))
    for(j in names(fitfun[[i]])) {
      if(!is.function(fitfun[[i]][[j]])) {
        stop(paste("the prior for block [", i, "][", j,
          "] must return a function!", sep = ""))
      }
    }
  }

  if(!is.null(data)) {
    if(!is.list(data))
      data <- list(data)
    if(is.null(names(data)))
      names(data) <- names(theta)
    if(!all(names(data) %in% names(theta)))
      stop("the 'data' list() names are different from theta!")
    for(i in ntheta) {
      if(!is.list(data[[i]]))
        data[[i]] <- rep(list(data[[i]]), length.out = length(theta[[i]]))
      if(is.null(names(data[[i]])))
        names(data[[i]]) <- names(theta[[i]])
      if(!all(names(data[[i]]) %in% names(theta[[i]])))
        stop(paste("data list() names for parameter", i, "different from the theta list() names!"))
    }
  }

  if(!is.null(priors)) {
    if(!is.list(priors))
      priors <- rep(list(priors), length.out = k)
    if(is.null(names(priors)))
      names(priors) <- names(theta)
    if(!all(names(priors) %in% names(theta)))
      stop("the 'priors' list() names are different from theta!")
    for(i in ntheta) {
      if(!is.list(priors[[i]]))
        priors[[i]] <- rep(list(priors[[i]]), length.out = length(theta[[i]]))
      if(!is.null(names(priors)))
        names(priors[[i]]) <- names(theta[[i]])
      if(!all(names(priors[[i]]) %in% names(theta[[i]])))
        stop(paste("priors list() names for parameter", i, "different from the theta list() names!"))
      for(j in names(priors[[i]])) {
        if(!is.function(priors[[i]][[j]])) {
          stop(paste("the prior for block [", i, "][", j,
            "] must return a function!", sep = ""))
        }
      }
    }
  }

  if(burnin < 1) burnin <- 1
  iterthin <- as.integer(seq(burnin, n.iter, by = thin))

  theta.save <- vector(mode = "list", length = length(theta))
  names(theta.save) <- names(theta)
  eta <- theta.save
  for(i in ntheta) {
    eta[[i]] <- 0
    for(j in names(theta[[i]]))
      eta[[i]] <- fitfun[[i]][[j]](data[[i]][[j]], theta[[i]][[j]])
  }
  for(i in ntheta) {
    theta.save[[i]] <- vector(mode = "list", length = length(theta[[i]]))
    names(theta.save[[i]]) <- names(theta[[i]])
    for(j in names(theta[[i]])) {
      p0 <- propose[[i]][[j]](fun, theta, id = c(i, j),
        prior = priors[[i]][[j]], data = data[[i]][[j]], eta = eta,
        iteration = 1, n.iter = n.iter, ...)
      if(!is.list(p0)) {
        stop(paste("the propose() function for block [", i, "][", j,
          "] must return a named list()!", sep = ""))
      }
      if(is.null(p0$alpha)) {
        stop(paste("the propose() function for block [", i, "][", j,
          "] must return the acceptance probability 'alpha'!", sep = ""))
      }
      if(is.null(p0$parameters)) {
        stop(paste("the propose() function for block [", i, "][", j,
          "] must return a vector 'parameters'!", sep = ""))
      }
      if(length(p0$parameters) != length(theta[[i]][[j]])) {
        stop(paste("the propose() function for block [", i, "][", j,
          "] must return a vector 'parameters' with the same length of initial parameters in theta!",
          sep = ""))
      }
      if(is.null(names(p0$parameters)))
        names(p0$parameters) <- paste("[", 1:length(p0$parameters), "]", sep = "")
      if(length(p0$parameters) < 2)
        names(p0$parameters) <- NULL
      theta.save[[i]][[j]] <- list(
        "samples" = matrix(NA, nrow = length(iterthin), ncol = length(p0$parameters)),
        "alpha" = rep(NA, length = length(iterthin)),
        "accepted" = rep(NA, length = length(iterthin))
      )
      colnames(theta.save[[i]][[j]]$samples) <- names(p0$parameters)
    }
  }
  ll <- rep(NA, length = length(iterthin))

  nstep <- step
  step <- floor(n.iter / step)

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
      elapsed <- if(elapsed > 60) {
        paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
      } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
      cat(p, rt, elapsed, sep = " ")
      if(.Platform$OS.type != "unix") flush.console()
    }
  }

  ptm <- proc.time()
  for(iter in 1:n.iter) {
    if(save <- iter %in% iterthin)
      js <- which(iterthin == iter)
    for(i in ntheta) {
      for(j in names(theta[[i]])) {
        ## Get proposed states.
        state <- propose[[i]][[j]](fun, theta, id = c(i, j),
          prior = priors[[i]][[j]], data = data[[i]][[j]], eta = eta,
          iteration = iter, n.iter = n.iter, ...)

        ## If accepted, set current state to proposed state.
        accepted <- if(is.na(state$alpha)) FALSE else log(runif(1)) <= state$alpha

        if(accepted) {
          eta[[i]] <- eta[[i]] - fitfun[[i]][[j]](data[[i]][[j]], theta[[i]][[j]])
          theta[[i]][[j]] <- state$parameters
          eta[[i]] <- eta[[i]] + fitfun[[i]][[j]](data[[i]][[j]], theta[[i]][[j]])
        }

        ## Save the samples.
        if(save) {
          theta.save[[i]][[j]]$samples[js, ] <- theta[[i]][[j]]
          theta.save[[i]][[j]]$alpha[js] <- min(c(exp(state$alpha), 1), na.rm = TRUE)
          theta.save[[i]][[j]]$accepted[js] <- accepted
          ll[js] <- if(!is.null(logLik)) {
            logLik(eta)
          } else sum(do.call(fun, c(theta, list(...))[names(formals(fun))]), na.rm = TRUE)
        }
      }
    }

    if(verbose) barfun(ptm, n.iter, iter, step, nstep)
  }

  if(verbose) cat("\n")

  for(i in ntheta) {
    theta.save[[i]] <- lapply(theta.save[[i]], function(x) { do.call("cbind", x) })
    for(j in names(theta.save[[i]])) {
      cn <- i
      if(length(theta.save[[i]]) > 1 | !simplify)
        cn <- paste(cn, j, sep = ".")
      cn2 <- colnames(theta.save[[i]][[j]])
      for(k in seq_along(cn2)) {
        cn2[k] <- if(grepl("[", cn2[k], fixed = TRUE)) {
          paste(if(strsplit(cn2[k], "[", fixed = TRUE)[[1]][1] != "") "." else NULL, cn2[k], sep = "")
        } else paste(if(cn2[k] != "") "." else NULL, cn2[k], sep = "")
      }
      cn <- paste(cn, cn2, sep = "")
      colnames(theta.save[[i]][[j]]) <- cn
    }
    theta.save[[i]] <- if(length(theta.save[[i]]) < 2 & simplify) {
      theta.save[[i]][[1]]
    } else do.call("cbind", theta.save[[i]])
  }
  theta.save <- if(length(theta.save) < 2) theta.save[[1]] else do.call("cbind", theta.save)
  theta.save <- cbind(theta.save, "logLik" = ll)

  theta.save <- mcmc(theta.save, start = burnin, end = n.iter, thin = thin)
  class(theta.save) <- c("gmcmc", class(theta.save))

  return(theta.save)
}


DIC <- function(object, ...)
{
  UseMethod("DIC")
}

DIC.gmcmc <- function(object, ...)
{
  object <- list(object, ...)
  rval <- NULL
  for(j in seq_along(object)) {
    dev <- -2 * object[[j]][, "logLik"]
    pd <- var(dev, na.rm = TRUE) / 2
    rval <- rbind(rval, c("DIC" = mean(dev, na.rm = TRUE), "pd" = pd))
  }
  Call <- match.call()
  row.names(rval) <- if(nrow(rval) > 1) as.character(Call[-1L]) else ""
  rval
}

plot.gmcmc <- function(x, ...)
{
  x <- na.omit(x)
  class(x) <- "mcmc"
  plot(x, ...)
}


gmcmc_mvnorm <- function(fun, theta, id, prior, ...)
{
  require("mvtnorm")

  args <- list(...)
  iteration <- args$iteration

  if(is.null(attr(theta[[id[1]]][[id[2]]], "scale")))
    attr(theta[[id[1]]][[id[2]]], "scale") <- 1
  if(is.null(attr(theta[[id[1]]][[id[2]]], "P")))
    attr(theta[[id[1]]][[id[2]]], "sigma") <- diag(length(theta[[id[1]]][[id[2]]]))

  if(iteration < floor(0.15 * args$n.iter)) {
    scale <- attr(theta[[id[1]]][[id[2]]], "scale")
    sigma <- attr(theta[[id[1]]][[id[2]]], "sigma")

    k <- 1
    do <- TRUE
    while(do & k < 100) {
      ll0 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
      p0 <- if(is.null(prior)) {
        sum(dnorm(theta[[id[1]]][[id[2]]], sd = 1000, log = TRUE), na.rm = TRUE)
      } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

      theta[[id[1]]][[id[2]]] <- drop(rmvnorm(n = 1, mean = theta[[id[1]]][[id[2]]], sigma = scale * sigma))

      ll1 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
      p1 <- if(is.null(prior)) {
        sum(dnorm(theta[[id[1]]][[id[2]]], sd = 1000, log = TRUE), na.rm = TRUE)
      } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

      alpha <- drop((ll1 + p1) - (ll0 + p0))

      if(!is.finite(alpha)) next

      accepted <- if(is.na(alpha)) FALSE else log(runif(1)) <= alpha

      scale <- if(alpha < log(0.23)) {
        scale - 0.1 * scale
      } else {
        scale + 0.1 * scale
      }

      if(accepted) do <- FALSE

      k <- k + 1
    }

    attr(theta[[id[1]]][[id[2]]], "scale") <- scale
    attr(theta[[id[1]]][[id[2]]], "sigma") <- sigma
  }

  scale <- attr(theta[[id[1]]][[id[2]]], "scale")
  sigma <- attr(theta[[id[1]]][[id[2]]], "sigma")

  ll0 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p0 <- if(is.null(prior)) {
    sum(dnorm(theta[[id[1]]][[id[2]]], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  theta[[id[1]]][[id[2]]] <- drop(rmvnorm(n = 1, mean = theta[[id[1]]][[id[2]]], sigma = scale * sigma))

  ll1 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p1 <- if(is.null(prior)) {
    sum(dnorm(theta[[id[1]]][[id[2]]], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  alpha <- drop((ll1 + p1) - (ll0 + p0))

  rval <- list("parameters" = theta[[id[1]]][[id[2]]], "alpha" = alpha)

  attr(rval$parameters, "scale") <- scale
  attr(rval$parameters, "sigma") <- sigma

  rval
}


gmcmc_unislice <- function(fun, theta, id, prior, j, ...,
  w = 1, m = 30, lower = -Inf, upper = +Inf)
{
  args <- list(...)

  x0 <- theta[[id[1]]][[id[2]]][j]
  gL <- gR <- theta

  fun2 <- function(theta) {
    ll <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
    lp <- if(is.null(prior)) {
      sum(dnorm(theta[[id[1]]][[id[2]]], sd = 1000, log = TRUE), na.rm = TRUE)
    } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)
    return(ll + lp)
  }

  gx0 <- fun2(theta)

  ## Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)

  ## Find the initial interval to sample from.
  ## w <- w * abs(x0) ## FIXME???
  u <- runif(1, 0, w)
  gL[[id[1]]][[id[2]]][j] <- theta[[id[1]]][[id[2]]][j] - u
  gR[[id[1]]][[id[2]]][j] <- theta[[id[1]]][[id[2]]][j] + (w - u)  ## should guarantee that g[j] is in [L, R], even with roundoff

  ## Expand the interval until its ends are outside the slice, or until
  ## the limit on steps is reached.
  if(is.infinite(m)) {
    repeat {
      if(gL[[id[1]]][[id[2]]][j] <= lower) break
      if(fun2(gL) <= logy) break
      gL[[id[1]]][[id[2]]][j] <- gL[[id[1]]][[id[2]]][j] - w
    }
    repeat {
      if(gR[[id[1]]][[id[2]]][j] >= upper) break
      if(fun2(gR) <= logy) break
      gR[[id[1]]][[id[2]]][j] <- gR[[id[1]]][[id[2]]][j] + w
    }
  } else {
    if(m > 1) {
      J <- floor(runif(1, 0, m))
      K <- (m - 1) - J
      while(J > 0) {
        if(gL[[id[1]]][[id[2]]][j] <= lower) break
        if(fun2(gL) <= logy) break
        gL[[id[1]]][[id[2]]][j] <- gL[[id[1]]][[id[2]]][j] - w
        J <- J - 1
      }
      while(K > 0) {
        if(gR[[id[1]]][[id[2]]][j] >= upper) break
        if(fun2(gR) <= logy) break
        gR[[id[1]]][[id[2]]][j] <- gR[[id[1]]][[id[2]]][j] + w
        K <- K - 1
      }
    }
  }

  ## Shrink interval to lower and upper bounds.
  if(gL[[id[1]]][[id[2]]][j] < lower) {
    gL[[id[1]]][[id[2]]][j] <- lower
  }
  if(gR[[id[1]]][[id[2]]][j] > upper) {
    gR[[id[1]]][[id[2]]][j] <- upper
  }

  ## Sample from the interval, shrinking it on each rejection.
  repeat {
    theta[[id[1]]][[id[2]]][j] <- runif(1, gL[[id[1]]][[id[2]]][j], gR[[id[1]]][[id[2]]][j])

    gx1 <- fun2(theta)

    if(gx1 >= logy) break

    if(theta[[id[1]]][[id[2]]][j] > x0) {
      gR[[id[1]]][[id[2]]][j] <- theta[[id[1]]][[id[2]]][j]
    } else {
      gL[[id[1]]][[id[2]]][j] <- theta[[id[1]]][[id[2]]][j]
    }
  }

  ## Return the point sampled
  return(theta[[id[1]]][[id[2]]][j])
}


gmcmc_slice <- function(fun, theta, id, prior, ...)
{
  for(j in seq_along(theta[[id[1]]][[id[2]]]))
    theta[[id[1]]][[id[2]]][j] <- gmcmc_unislice(fun, theta, id, prior, j, ...)
  return(list("parameters" = theta[[id[1]]][[id[2]]], "alpha" = log(1)))
}


gmcmc_iwls <- function(family, theta, id, prior, eta, response, data, ...)
{
  require("mvtnorm")

  theta <- theta[[id[1]]][[id[2]]]
  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$get.mu(data$X, theta)

  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute weights.
  weights <- family$weights[[id[1]]](response, peta)

  ## Score.
  score <- family$score[[id[1]]](response, peta)

  ## Compute working observations.
  z <- eta[[id[1]]] + 1 / weights * score

  ## Compute old log likelihood and old log coefficients prior.
  pibeta <- family$loglik(response, peta)
  p1 <- data$prior(theta)

  ## Compute partial predictor.
  eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")

  ## Compute mean and precision.
  S <- 0
  XW <- t(data$X * weights)
  P <- if(data$fixed) {
    if(k <- ncol(data$X) < 2) {
      1 / (XW %*% data$X)
    } else chol2inv(chol(XW %*% data$X))
  } else {
    for(j in seq_along(data$S))
      S <- S + 1 / get.par(theta, "tau2")[j] * data$S[[j]]
    chol2inv(chol(XW %*% data$X + S))
  }
  P[P == Inf] <- 0
  M <- P %*% (XW %*% (z - eta[[id[1]]]))

  ## Save old coefficients
  g0 <- drop(get.par(theta, "gamma"))

  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = M, sigma = P))

  ## Compute log priors.
  p2 <- data$prior(c("g" = g, get.par(theta, "tau2")))
  qbetaprop <- dmvnorm(g, mean = M, sigma = P, log = TRUE)

  ## Compute fitted values.        
  attr(theta, "fitted.values") <- data$get.mu(data$X, g)

  ## Set up new predictor.
  eta[[id[1]]] <- eta[[id[1]]] + attr(theta, "fitted.values")

  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute new log likelihood.
  pibetaprop <- family$loglik(response, peta)

  ## Compute new weights
  weights <- family$weights[[id[1]]](response, peta)

  ## New score.
  score <- family$score[[id[1]]](response, peta)

  ## New working observations.
  z <- eta[[id[1]]] + 1 / weights * score

  ## Compute mean and precision.
  XW <- t(data$X * weights)
  P2 <- if(data$fixed) {
    if(k < 2) {
      1 / (XW %*% data$X)
    } else chol2inv(chol(XW %*% data$X))
  } else {
    chol2inv(L <- chol(P0 <- XW %*% data$X + S))
  }
  P2[P2 == Inf] <- 0
  M2 <- P2 %*% (XW %*% (z - (eta[[id[1]]] - attr(theta, "fitted.values"))))

  ## Get the log prior.
  qbeta <- dmvnorm(g0, mean = M2, sigma = P2, log = TRUE)

  ## Sample variance parameter.
  if(!data$fixed & is.null(data$sp)) {
    if(!data$fixed & is.null(data$sp)) {
      tau2 <- NULL
      for(j in seq_along(data$S)) {
        a <- data$rank[j] / 2 + data$a
        b <- 0.5 * crossprod(g, data$S[[j]]) %*% g + data$b
        tau2 <- c(tau2, 1 / rgamma(1, a, b))
      }
      theta <- set.par(theta, tau2, "tau2")
    }
  }
  theta <- set.par(theta, g, "gamma")

  ## Compute acceptance probablity.
  alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

  return(list("parameters" = theta, "alpha" = alpha))
}

