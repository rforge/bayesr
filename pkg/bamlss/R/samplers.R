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


gmcmc <- function(fun, theta, priors = NULL,
  propose = NULL, data = NULL, n.iter = 12000,
  burnin = 2000, thin = 10, verbose = TRUE, step = 20, ...)
{
  require("coda")

  if(!is.list(theta))
    theta <- list(theta)
  if(is.null(names(theta)))
    names(theta) <- paste("theta[", 1:length(theta), "]", sep = "")
  ntheta <- names(theta)
  k <- length(theta)
  for(j in 1:k) {
    if(is.null(names(theta[[j]])))
      names(theta[[j]]) <- paste("p[", 1:length(theta[[j]]), "]", sep = "")
  }

  if(is.null(propose))
    propose <- gmcmc_propose_default
  if(!is.list(propose))
    propose <- list(propose)
  propose <- rep(propose, length.out = k)
  if(is.null(names(propose)))
    names(propose) <- names(theta)
  if(!all(names(propose) %in% names(theta)))
    stop("the 'propose' list() names are different from theta!")

  if(!is.null(data)) {
    data <- as.list(data)
    if(length(data) != k)
      stop("the data list() must have the same length as theta!")
    if(is.null(names(data)))
      names(data) <- names(theta)
    if(!all(names(data) %in% names(theta)))
      stop("the 'data' list() names are different from theta!")
  }

  if(!is.null(priors)) {
    priors <- as.list(priors)
    if(length(priors) != k)
      stop("the data list() must have the same length as theta!")
    if(is.null(names(priors)))
      names(priors) <- names(theta)
    if(!all(names(priors) %in% names(theta)))
      stop("the 'data' list() names are different from theta!")
    for(j in ntheta) {
      if(!is.function(priors[[j]])) {
        stop(paste("the prior for block ", j,
          " must return a function!", sep = ""))
      }
    }
  }

  iterthin <- as.integer(seq(burnin, n.iter, by = thin))

  theta.save <- vector(mode = "list", length = length(theta))
  names(theta.save) <- names(theta)
  for(j in ntheta) {
    p0 <- propose[[j]](fun, theta, id = j,
      prior = priors[[j]], data = data[[j]],
      iteration = 1, n.iter = n.iter, ...)
    if(!is.list(p0)) {
      stop(paste("the propose() function for block ", j,
        " must return a named list()!", sep = ""))
    }
    if(is.null(p0$alpha)) {
      stop(paste("the propose() function for block ", j,
        " must return the acceptance probability 'alpha'!", sep = ""))
    }
    if(is.null(p0$parameters)) {
      stop(paste("the propose() function for block ", j,
        " must return a vector 'parameters'!", sep = ""))
    }
    if(length(p0$parameters) != length(theta[[j]])) {
      stop(paste("the propose() function for block ", j,
        " must return a vector 'parameters' with the same length of initial parameters in theta!",
        sep = ""))
    }
    psave <- if(is.null(p0$save)) p0$parameters else p0$save
    if(is.null(names(psave)))
      names(psave) <- paste("p[", 1:length(psave), "]", sep = "")
    theta.save[[j]] <- list(
      "samples" = matrix(NA, nrow = length(iterthin), ncol = length(psave)),
      "alpha" = rep(NA, length = length(iterthin)),
      "accepted" = rep(NA, length = length(iterthin))
    )
    colnames(theta.save[[j]]$samples) <- names(psave)
  }

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
      cat(p, rt, sep = " ")
      if(.Platform$OS.type != "unix") flush.console()
    }
  }

  ptm <- proc.time()

  for(i in 1:n.iter) {
    if(save <- i %in% iterthin)
      js <- which(iterthin == i)
    for(j in ntheta) {
      ## Get proposed states.
      state <- propose[[j]](fun, theta, id = j,
        prior = priors[[j]], data = data[[j]],
        iteration = i, n.iter = n.iter, ...)

      ## If accepted, set current state to proposed state.
      accepted <- if(is.na(state$alpha)) FALSE else log(runif(1)) <= state$alpha

      if(accepted)
        theta[[j]] <- state$parameters

      ## Save the samples.
      if(save) {
        theta.save[[j]]$samples[js, ] <- unlist(if(is.null(state$save)) state$parameters else state$save)
        theta.save[[j]]$alpha[js] <- min(c(exp(state$alpha), 1), na.rm = TRUE)
        theta.save[[j]]$accepted[js] <- accepted
      }
    }

    if(verbose) barfun(ptm, n.iter, i, step, nstep)
  }

  if(verbose) cat("\n")

  theta.save <- lapply(theta.save, function(x) { do.call("cbind", x) })
  if(length(theta.save) > 1L) {
    for(j in seq_along(theta.save))
      colnames(theta.save[[j]]) <- paste(ntheta[j], colnames(theta.save[[j]]), sep = ".")
    theta.save <- do.call("cbind", theta.save)
  } else {
    theta.save <- theta.save[[1L]]
  }
  theta.save <- mcmc(theta.save, start = burnin, end = n.iter, thin = thin)

  return(theta.save)
}


gmcmc_propose_default <- function(fun, theta, id, prior, ...)
{
  require("mvtnorm")

  args <- list(...)
  iteration <- args$iteration

  if(is.null(attr(theta[[id]], "scale")))
    attr(theta[[id]], "scale") <- 1
  if(is.null(attr(theta[[id]], "P")))
    attr(theta[[id]], "sigma") <- diag(length(theta[[id]]))

  if(iteration < floor(0.15 * args$n.iter)) {
    scale <- attr(theta[[id]], "scale")
    sigma <- attr(theta[[id]], "sigma")

    k <- 1
    do <- TRUE
    while(do & k < 100) {
      ll0 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
      p0 <- if(is.null(prior)) {
        sum(dnorm(theta[[id]], sd = 1000, log = TRUE), na.rm = TRUE)
      } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

      theta[[id]] <- drop(rmvnorm(n = 1, mean = theta[[id]], sigma = scale * sigma))

      ll1 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
      p1 <- if(is.null(prior)) {
        sum(dnorm(theta[[id]], sd = 1000, log = TRUE), na.rm = TRUE)
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

    attr(theta[[id]], "scale") <- scale
    attr(theta[[id]], "sigma") <- sigma
  }

  scale <- attr(theta[[id]], "scale")
  sigma <- attr(theta[[id]], "sigma")

  ll0 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p0 <- if(is.null(prior)) {
    sum(dnorm(theta[[id]], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  theta[[id]] <- drop(rmvnorm(n = 1, mean = theta[[id]], sigma = scale * sigma))

  ll1 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p1 <- if(is.null(prior)) {
    sum(dnorm(theta[[id]], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  alpha <- drop((ll1 + p1) - (ll0 + p0))

  rval <- list("parameters" = theta[[id]], "alpha" = alpha)
  attr(rval$parameters, "scale") <- scale
  attr(rval$parameters, "sigma") <- sigma

  rval
}


gmcmc_unislice <- function(fun, theta, id, prior, j, ...,
  w = 1, m = 30, lower = -Inf, upper = +Inf)
{
  args <- list(...)

  x0 <- theta[[id]][j]
  gL <- gR <- theta

  fun2 <- function(theta) {
    ll <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
    lp <- if(is.null(prior)) {
      sum(dnorm(theta[[id]], sd = 1000, log = TRUE), na.rm = TRUE)
    } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)
    return(ll + lp)
  }

  gx0 <- fun2(theta)

  ## Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)

  ## Find the initial interval to sample from.
  ## w <- w * abs(x0) ## FIXME???
  u <- runif(1, 0, w)
  gL[[id]][j] <- theta[[id]][j] - u
  gR[[id]][j] <- theta[[id]][j] + (w - u)  ## should guarantee that g[j] is in [L, R], even with roundoff

  ## Expand the interval until its ends are outside the slice, or until
  ## the limit on steps is reached.
  if(is.infinite(m)) {
    repeat {
      if(gL[[id]][j] <= lower) break
      if(fun2(gL) <= logy) break
      gL[[id]][j] <- gL[[id]][j] - w
    }
    repeat {
      if(gR[[id]][j] >= upper) break
      if(fun2(gR) <= logy) break
      gR[[id]][j] <- gR[[id]][j] + w
    }
  } else {
    if(m > 1) {
      J <- floor(runif(1, 0, m))
      K <- (m - 1) - J
      while(J > 0) {
        if(gL[[id]][j] <= lower) break
        if(fun2(gL) <= logy) break
        gL[[id]][j] <- gL[[id]][j] - w
        J <- J - 1
      }
      while(K > 0) {
        if(gR[[id]][j] >= upper) break
        if(fun2(gR) <= logy) break
        gR[[id]][j] <- gR[[id]][j] + w
        K <- K - 1
      }
    }
  }

  ## Shrink interval to lower and upper bounds.
  if(gL[[id]][j] < lower) {
    gL[[id]][j] <- lower
  }
  if(gR[[id]][j] > upper) {
    gR[[id]][j] <- upper
  }

  ## Sample from the interval, shrinking it on each rejection.
  repeat {
    theta[[id]][j] <- runif(1, gL[[id]][j], gR[[id]][j])

    gx1 <- fun2(theta)

    if(gx1 >= logy) break

    if(theta[[id]][j] > x0) {
      gR[[id]][j] <- theta[[id]][j]
    } else {
      gL[[id]][j] <- theta[[id]][j]
    }
  }

  ## Return the point sampled
  return(theta[[id]][j])
}


gmcmc_slice <- function(fun, theta, id, prior, ...)
{
  for(j in seq_along(theta[[id]]))
    theta[[id]][j] <- gmcmc_unislice(fun, theta, id, prior, j, ...)
  return(list("parameters" = theta[[id]], "alpha" = log(1)))
}


gmcmc_fs <- function(fun, theta, id, prior, ...)
{
  args <- list(...)
  theta[[id]] <- theta[[id]] + args$I[[id]](theta[[id]]) %*% args$S[[id]](theta[[id]])
}

