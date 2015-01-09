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

  theta <- smooths <- propose2 <- list()
  for(i in names(x)) {
    theta[[i]] <- smooths[[i]] <- propose2[[i]] <- list()
    nt <- NULL
    for(j in seq_along(x[[i]]$smooth)) {
      theta[[i]][[j]] <- x[[i]]$smooth[[j]]$state$parameters
      nt <- c(nt, x[[i]]$smooth[[j]]$label)
      smooths[[i]][[j]] <- x[[i]]$smooth[[j]]
      propose2[[i]][[j]] <- if(!is.null(propose)) propose else x[[i]]$smooth[[j]]$propose
    }
    names(theta[[i]]) <- names(smooths[[i]]) <- names(propose2[[i]]) <- nt
  }
  rm(x)

  samps <- gmcmc(fun = family, theta = theta, data = smooths,
    propose = propose2, n.iter = n.iter, burnin = burnin, ...)
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
  for(i in 1:k) {
    if(!is.list(theta[[i]]))
      theta[[i]] <- list(theta[[i]])
    if(is.null(names(theta[[i]])))
      names(theta[[i]]) <- paste("term[", 1:length(theta[[i]]), "]", sep = "")
    for(j in seq_along(theta[[i]])) {
      if(is.null(names(theta[[i]][[j]])))
        names(theta[[i]][[j]]) <- paste("p[", 1:length(theta[[i]][[j]]), "]", sep = "")
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

  iterthin <- as.integer(seq(burnin, n.iter, by = thin))

  theta.save <- vector(mode = "list", length = length(theta))
  names(theta.save) <- names(theta)
  for(i in ntheta) {
    theta.save[[i]] <- vector(mode = "list", length = length(theta[[i]]))
    names(theta.save[[i]]) <- names(theta[[i]])
    for(j in names(theta[[i]])) {
      p0 <- propose[[i]][[j]](fun, theta, id = c(i, j),
        prior = priors[[i]][[j]], data = data[[i]][[j]],
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
      psave <- if(is.null(p0$save)) p0$parameters else p0$save
      if(is.null(names(psave)))
        names(psave) <- paste("p[", 1:length(psave), "]", sep = "")
      theta.save[[i]][[j]] <- list(
        "samples" = matrix(NA, nrow = length(iterthin), ncol = length(psave)),
        "alpha" = rep(NA, length = length(iterthin)),
        "accepted" = rep(NA, length = length(iterthin))
      )
      colnames(theta.save[[i]][[j]]$samples) <- names(psave)
    }
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

  for(iter in 1:n.iter) {
    if(save <- iter %in% iterthin)
      js <- which(iterthin == iter)
    for(i in ntheta) {
      for(j in names(theta[[i]])) {
        ## Get proposed states.
        state <- propose[[i]][[j]](fun, theta, id = c(i, j),
          prior = priors[[i]][[j]], data = data[[i]][[j]],
          iteration = iter, n.iter = n.iter, ...)

        ## If accepted, set current state to proposed state.
        accepted <- if(is.na(state$alpha)) FALSE else log(runif(1)) <= state$alpha

        if(accepted)
          theta[[i]][[j]] <- state$parameters

        ## Save the samples.
        if(save) {
          theta.save[[i]][[j]]$samples[js, ] <- unlist(if(is.null(state$save)) state$parameters else state$save)
          theta.save[[i]][[j]]$alpha[js] <- min(c(exp(state$alpha), 1), na.rm = TRUE)
          theta.save[[i]][[j]]$accepted[js] <- accepted
        }
      }
    }

    if(verbose) barfun(ptm, n.iter, iter, step, nstep)
  }

  if(verbose) cat("\n")

  for(i in ntheta) {
    theta.save[[i]] <- lapply(theta.save[[i]], function(x) { do.call("cbind", x) })
    for(j in names(theta.save[[i]]))
      colnames(theta.save[[i]][[j]]) <- paste(i, j, colnames(theta.save[[i]][[j]]), sep = ".")
    theta.save[[i]] <- do.call("cbind", theta.save[[i]])
  }
  theta.save <- do.call("cbind", theta.save)

  theta.save <- mcmc(theta.save, start = burnin, end = n.iter, thin = thin)

  return(theta.save)
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


gmcmc_iwls <- function(fun, theta, id, prior, ...)
{
  args <- list(...)
print(id)
stop()
  propose_iwls0(x = args$data, family = fun, response, eta, id = id[1], ...)
}

