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
  propose = "iwls", cores = NULL, chains = NULL, ...)
{
  family <- attr(x, "family")
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  np <- length(nx)
  response <- attr(x, "response.vec")

  if(is.character(propose)) {
    propose <- if(grepl("gmcmc_", propose)) {
      eval(parse(text = propose))
    } else {
      if(grepl("sm.", propose))
        eval(parse(text = paste("gmcmc", propose, sep = "_")))
      else
        eval(parse(text = paste("gmcmc", paste("sm", propose, sep = "."), sep = "_")))
    }
  }

  theta <- smooths <- propose2 <- fitfun <- list()
  for(i in names(x)) {
    theta[[i]] <- smooths[[i]] <- propose2[[i]] <- fitfun[[i]] <- list()
    nt <- NULL
    for(j in seq_along(x[[i]]$smooth)) {
      theta[[i]][[j]] <- x[[i]]$smooth[[j]]$state$parameters
      attr(theta[[i]][[j]], "fitted.values") <- x[[i]]$smooth[[j]]$state$fitted.values
      attr(theta[[i]][[j]], "hess") <- x[[i]]$smooth[[j]]$state$hessian
      x[[i]]$smooth[[j]]$state$fitted.values <- NULL
      x[[i]]$smooth[[j]]$XW <- t(x[[i]]$smooth[[j]]$X)
      x[[i]]$smooth[[j]]$XWX <- crossprod(x[[i]]$smooth[[j]]$X)
      x[[i]]$smooth[[j]]$fit.reduced <- as.numeric(rep(0, nrow(x[[i]]$smooth[[j]]$X)))
      nt <- c(nt, x[[i]]$smooth[[j]]$label)
      smooths[[i]][[j]] <- x[[i]]$smooth[[j]]
      if(!is.null(x[[i]]$smooth[[j]]$xt$propose))
        x[[i]]$smooth[[j]]$propose <- x[[i]]$smooth[[j]]$xt$propose
      propose2[[i]][[j]] <- if(is.null(x[[i]]$smooth[[j]]$propose)) propose else x[[i]]$smooth[[j]]$propose
      fitfun[[i]][[j]] <- function(x, p) {
        attr(p, "fitted.values")
      }
    }
    names(theta[[i]]) <- names(smooths[[i]]) <- names(propose2[[i]]) <- names(fitfun[[i]]) <- nt
  }
  rm(x)

  logLik <- function(eta) {
    family$loglik(response, family$map2par(eta))
  }

  n <- if(!is.null(dim(response))) nrow(response) else length(response)
  zworking <- as.numeric(rep(0, length = n))
  resids <- as.numeric(rep(0, length = n))

  samps <- gmcmc(fun = family, theta = theta, fitfun = fitfun, data = smooths,
    propose = propose2, logLik = logLik, n.iter = n.iter, burnin = burnin, thin = thin,
    response = response, simplify = FALSE, zworking = zworking, resids = resids,
    cores = cores, chains = chains, combine = FALSE, ...)

  samps
}


gmcmc <- function(fun, theta, priors = NULL, propose = NULL,
  fitfun = NULL, logLik = NULL, data = NULL, attr.copy = "hess",
  n.iter = 12000, burnin = 2000, thin = 10, verbose = TRUE, step = 20,
  simplify = TRUE, chains = NULL, cores = NULL, combine = TRUE, sleep = 1,
  compile = FALSE, ...)
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
    if(is.list(theta[[i]])) {
      if(is.null(names(theta[[i]])) & length(theta[[i]]))
        names(theta[[i]]) <- paste("term[", 1:length(theta[[i]]), "]", sep = "")
      if(length(theta[[i]])) {
        for(j in seq_along(theta[[i]])) {
          if(is.null(names(theta[[i]][[j]])) & length(theta[[i]][[j]]))
            names(theta[[i]][[j]]) <- paste("[", 1:length(theta[[i]][[j]]), "]", sep = "")
        }
      }
    }
  }
  class(theta) <- c("gmcmc.theta", "list")

  parse_input <- function(input = NULL, default, type) {
    ninput <- deparse(substitute(input), backtick = TRUE, width.cutoff = 500)
    if(is.null(input))
      input <- default
    if(!is.list(input))
      input <- rep(list(input), length.out = k)
    if(is.null(names(input)))
      names(input) <- names(theta)
    if(!all(names(input) %in% names(theta)))
      stop(paste("the '", ninput,"' list() names are different from theta!", sep = ""))
    for(i in names(theta)) {
      if(!is.list(input[[i]]) & is.list(theta[[i]]))
        input[[i]] <- rep(list(input[[i]]), length.out = length(theta[[i]]))
      if(is.list(theta[[i]])) {
        if(is.null(names(input[[i]])))
          names(input[[i]]) <- names(theta[[i]])
        if(!all(names(input[[i]]) %in% names(theta[[i]])))
          stop(paste(ninput, "list() names for parameter", i, "different from the theta list() names!"))
        if(type != "NA") {
          for(j in names(input[[i]])) {
            if(!inherits(input[[i]][[j]], type)) {
              stop(paste("the '", ninput, "' object for block [", i, "][", j,
                "] is not a ", type,"!", sep = ""), call. = FALSE)
            }
          }
        }
      } else {
        if(type != "NA") {
          if(!inherits(input[[i]], type)) {
            stop(paste("the '", ninput, "' object for block [", i,
              "] is not a ", type,"!", sep = ""), call. = FALSE)
          }
        }
      }
    }
    return(input)
  }

  if(is.null(fitfun)) {
    fitfun <- function(x, p) {
      if(is.null(x)) {
        return(sum(p))
      } else return(drop(x %*% p))
    }
  }

  propose <- parse_input(propose, gmcmc_mvnorm, "function")
  fitfun <- parse_input(fitfun, NULL, "function")
  if(!is.null(data))
    data <- parse_input(data, NULL, "NA")
  if(!is.null(priors))
    priors <- parse_input(priors, NULL, "function")

  if(burnin < 1) burnin <- 1
  iterthin <- as.integer(seq(burnin, n.iter, by = thin))

  rho <- new.env()
  theta.save <- vector(mode = "list", length = length(theta))
  names(theta.save) <- names(theta)
  eta <- theta.save
  for(i in ntheta) {
    eta[[i]] <- 0
    if(is.list(theta[[i]])) {
      for(j in names(theta[[i]]))
        eta[[i]] <- eta[[i]] + fitfun[[i]][[j]](data[[i]][[j]], theta[[i]][[j]])
    } else {
      eta[[i]] <- eta[[i]] + fitfun[[i]](data[[i]], theta[[i]])
    }
  }

  for(i in ntheta) {
    if(is.list(theta[[i]])) {
      theta.save[[i]] <- vector(mode = "list", length = length(theta[[i]]))
      names(theta.save[[i]]) <- names(theta[[i]])
      for(j in names(theta[[i]])) {
        p0 <- propose[[i]][[j]](fun, theta, id = c(i, j),
          prior = priors[[i]][[j]], data = data[[i]][[j]], eta = eta,
          iteration = 1, n.iter = n.iter, burnin = burnin, rho = rho, ...)
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
        if(length(p0$parameters) < 2 & FALSE) ## FIXME!!!
          names(p0$parameters) <- NULL
        theta.save[[i]][[j]] <- list(
          "samples" = matrix(NA, nrow = length(iterthin), ncol = length(p0$parameters)),
          "alpha" = rep(NA, length = length(iterthin)),
          "accepted" = rep(NA, length = length(iterthin))
        )
        colnames(theta.save[[i]][[j]]$samples) <- names(p0$parameters)
        if(!is.null(p0$extra)) {
          theta.save[[i]][[j]]$extra <- matrix(NA, nrow = length(iterthin), ncol = length(p0$extra))
          colnames(theta.save[[i]][[j]]$extra) <- names(p0$extra)
        }
        names(theta[[i]][[j]]) <- names(p0$parameters)
        if(!is.null(attr.copy)) {
          for(a in attr.copy)
            attr(theta[[i]][[j]], a) <- attr(p0$parameters, a)
        }
      }
    } else {
      p0 <- propose[[i]](fun, theta, id = i,
        prior = priors[[i]], data = data[[i]], eta = eta,
        iteration = 1, n.iter = n.iter, burnin = burnin, rho = rho, ...)
      if(!is.list(p0)) {
        stop(paste("the propose() function for block [", i,
          "] must return a named list()!", sep = ""))
      }
      if(is.null(p0$alpha)) {
        stop(paste("the propose() function for block [", i,
          "] must return the acceptance probability 'alpha'!", sep = ""))
      }
      if(is.null(p0$parameters)) {
        stop(paste("the propose() function for block [", i,
          "] must return a vector 'parameters'!", sep = ""))
      }
      if(length(p0$parameters) != length(theta[[i]])) {
        stop(paste("the propose() function for block [", i,
          "] must return a vector 'parameters' with the same length of initial parameters in theta!",
          sep = ""))
      }
      if(is.null(names(p0$parameters)))
        names(p0$parameters) <- paste("[", 1:length(p0$parameters), "]", sep = "")
      if(length(p0$parameters) < 2 & FALSE) ## FIXME!!!
        names(p0$parameters) <- NULL
      theta.save[[i]] <- list(
        "samples" = matrix(NA, nrow = length(iterthin), ncol = length(p0$parameters)),
        "alpha" = rep(NA, length = length(iterthin)),
        "accepted" = rep(NA, length = length(iterthin))
      )
      colnames(theta.save[[i]]$samples) <- names(p0$parameters)
      if(!is.null(p0$extra)) {
        theta.save[[i]]$extra <- matrix(NA, nrow = length(iterthin), ncol = length(p0$extra))
        colnames(theta.save[[i]]$extra) <- names(p0$extra)
      }
      names(theta[[i]]) <- names(p0$parameters)
      if(!is.null(attr.copy)) {
        for(a in attr.copy)
          attr(theta[[i]], a) <- attr(p0$parameters, a)
      }
    }
  }

  ll <- rep(NA, length = length(iterthin))

  nstep <- step
  step <- floor(n.iter / step)

  barfun <- function(ptm, n.iter, i, step, nstep, start = TRUE) {
    if(i == 10 & start) {
      cat("\r")
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

  sampler <- function(...) {
    cat2("Starting the sampler...")
    ptm <- proc.time()
    for(iter in 1:n.iter) {
      if(save <- iter %in% iterthin)
        js <- which(iterthin == iter)
      for(i in ntheta) {
        if(is.list(theta[[i]])) {
          for(j in names(theta[[i]])) {
            ## Get proposed states.
            state <- propose[[i]][[j]](fun, theta, id = c(i, j),
              prior = priors[[i]][[j]], data = data[[i]][[j]], eta = eta,
              iteration = iter, n.iter = n.iter, burnin = burnin, rho = rho, ...)

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
              if(!is.null(state$extra)) {
                theta.save[[i]][[j]]$extra[js, ] <- state$extra
              }
              ll[js] <- if(!is.null(logLik)) {
                logLik(eta)
              } else sum(do.call(fun, c(theta, list(...))[names(formals(fun))]), na.rm = TRUE)
            }
          }
        } else {
          ## Get proposed states.
          state <- propose[[i]](fun, theta, id = i,
            prior = priors[[i]], data = data[[i]], eta = eta,
            iteration = iter, n.iter = n.iter, burnin = burnin, rho = rho, ...)

          ## If accepted, set current state to proposed state.
          accepted <- if(is.na(state$alpha)) FALSE else log(runif(1)) <= state$alpha

          if(accepted) {
            eta[[i]] <- eta[[i]] - fitfun[[i]](data[[i]], theta[[i]])
            theta[[i]] <- state$parameters
            eta[[i]] <- eta[[i]] + fitfun[[i]](data[[i]], theta[[i]])
          }

          ## Save the samples.
          if(save) {
            theta.save[[i]]$samples[js, ] <- theta[[i]]
            theta.save[[i]]$alpha[js] <- min(c(exp(state$alpha), 1), na.rm = TRUE)
            theta.save[[i]]$accepted[js] <- accepted
            if(!is.null(state$extra)) {
              theta.save[[i]]$extra[js, ] <- state$extra
            }
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
      theta.save[[i]] <- if(is.list(theta[[i]])) {
        lapply(theta.save[[i]], function(x) { do.call("cbind", x) })
      } else do.call("cbind", theta.save[[i]])
      if(is.list(theta[[i]])) {
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
      } else {
        cn <- colnames(theta.save[[i]])
        sep <- c(".", "")[as.integer(grepl("[", cn, fixed = TRUE)) + 1L]
        for(j in seq_along(cn)) {
          cn[j] <- paste(if(length(theta) > 1) i else NULL,
            if(length(theta[[i]]) > 1) cn[j] else NULL,
            sep = if(length(theta) > 1) sep[j] else "")
        }
        colnames(theta.save[[i]]) <- cn
      }
      theta.save[[i]] <- if(length(theta.save[[i]]) < 2 & simplify) {
        if(is.list(theta[[i]])) theta.save[[i]][[1]] else theta.save[[i]]
      } else {
        if(is.list(theta[[i]])) do.call("cbind", theta.save[[i]]) else theta.save[[i]]
      }
    }
    theta.save <- if(length(theta.save) < 2) theta.save[[1]] else do.call("cbind", theta.save)
    theta.save <- cbind(theta.save, "logLik" = ll)

    theta.save <- mcmc(theta.save, start = burnin, end = n.iter, thin = thin)
    class(theta.save) <- c("gmcmc", class(theta.save))

    return(theta.save)
  }

  if(compile) {
    require("compiler")
    sampler <- cmpfun(sampler)
  }

  sampling <- function(chains, ...) {
    if(is.null(chains))  chains <- 1
    if(chains < 2) {
      return(sampler(...))
    } else {
      samps <- list()
      for(j in 1:chains)
        samps[[j]] <- sampler(...)
      samps <- as.mcmc.list(samps)
      return(samps)
    }
  }

  if(is.null(cores)) cores <- 1
  if(cores < 2) {
    rval <- sampling(chains, ...)
  } else {
    require("parallel")
    parallel_fun <- function(j) {
      if(j > 1 & !is.null(sleep)) Sys.sleep(sleep)
      sampling(chains, ...)
    }
    rval <- mclapply(1:cores, parallel_fun, mc.cores = cores)
    if(!inherits(rval, "mcmc.list") & inherits(rval[[1L]], "mcmc.list"))
      rval <- as.mcmc.list(do.call("c", rval))
    if(!inherits(rval, "mcmc.list"))
      rval <- as.mcmc.list(rval)
    if(combine)
      rval <- combine_chains(rval)
  }

  return(rval)
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

"[.gmcmc.theta" <- function(x, i) {
  if(length(i) < 2)
    return(x[[i]])
  else
    return(x[[i[1]]][[i[2]]])
}

"[<-.gmcmc.theta" <- function(x, i, value) {
  if(length(i) < 2)
    x[[i]] <- value
  else
    x[[i[1]]][[i[2]]] <- value
  return(x)
}

gmcmc_mvnorm <- function(fun, theta, id, prior, ...)
{
  require("mvtnorm")

  args <- list(...)
  iteration <- args$iteration

  if(is.null(attr(theta[id], "scale")))
    attr(theta[id], "scale") <- 1
  if(is.null(attr(theta[id], "P")))
    attr(theta[id], "sigma") <- diag(length(theta[id]))

  if(iteration < floor(0.15 * args$n.iter)) {
    scale <- attr(theta[id], "scale")
    sigma <- attr(theta[id], "sigma")

    k <- 1
    do <- TRUE
    while(do & k < 100) {
      ll0 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
      p0 <- if(is.null(prior)) {
        sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
      } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

      theta[id] <- drop(rmvnorm(n = 1, mean = theta[id], sigma = scale * sigma))

      ll1 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
      p1 <- if(is.null(prior)) {
        sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
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

    attr(theta[id], "scale") <- scale
    attr(theta[id], "sigma") <- sigma
  }

  scale <- attr(theta[id], "scale")
  sigma <- attr(theta[id], "sigma")

  ll0 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p0 <- if(is.null(prior)) {
    sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  theta[id] <- drop(rmvnorm(n = 1, mean = theta[id], sigma = scale * sigma))

  ll1 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p1 <- if(is.null(prior)) {
    sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  alpha <- drop((ll1 + p1) - (ll0 + p0))

  rval <- list("parameters" = theta[id], "alpha" = alpha)

  attr(rval$parameters, "scale") <- scale
  attr(rval$parameters, "sigma") <- sigma

  rval
}


gmcmc_unislice <- function(fun, theta, id, prior, j, ...,
  w = 1, m = 30, lower = -Inf, upper = +Inf)
{
  args <- list(...)

  x0 <- theta[id][j]
  gL <- gR <- theta

  fun2 <- function(theta) {
    ll <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
    lp <- if(is.null(prior)) {
      sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
    } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)
    return(ll + lp)
  }

  gx0 <- fun2(theta)

  ## Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)

  ## Find the initial interval to sample from.
  ## w <- w * abs(x0) ## FIXME???
  u <- runif(1, 0, w)
  gL[id][j] <- theta[id][j] - u
  gR[id][j] <- theta[id][j] + (w - u)  ## should guarantee that g[j] is in [L, R], even with roundoff

  ## Expand the interval until its ends are outside the slice, or until
  ## the limit on steps is reached.
  if(is.infinite(m)) {
    repeat {
      if(gL[id][j] <= lower) break
      if(fun2(gL) <= logy) break
      gL[id][j] <- gL[id][j] - w
    }
    repeat {
      if(gR[id][j] >= upper) break
      if(fun2(gR) <= logy) break
      gR[id][j] <- gR[id][j] + w
    }
  } else {
    if(m > 1) {
      J <- floor(runif(1, 0, m))
      K <- (m - 1) - J
      while(J > 0) {
        if(gL[id][j] <= lower) break
        if(fun2(gL) <= logy) break
        gL[id][j] <- gL[id][j] - w
        J <- J - 1
      }
      while(K > 0) {
        if(gR[id][j] >= upper) break
        if(fun2(gR) <= logy) break
        gR[id][j] <- gR[id][j] + w
        K <- K - 1
      }
    }
  }

  ## Shrink interval to lower and upper bounds.
  if(gL[id][j] < lower) {
    gL[id][j] <- lower
  }
  if(gR[id][j] > upper) {
    gR[id][j] <- upper
  }

  ## Sample from the interval, shrinking it on each rejection.
  repeat {
    theta[id][j] <- runif(1, gL[id][j], gR[id][j])

    gx1 <- fun2(theta)

    if(gx1 >= logy) break

    if(theta[id][j] > x0) {
      gR[id][j] <- theta[id][j]
    } else {
      gL[id][j] <- theta[id][j]
    }
  }

  ## Return the point sampled
  return(theta[id][j])
}


gmcmc_slice <- function(fun, theta, id, prior, ...)
{
  for(j in seq_along(theta[id]))
    theta[id][j] <- gmcmc_unislice(fun, theta, id, prior, j, ...)
  return(list("parameters" = theta[id], "alpha" = log(1)))
}


gmcmc_sm.iwlsC <- function(family, theta, id, prior,
  eta, response, data, zworking, resids, rho, ...)
{
  rval <- .Call("gmcmc_iwls", family, theta, id, eta, response, data, zworking, resids, rho)
  data$state$parameters <- as.numeric(rval$parameters)
  names(data$state$parameters) <- names(rval$parameters)
  rval$extra <- c("edf" = data$edf(data))
  rval
}

gmcmc_sm.iwls <- function(family, theta, id, prior, eta, response, data, ...)
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
  eta2 <- eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")

  ## Compute reduced residuals.
  e <- z - eta2
  xbin.fun(data$xbin.sind, weights, e, data$weights, data$rres, data$xbin.order)

  ## Compute mean and precision.
  XWX <- crossprod(data$X, data$X * data$weights)
  S <- 0
  P <- if(data$fixed) {
    if((k <- ncol(data$X)) < 2) {
      1 / XWX
    } else chol2inv(L1 <- chol(P01 <- XWX))
  } else {
    for(j in seq_along(data$S))
      S <- S + 1 / get.par(theta, "tau2")[j] * data$S[[j]]
    chol2inv(L1 <- chol(P01 <- XWX + S))
  }

  P[P == Inf] <- 0
  M <- P %*% crossprod(data$X, data$rres)

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

  ## Compute reduced residuals.
  e <- z - eta2
  xbin.fun(data$xbin.sind, weights, e, data$weights, data$rres, data$xbin.order)

  ## Compute mean and precision.
  XWX <- crossprod(data$X, data$X * data$weights)
  P2 <- if(data$fixed) {
    if(k < 2) {
      1 / (XWX)
    } else chol2inv(L2 <- chol(P02 <- XWX))
  } else {
    chol2inv(L2 <- chol(P02 <- XWX + S))
  }
  P2[P2 == Inf] <- 0
  M2 <- P2 %*% crossprod(data$X, data$rres)

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
  data$state$parameters <- as.numeric(theta)
  names(data$state$parameters) <- names(theta)

  ## Compute acceptance probablity.
  alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

  return(list("parameters" = theta, "alpha" = alpha, "extra" = c("edf" = data$edf(data))))
}


gmcmc_sm.mvn <- function(family, theta, id, prior, eta, response, data, ...)
{
  require("mvtnorm")

  theta <- theta[[id[1]]][[id[2]]]

  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$get.mu(data$X, theta)

  ll1 <- family$loglik(response, family$map2par(eta))
  p1 <- data$prior(theta)

  eta2 <- eta
  eta2[[id[1]]] <- eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")
  if(!is.null(data$state)) {
    if(!is.null(data$state$hessian))
      attr(theta, "hess") <- data$state$hessian
    if(is.null(attr(theta, "hess"))) {
      g <- get.par(data$state$parameters, "gamma")
      tau2 <- if(!data$fixed) get.par(data$state$parameters, "tau2") else NULL

      lp <- function(g) {
        eta[[id[1]]] <- eta[[id[1]]] + data$get.mu(data$X, g)
        family$loglik(response, family$map2par(eta)) + data$prior(c(g, tau2))
      }

      if(is.null(family$gradient[[id[1]]])) {
        gfun <- NULL
      } else {
        gfun <- list()
        gfun[[id[1]]] <- function(g, y, eta, x, ...) {
          gg <- family$gradient[[id[1]]](g, y, eta, x)
          if(!is.null(data$grad)) {
            gg <- gg + data$grad(score = NULL, c(g, tau2), full = FALSE)
          }
          drop(gg)
       }
      }

      if(is.null(family$hessian[[id[1]]])) {
        hfun <- NULL
      } else {
        hfun <- list()
        hfun[[id[1]]] <- function(g, y, eta, x, ...) {
          hg <- family$hessian[[id[1]]](g, y, eta, x)
          if(!is.null(data$hess)) {
            hg <- hg + data$hess(score = NULL, c(g, tau2), full = FALSE)
          }
          hg
        }
      }

      g.grad <- grad(fun = lp, theta = g, id = id[1], prior = NULL,
        args = list("gradient" = gfun, "x" = data, "y" = response, "eta" = eta))

      g.hess <- hess(fun = lp, theta = g, id = id[1], prior = NULL,
        args = list("gradient" = gfun, "hessian" = hfun, "x" = data, "y" = response, "eta" = eta))

      attr(theta, "hess") <- matrix_inv(g.hess)
    }
  }

  g <- drop(rmvnorm(n = 1, mean = drop(get.par(theta, "gamma")), sigma = attr(theta, "hess")))
  theta <- set.par(theta, g, "gamma")

  attr(theta, "fitted.values") <- data$get.mu(data$X, theta)

  eta[[id[1]]] <- eta[[id[1]]] + attr(theta, "fitted.values")
  ll2 <- family$loglik(response, family$map2par(eta))
  p2 <- data$prior(theta)

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

  ## Compute acceptance probablity.
  alpha <- drop((ll2 + p2) - (ll1 + p1))

  data$state$parameters <- as.numeric(theta)
  names(data$state$parameters) <- names(theta)

  rval <- list("parameters" = theta, "alpha" = alpha, "extra" = c("edf" = data$edf(data)))
  rval
}


gmcmc_sm.newton <- function(family, theta, id, prior, eta, response, data, ...)
{
  require("mvtnorm")

  theta <- theta[id]
  g <- get.par(theta, "g")
  tau2 <- if(!data$fixed) get.par(theta, "tau2") else NULL
  nu <- if(is.null(data$nu)) 0.1 else data$nu

  pibeta <- family$loglik(response, family$map2par(eta))
  p1 <- data$prior(theta)

  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$get.mu(data$X, theta)

  eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")

  lp <- function(g) {
    eta[[id[1]]] <- eta[[id[1]]] + data$get.mu(data$X, g)
    family$loglik(response, family$map2par(eta)) + data$prior(c(g, tau2))
  }

  if(is.null(family$gradient[[id[1]]])) {
    gfun <- NULL
  } else {
    gfun <- list()
    gfun[[id[1]]] <- function(g, y, eta, x, ...) {
      gg <- family$gradient[[id[1]]](g, y, eta, x)
      if(!is.null(data$grad)) {
        gg <- gg + data$grad(score = NULL, c(g, tau2), full = FALSE)
      }
      drop(gg)
    }
  }

  if(is.null(family$hessian[[id[1]]])) {
    hfun <- NULL
  } else {
    hfun <- list()
    hfun[[id[1]]] <- function(g, y, eta, x, ...) {
      hg <- family$hessian[[id[1]]](g, y, eta, x)
      if(!is.null(data$hess)) {
        hg <- hg + data$hess(score = NULL, c(g, tau2), full = FALSE)
      }
      hg
    }
  }

  g.grad <- grad(fun = lp, theta = g, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "x" = data, "y" = response, "eta" = eta))

  g.hess <- hess(fun = lp, theta = g, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "hessian" = hfun, "x" = data, "y" = response, "eta" = eta))

  Sigma <- matrix_inv(g.hess)
  mu <- drop(g + nu * Sigma %*% g.grad)

  g2 <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma))
  names(g2) <- names(g)

  p2 <- data$prior(c("g" = g2, tau2))
  qbetaprop <- dmvnorm(g2, mean = mu, sigma = Sigma, log = TRUE)

  g.grad2 <- grad(fun = lp, theta = g2, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "x" = data, "y" = response, "eta" = eta))

  g.hess2 <- hess(fun = lp, theta = g2, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "hessian" = hfun, "x" = data, "y" = response, "eta" = eta))

  Sigma2 <- matrix_inv(g.hess2)
  mu2 <- drop(g2 + nu * Sigma2 %*% g.grad2)

  theta <- set.par(theta, g2, "g")

  attr(theta, "fitted.values") <- data$get.mu(data$X, g2)
  eta[[id[1]]] <- eta[[id[1]]] + attr(theta, "fitted.values")

  pibetaprop <- family$loglik(response, family$map2par(eta))

  qbeta <- dmvnorm(matrix(g, nrow = 1), mean = mu2, sigma = Sigma2, log = TRUE)

  alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

  ## Sample variance parameter.
  if(!data$fixed & is.null(data$sp)) {
    if(!data$fixed & is.null(data$sp)) {
      tau2 <- NULL
      for(j in seq_along(data$S)) {
        a <- data$rank[j] / 2 + data$a
        b <- 0.5 * crossprod(g2, data$S[[j]]) %*% g2 + data$b
        tau2 <- c(tau2, 1 / rgamma(1, a, b))
      }
      theta <- set.par(theta, tau2, "tau2")
    }
  }

  data$state$parameters <- as.numeric(theta)
  names(data$state$parameters) <- names(theta)

  rval <- list("parameters" = theta, "alpha" = alpha, "extra" = c("edf" = data$edf(data)))
}


gmcmc_logPost <- function(g, x, family, response = NULL, eta = NULL, id, ll = NULL)
{
  if(is.null(ll)) {
    eta[[id[1]]] <- eta[[id[1]]] + x$get.mu(x$X, g)
    ll <- family$loglik(response, family$map2par(eta))
  }
  lp <- x$prior(g)
  return(ll + lp)
}


gmcmc_sm.slice <- function(family, theta, id, prior, eta, response, data, ...)
{
  theta <- theta[[id[1]]][[id[2]]]

  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$get.mu(data$X, theta)

  ## Remove fitted values.
  eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")

  ## Sample coefficients.
  data$state$parameters <- set.par(data$state$parameters, get.par(theta, "g"), "g")
  for(j in seq_along(get.par(theta, "gamma"))) {
    data$state$parameters <- uni.slice(data$state$parameters, data, family, response,
      eta, id[1], j, logPost = gmcmc_logPost)
  }

  ## New fitted values.
  fit <- data$get.mu(data$X, data$state$parameters)

  ## Sample variance parameter.
  if(!data$fixed & is.null(data$sp)) {
    if(length(i <- grep("tau2", names(theta)))) {
      eta[[id[1]]] <- eta[[id[1]]] + fit
      ll <- family$loglik(response, family$map2par(eta))
      for(j in i) {
        data$state$parameters <- uni.slice(data$state$parameters, data, family, NULL,
          NULL, id[1], j, logPost = gmcmc_logPost, lower = 0, ll = ll)
      }
    }
  }

  ## New theta.
  theta <- data$state$parameters
  attr(theta, "fitted.values") <- fit

  return(list("parameters" = theta, "alpha" = log(1), "extra" = c("edf" = data$edf(data))))
}


gmcmc_mvnorm2 <- function(fun, theta, id, prior, ...)
{
  require("mvtnorm")

  args <- list(...)
  iteration <- args$iteration

  ll0 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p0 <- if(is.null(prior)) {
    sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  adapt <- if(is.null(args$adapt)) args$burnin else args$adapt 
  if(iteration <= adapt) {
    eps <- attr(theta[id], "eps")
    if(is.null(eps))
      eps <- 1
    if(eps > 0.001) {
      objfun <- function(par) {
        theta[id] <- par
        ll <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
        lp <- if(is.null(prior)) {
          sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
        } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)
        return(-1 * (ll + lp))
      }
      start <- attr(theta[id], "mode")
      if(is.null(start))
        start <- theta[id]
      opt <- optim(start, fn = objfun, method = "L-BFGS-B", hessian = TRUE,
        control = list(maxit = 100))
      theta[id] <- opt$par
      attr(theta[id], "sigma") <- matrix_inv(diag(diag(opt$hessian)))
      attr(theta[id], "scale") <- 1
      attr(theta[id], "eps") <- mean(abs((start - opt$par) / start), na.rm = TRUE)
      attr(theta[id], "mode") <- opt$par
    }
  }

  theta.attr <- attributes(theta[id])

  scale <- theta.attr$scale
  sigma <- theta.attr$sigma
  sigma <- scale * sigma

  theta[id] <- drop(rmvnorm(n = 1, mean = theta[id], sigma = sigma))

  ll1 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p1 <- if(is.null(prior)) {
    sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  alpha <- drop((ll1 + p1) - (ll0 + p0))

  rval <- list("parameters" = theta[id], "alpha" = alpha)
  attributes(rval$parameters) <- theta.attr

  rval
}


eval_fun <- function(fun, theta, prior, id, args)
{
  ll <- if(is.list(theta)) {
    sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  } else sum(fun(theta), na.rm = TRUE)
  lp <- if(is.null(prior)) {
    if(is.list(theta)) {
      if(is.list(theta[[id[1]]])) {
        sum(dnorm(theta[[id[1]]][[id[2]]], sd = 1000, log = TRUE), na.rm = TRUE)
      } else {
        sum(dnorm(theta[[id[1]]], sd = 1000, log = TRUE), na.rm = TRUE)
      }
    } else {
      sum(dnorm(theta, sd = 1000, log = TRUE), na.rm = TRUE)
    }
  } else {
    if(is.list(theta)) {
      sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)
    } else {
      sum(prior(theta), na.rm = TRUE)
    }
  }
  ll + lp
}


grad <- function(fun, theta, prior, id, args, eps = .Machine$double.eps^0.25)
{
  if(!is.null(args$gradient[[id[1]]])) {
    gfun <- args$gradient[[id[1]]]
    args[[names(formals(gfun))[1]]] <- theta
    grad.theta <- do.call(gfun, args[names(formals(gfun))])
  } else {
    gfun <- function(theta, j) {
      theta2 <- theta
      if(is.list(theta)) {
        if(is.list(theta[[id[1]]])) {
          theta2[[id[1]]][[id[2]]][j] <- theta2[[id[1]]][[id[2]]][j] + eps
          theta[[id[1]]][[id[2]]][j] <- theta[[id[1]]][[id[2]]][j] - eps
        } else {
          theta2[[id[1]]][j] <- theta2[[id[1]]][j] + eps
          theta[[id[1]]][j] <- theta[[id[1]]][j] - eps
        }
      } else {
        theta2[j] <- theta2[j] + eps
        theta[j] <- theta[j] - eps
      }
      lp1 <- eval_fun(fun, theta2, prior, id, args)
      lp2 <- eval_fun(fun, theta, prior, id, args)
      drop((lp1 - lp2) / (2 * eps))
    }
    k <- if(is.list(theta)) {
      if(is.list(theta[[id[1]]])) {
        length(theta[[id[1]]][[id[2]]])
      } else length(theta[[id[1]]])
    } else length(theta)
    grad.theta <- rep(NA, k)
    for(j in 1:k) {
      grad.theta[j] <- gfun(theta, j)
    }
  }
  names(grad.theta) <- if(is.list(theta)) {
    if(is.list(theta[[id[1]]])) {
      names(theta[[id[1]]][[id[2]]])
    } else names(theta[[id[1]]])
  } else names(theta)
  grad.theta
}

hess <- function(fun, theta, prior, id, args, diag = FALSE)
{
  if(!is.null(args$hessian)) {
    hfun <- args$hessian[[id[1]]]
    args[[names(formals(hfun))[1]]] <- theta
    H <- do.call(hfun, args[names(formals(hfun))])
  } else {
    theta2 <- theta
    objfun <- function(par) {
      if(is.list(theta)) {
        if(is.list(theta[[id[1]]])) {
          theta2[[id[1]]][[id[2]]] <- par
        } else theta2[[id[1]]] <- par
        return(eval_fun(fun, theta2, prior, id, args))
      } else {
        return(eval_fun(fun, par, prior, id, args))
      }
    }
    gfun <- NULL
    if(!is.null(args$gradient)) {
      gfun2 <- args$gradient[[id[1]]]
      gfun <- function(par) {
        if(is.list(theta)) {
          if(is.list(theta[[id[1]]])) {
            theta2[[id[1]]][[id[2]]] <- par
          } else theta2[[id[1]]] <- par
          return(do.call(gfun2, c(theta2, args)[names(formals(gfun2))]))
        } else {
          args[[names(formals(gfun2))[1]]] <- par
          return(do.call(gfun2, args[names(formals(gfun2))]))
        }
      }
    }
    start <- if(is.list(theta)) {
      if(is.list(theta[[id[1]]])) {
        theta[[id[1]]][[id[2]]]
      } else theta[[id[1]]]
    } else theta
    H <- -1 * optimHess(start, fn = objfun, gr = gfun, control = list(fnscale = -1))
  }
  if(diag)
    H <- diag(diag(H))
  nh <- if(is.list(theta)) {
    if(is.list(theta[[id[1]]])) {
      names(theta[[id[1]]][[id[2]]])
    } else names(theta[[id[1]]])
  } else names(theta)
  rownames(H) <- colnames(H) <- nh
  H
}


gmcmc_newton <- function(fun, theta, id, prior, ...)
{
  require("mvtnorm")

  args <- list(...)

  adapt <- if(is.null(args$adapt)) {
    if(args$burnin < 1) 1 else args$burnin
  } else args$adapt
  eps0 <- if(is.null(args$eps0)) 1e-04 else args$eps0
  if(args$iteration <= adapt) {
    eps <- attr(theta[id], "eps")
    if(is.null(eps))
      eps <- 1
    if(is.na(eps))
      eps <- 1
    if(eps > eps0) {
      theta2 <- theta
      objfun <- function(par, ...) {
        theta2[id] <- par
        -1 * eval_fun(fun, theta2, prior, id, args)
      }
      gfun <- NULL
      if(!is.null(args$gradient)) {
        gfun2 <- args$gradient[[id[1]]]
        gfun <- function(par, ...) {
          theta2[id] <- par
          -1 * do.call(gfun2, c(theta2, args)[names(formals(gfun2))])
        }
      }

      start <- attr(theta[id], "mode")
      if(is.null(start))
        start <- theta[id]

      opt <- optim(start, fn = objfun, gr = gfun, method = "L-BFGS-B")
      
      rval <- list("parameters" = opt$par, "alpha" = log(10))
      attr(rval$parameters, "eps") <- mean(abs((start - opt$par) / start), na.rm = TRUE)
      attr(rval$parameters, "mode") <- opt$par
      mode(rval$parameters) <- "numeric"
    } else {
      rval <- list("parameters" = attr(theta[id], "mode"), "alpha" = log(10))
    }
  } else {
    p.old <- eval_fun(fun, theta, prior, id, args)
    theta2 <- theta

    grad.theta <- grad(fun, theta, prior, id, args)
    hess.theta <- hess(fun, theta, prior, id, args, diag = FALSE)

    Sigma <- matrix_inv(hess.theta)
    mu <- drop(theta[id] + Sigma %*% grad.theta)

    q.prop <- dmvnorm(matrix(theta[id], nrow = 1), mean = mu, sigma = Sigma, log = TRUE)

    theta2[id] <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma))

    grad.theta2 <- grad(fun, theta2, prior, id, args)
    hess.theta2 <- hess(fun, theta2, prior, id, args, diag = FALSE)

    Sigma2 <- matrix_inv(hess.theta2)
    mu2 <- drop(theta2[id] + Sigma2 %*% grad.theta2)

    p.prop <- eval_fun(fun, theta2, prior, id, args)

    q.old <- dmvnorm(matrix(theta[id], nrow = 1), mean = mu2, sigma = Sigma2, log = TRUE)

    alpha <- (p.prop - p.old) + (q.old - q.prop)

    rval <- list("parameters" = theta2[id], "alpha" = alpha)
  }

  rval
}


null.sampler <- function(x, n.samples = 500, criterion = c("AICc", "BIC", "AIC"), ...)
{
  family <- attr(x, "family")
  nx <- family$names
  np <- length(nx)
  par <- make_par(x, add.tau2 = TRUE)
  nh <- names(par$par)
  response <- attr(x, "response.vec")
  eta <- get.eta(x)
  if(n.samples < 1) n.samples <- 1

  if(n.samples > 1) {
    require("mvtnorm")
    no.special <- TRUE
    for(j in 1:np) {
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        if(inherits(x[[nx[j]]]$smooth[[sj]], "special"))
          no.special <- FALSE
      }
    }
    if(!is.null(family$hessian) & is.null(attr(x, "hessian")) & no.special) {
      hessian <- list(); i <- 1
      nh3 <- NULL
      for(j in 1:np) {
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          nh2 <- grep(paste("p", j, ".t", sj, ".", sep = ""), nh, fixed = TRUE, value = TRUE)
          nhg <- nh2[!grepl("tau2", nh2)]
          g <- get.state(x[[nx[j]]]$smooth[[sj]], "gamma")
          hessian[[i]] <- family$hessian[[nx[j]]](g, response, eta, x[[nx[j]]]$smooth[[sj]])
          nh3 <- c(nh3, nhg)
          i <- i + 1
        }
      }
      require("Matrix")
      hessian <- -1 * as.matrix(do.call("bdiag", hessian))
      rownames(hessian) <- colnames(hessian) <- nh3
    } else {
      hessian <- if(is.null(attr(x, "hessian"))) {
        opt0(x, hessian = TRUE, verbose = FALSE, ...)
      } else attr(x, "hessian")
    }
    hessian <- matrix_inv(-1 * hessian)
  }

  criterion <- match.arg(criterion)

  samps <- matrix(NA, n.samples, length(nh))
  colnames(samps) <- nh
  sn <- NULL

  edf <- 0
  edf_sm <- edf_sm_n <- NULL

  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      nh2 <- grep(paste("p", j, ".t", sj, ".", sep = ""), nh, fixed = TRUE, value = TRUE)
      nhg <- nh2[!grepl("tau2", nh2)]
      g <- get.state(x[[nx[j]]]$smooth[[sj]], "gamma")
      if(n.samples > 1) {
        sigma <- hessian[nhg, nhg, drop = FALSE]
        if(any(eigen(sigma, symmetric = TRUE)$values <= 0)) {
          require("Matrix")
          sigma2 <- try(nearPD(sigma)$mat, silent = TRUE)
          if(inherits(sigma2, "try-error")) {
            sigma2 <- diag(sigma)
            sigma2 <- if(length(sigma2) < 2) matrix(sigma2, 1, 1) else diag(sigma2)
          }
          sigma <- as.matrix(sigma2)
        }
        if(length(sigma) < 2) sigma <- matrix(sigma, 1, 1)
        g <- rmvnorm(n = n.samples, mean = g, sigma = sigma)
        colnames(g) <- nhg
        samps[, nhg] <- g
      } else samps[1L, nhg] <- g
      sn <- c(sn, paste(nx[j], x[[nx[j]]]$smooth[[sj]]$label, paste("g", 1:length(nhg), sep = ""), sep = "."))
      if(!x[[nx[j]]]$smooth[[sj]]$fixed) {
        nhtau2 <- nh2[grepl("tau2", nh2)]
        tedf <- x[[nx[j]]]$smooth[[sj]]$edf(x[[nx[j]]]$smooth[[sj]], type = 2)
        edf_sm <- cbind(edf_sm, rep(tedf, length = n.samples))
        edf_sm_n <- c(edf_sm_n, paste(nx[j], x[[nx[j]]]$smooth[[sj]]$label, "edf", sep = "."))
        samps[1L, nhtau2] <- get.state(x[[nx[j]]]$smooth[[sj]], "tau2") ##tedf
        edf <- edf + tedf
        sn <- c(sn, paste(nx[j], x[[nx[j]]]$smooth[[sj]]$label, paste("tau2", 1:length(nhtau2), sep = ""), sep = "."))
      } else {
        edf <- edf + ncol(x[[nx[j]]]$smooth[[sj]]$X)
      }
    }
  }

  colnames(edf_sm) <- edf_sm_n
  colnames(samps) <- sn
  samps <- cbind(samps, edf_sm)
  IC <- try(as.numeric(get.ic(family, response, family$map2par(eta), edf, length(eta[[1L]]), criterion)))
  if(inherits(IC, "try-error"))
    IC <- NA
  samps <- cbind(samps, IC, edf)
  colnames(samps)[(ncol(samps) - 1):ncol(samps)] <- c(criterion, "save.edf")

  as.mcmc(samps)
}

get.sigma <- function(x, family, response, eta, id)
{
  if(is.null(x$state$hessian)) {
    if(!is.null(family$hessian[[id]])) {
      xhess <- family$hessian[[id]](get.state(x, "g"), response, eta, x)
      xhess <- xhess + x$hess(score = NULL, x$state$parameters, full = FALSE)
    } else {
      eta[[id]] <- eta[[id]] - fitted(x$state)

      tau2 <- get.state(x, "tau2")

      lp <- function(g) {
        eta[[id]] <- eta[[id]] + x$get.mu(x$X, g)
        family$loglik(response, family$map2par(eta)) + x$prior(c(g, tau2))
      }

      if(is.null(family$gradient[[id]])) {
        if(!is.null(family$score[[id]]) & !is.null(x$grad)) {
          gfun <- function(g, ...) {
            eta[[id]] <- eta[[id]] + x$get.mu(x$X, g)
            x$grad(family$score[[id]](response, family$map2par(eta)), g, full = FALSE)
          }
        } else gfun <- NULL
      } else {
        gfun <- function(g, ...) {
          gg <- family$gradient[[id]](g, response, eta, x)
          if(!is.null(x$grad)) {
            gg <- gg + x$grad(score = NULL, c(g, tau2), full = FALSE)
          }
          drop(gg)
        }
      }

      xhess <- -1 * optimHess(get.state(x, "g"), fn = lp, gr = gfun, control = list("fnscale" = -1))
    }

    if(length(xhess) < 2) {
      xhess <- matrix(1 / xhess, 1, 1)
    } else {
      ## xhess <- chol2inv(chol(xhess))
      xhess <- diag(1 / diag(xhess))
    }

    return(xhess)
  } else return(x$state$hessian)
}


cat2 <- function(x, sleep = 0.03) {
  x <- strsplit(x, "")[[1]]
  add <- ""
  for(i in seq_along(x)) {
    if(i > 1) {
      Sys.sleep(sleep)
      cat("\r")
    }
    cat(paste(x[1:i], collapse = ""))
  }
}


## Survival propose() functions.
gmcmc_surv_sm.newton <- function(family, theta, id, prior, eta, response, data, ...)
{
  require("mvtnorm")

  args <- list(...)
  eta_Surv_timegrid0 <- eta_Surv_timegrid

  theta <- theta[id]
  g <- get.par(theta, "g")
  tau2 <- if(!data$fixed) get.par(theta, "tau2") else NULL
  nu <- if(is.null(data$nu)) 0.1 else data$nu

  p.old <- family$loglik(response, family$map2par(eta)) + data$prior(theta)

  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$get.mu(data$X, theta)
  if(is.null(attr(theta, "fitted_timegrid")))
    attr(theta, "fitted_timegrid") <- data$get.mu_timegrid(theta)

  eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")
  eta_Surv_timegrid <<- eta_Surv_timegrid - attr(theta, "fitted_timegrid")

  lp <- function(g) {
    eta[[id[1]]] <- eta[[id[1]]] + data$get.mu(data$X, g)
    family$loglik(response, family$map2par(eta)) + data$prior(c(g, tau2))
  }

  if(is.null(family$gradient[[id[1]]])) {
    gfun <- NULL
  } else {
    gfun <- list()
    gfun[[id[1]]] <- function(g, y, eta, x, ...) {
      gg <- family$gradient[[id[1]]](g, y, eta, x)
      if(!is.null(data$grad)) {
        gg <- gg + data$grad(score = NULL, c(g, tau2), full = FALSE)
      }
      drop(gg)
    }
  }

  if(is.null(family$hessian[[id[1]]])) {
    hfun <- NULL
  } else {
    hfun <- list()
    hfun[[id[1]]] <- function(g, y, eta, x, ...) {
      hg <- family$hessian[[id[1]]](g, y, eta, x)
      if(!is.null(data$hess)) {
        hg <- hg + data$hess(score = NULL, c(g, tau2), full = FALSE)
      }
      hg
    }
  }

  g.grad <- grad(fun = lp, theta = g, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "x" = data, "y" = response, "eta" = eta))

  g.hess <- hess(fun = lp, theta = g, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "hessian" = hfun, "x" = data, "y" = response, "eta" = eta))

  Sigma <- matrix_inv(g.hess)
  mu <- drop(g + nu * Sigma %*% g.grad)

  q.prop <- dmvnorm(matrix(g, nrow = 1), mean = mu, sigma = Sigma, log = TRUE)

  g2 <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma))
  names(g2) <- names(g)

  g.grad2 <- grad(fun = lp, theta = g2, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "x" = data, "y" = response, "eta" = eta))

  g.hess2 <- hess(fun = lp, theta = g2, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "hessian" = hfun, "x" = data, "y" = response, "eta" = eta))

  Sigma2 <- matrix_inv(g.hess2)
  mu2 <- drop(g2 + nu * Sigma2 %*% g.grad2)

  theta <- set.par(theta, g2, "g")

  attr(theta, "fitted.values") <- data$get.mu(data$X, g2)
  attr(theta, "fitted_timegrid") <- data$get.mu_timegrid(g2)
  eta[[id[1]]] <- eta[[id[1]]] + attr(theta, "fitted.values")
  eta_Surv_timegrid <<- eta_Surv_timegrid + attr(theta, "fitted_timegrid")

  p.prop <- family$loglik(response, family$map2par(eta)) + data$prior(c(g2, tau2))

  q.old <- dmvnorm(matrix(g2, nrow = 1), mean = mu2, sigma = Sigma2, log = TRUE)

  alpha <- (p.prop - p.old) + (q.old - q.prop)

  ## Sample variance parameter.
  if(!data$fixed & is.null(data$sp)) {
    if(!data$fixed & is.null(data$sp)) {
      tau2 <- NULL
      for(j in seq_along(data$S)) {
        a <- data$rank[j] / 2 + data$a
        b <- 0.5 * crossprod(g2, data$S[[j]]) %*% g2 + data$b
        tau2 <- c(tau2, 1 / rgamma(1, a, b))
      }
      theta <- set.par(theta, tau2, "tau2")
    }
  }

  accepted <- if(is.na(alpha)) FALSE else log(runif(1)) <= alpha
  if(args$iteration < 2 | !accepted)
    eta_Surv_timegrid <<- eta_Surv_timegrid0

  alpha <- if(accepted) 1 else -Inf

  data$state$parameters <- as.numeric(theta)
  names(data$state$parameters) <- names(theta)

  rval <- list("parameters" = theta, "alpha" = alpha, "extra" = c("edf" = data$edf(data)))
}

gmcmc_surv_sm.mvn <- function(family, theta, id, prior, eta, response, data, ...)
{
  require("mvtnorm")

  args <- list(...)
  eta_Surv_timegrid0 <- eta_Surv_timegrid

  theta <- theta[id]

  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$get.mu(data$X, theta)
  if(is.null(attr(theta, "fitted_timegrid")))
    attr(theta, "fitted_timegrid") <- data$get.mu_timegrid(theta)

  ll1 <- family$loglik(response, family$map2par(eta))
  p1 <- data$prior(theta)

  eta2 <- eta
  eta2[[id[1]]] <- eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")
  eta_Surv_timegrid <<- eta_Surv_timegrid - attr(theta, "fitted_timegrid")

  if(is.null(attr(theta, "hess"))) {
    g_opt <- get.par(data$state$parameters, "gamma")
    tau2_opt <- get.par(data$state$parameters, "tau2")
    hess <- family$hessian$lambda(g_opt, response, eta, data)
    if(!is.null(data$hess))
      hess <- hess + data$hess(score = NULL, c(g_opt, tau2_opt), full = FALSE)
    hess <- matrix_inv(hess)
    attr(theta, "hess") <- hess
  }

  g <- drop(rmvnorm(n = 1, mean = drop(get.par(theta, "gamma")), sigma = attr(theta, "hess")))
  theta <- set.par(theta, g, "gamma")

  attr(theta, "fitted.values") <- data$get.mu(data$X, theta)
  attr(theta, "fitted_timegrid") <- data$get.mu_timegrid(theta)

  eta[[id[1]]] <- eta[[id[1]]] + attr(theta, "fitted.values")
  eta_Surv_timegrid <<- eta_Surv_timegrid + attr(theta, "fitted_timegrid")
  ll2 <- family$loglik(response, family$map2par(eta))
  p2 <- data$prior(theta)

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

  ## Compute acceptance probablity.
  alpha <- drop((ll2 + p2) - (ll1 + p1))

  accepted <- if(is.na(alpha)) FALSE else log(runif(1)) <= alpha

  if(args$iteration < 2 | !accepted)
    eta_Surv_timegrid <<- eta_Surv_timegrid0

  alpha <- if(accepted) 1 else -Inf

  data$state$parameters <- as.numeric(theta)
  names(data$state$parameters) <- names(theta)

  rval <- list("parameters" = theta, "alpha" = alpha, "extra" = c("edf" = data$edf(data)))
  rval
}

