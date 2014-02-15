## Some useful links:
## http://adv-r.had.co.nz/C-interface.html
## http://stackoverflow.com/questions/7457635/calling-r-function-from-c
## http://gallery.rcpp.org/articles/r-function-from-c++/
propose_default <- function(x, family,
  response, eta, id, rho, ...)
{
  .Call("do_propose", x, family, response, eta, id, rho)
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
  family <- gaussian.BayesR()
  eta <- list("mu" = rep(0, n), "sigma" = rep(0, n))
  id <- "mu"

  a <- propose_default(x, family, response, eta, id, new.env())
}


## Setup for IWLS sampler, handling
## sampling functions.
transformIWLS <- function(x, ...)
{
  call <- x$call; x$call <- NULL

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
      df <- sum(diag(chol2inv(chol(XX + 1 / tau2 * x$S[[1]])) %*% XX))
      return((value - df)^2)
    }
    le <- optimize(objfun, c(lower, upper), value = 1)$minimum
    ri <- optimize(objfun, c(lower, upper), value = ncol(x$X))$minimum
    return(c(le, ri))
  }

  x$interval <- if(is.null(x$xt$interval)) tau2interval(x) else x$xt$interval
  x$grid <- if(is.null(x$xt$grid)) 40 else x$xt$grid 

  if(is.null(x$state)) {
    x$p.save <- c("g", "tau2")
    x$state <- list()
    x$state$g <- rep(0, ncol(x$X)) ##runif(ncol(x$X), 0.001, 0.002)
    x$state$tau2 <- if(is.null(x$sp)) {
      if(x$fixed) 1e-20 else 10
    } else x$sp
    if(x$fixed) x$state$tau2 <- 1e-20
    x$s.colnames <- if(is.null(x$s.colnames)) {
      c(paste("c", 1:length(x$state$g), sep = ""),
        if(!x$fixed) "tau2" else "tau2")
    } else x$s.colnames
    x$np <- length(x$s.colnames)
    XX <- crossprod(x$X)
    x$state$edf <- sum(diag(chol2inv(chol(XX + if(x$fixed) 0 else 1 / x$state$tau2 * x$S[[1]])) %*% XX))
  }
  if(is.null(x$xt$adaptive))
    x$xt$adaptive <- TRUE

  if(is.null(x$propose)) {
    if(TRUE) {
      x$propose <- propose_default
    } else {
    require("mvtnorm")
    x$propose <- function(x, family, response, eta, id, ...) {
      ## Compute weights.
      weights <- family$weights[[id]](response, eta)

      ## Score.
      score <- family$score[[id]](response, eta)

      ## Compute working observations.
      z <- eta[[id]] + 1 / weights * score

      ## Compute old log likelihood and old log coefficients prior.
      pibeta <- family$loglik(response, eta)
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

      ## Compute new log likelihood.
      pibetaprop <- family$loglik(response, eta)

      ## Compute new weights
      weights <- family$weights[[id]](response, eta)

      ## New score.
      score <- family$score[[id]](response, eta)

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
    }}
  }

  ## Function for computing starting values with backfitting.
  if(is.null(x$update)) {
    x$update <- function(x, family, response, eta, id, ...) {
      ## Compute weights.
      weights <- family$weights[[id]](response, eta)

      ## Score.
      score <- family$score[[id]](response, eta)

      ## Compute working observations.
      z <- eta[[id]] + 1 / weights * score

      ## Compute partial predictor.
      eta[[id]] <- eta[[id]] - x$state$fit

      ## Compute mean and precision.
      XW <- t(x$X * weights)
      XWX <- XW %*% x$X

      if(is.null(x$optimize) | x$fixed | !is.null(x$sp)) {
        P <- if(x$fixed) {
          chol2inv(chol(XWX))
        } else chol2inv(chol(XWX + 1 / x$state$tau2 * x$S[[1]]))
        x$state$g <- drop(P %*% (XW %*% (z - eta[[id]])))
      } else {
        args <- list(...)
        edf0 <- args$edf - x$state$edf
        eta2 <- eta
        e <- z - eta[[id]]

        objfun <- function(tau2, ...) {
          P <- chol2inv(chol(XWX + 1 / tau2 * x$S[[1]]))
          g <- drop(P %*% (XW %*% e))
          fit <- drop(x$X %*% g)
          edf <- sum(diag(P %*% XWX))
          if(!is.null(x$xt$center)) {
            if(x$xt$center) edf <- edf - 1
          }
          eta2[[id]] <- eta2[[id]] + fit
          IC <- get.ic(family, response, eta2, edf0 + edf, length(e), x$criterion)
          return(IC)
        }

        ## x$state$tau2 <- optimize(objfun, interval = x$interval, grid = x$grid)$minimum
        x$state$tau2 <- optimize2(objfun, interval = x$interval, grid = x$grid)$minimum
        P <- chol2inv(chol(XWX + 1 / x$state$tau2 * x$S[[1]]))
        x$state$g <- drop(P %*% (XW %*% e))
      }

      ## Compute fitted values.      
      x$state$fit <- drop(x$X %*% x$state$g)
      x$state$edf <- sum(diag(P %*% XWX))
      if(!is.null(x$xt$center)) {
        if(x$xt$center) x$state$edf <- x$state$edf - 1
      }

      return(x$state)
    }
  }

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
samplerIWLS <- function(x, n.iter = 12000, thin = 10, burnin = 2000,
  verbose = TRUE, step = 20, svalues = TRUE, eps = .Machine$double.eps^0.25, maxit = 400,
  tdir = NULL, method = c("backfitting", "MCMC", "backfitting2", "backfitting3"),
  n.samples = 200, criterion = c("AIC", "BIC", "AICc"), lower = 1e-09, upper = 1e+04,
  optim.control = list(pgtol = 1e-04, maxit = 5), digits = 3, ...)
{
  family <- attr(x, "family")
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  response <- attr(x, "response.vec")
  criterion <- match.arg(criterion)
  method <- match.arg(method, several.ok = TRUE)

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
  
  ## Add accptance rate and fitted values vectors.
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
          obj$smooth[[j]]$s.alpha <- rep(0, nrow = n.save)
          obj$smooth[[j]]$s.samples <- matrix(0, nrow = n.save, ncol = obj$smooth[[j]]$np)
          obj$smooth[[j]]$state$fit <- rep(0, nrow(obj$smooth[[j]]$X))
          obj$smooth[[j]]$fxsp <- if(!is.null(obj$smooth[[j]]$sp)) TRUE else FALSE
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
    eta[[j]] <- rep(0.1, length(response))

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
      formatC(round(x, digits), width = width)
    }

    ## Backfitting main function.
    backfit <- function(x, eta, verbose) {
      edf <- get_edf(x)

      eps0 <- iter <- 1
      while(eps0 > eps & iter < maxit) {
        eta0 <- eta
        ## Cycle through all parameters
        for(j in 1:np) {
          ## And all terms.
          for(sj in seq_along(x[[nx[j]]]$smooth)) {
            ## Get updated parameters.
            p.state <- x[[nx[j]]]$smooth[[sj]]$update(x[[nx[j]]]$smooth[[sj]],
              family, response, eta, nx[j], edf = edf)

            ## Compute equivalent degrees of freedom.
            edf <- edf - x[[nx[j]]]$smooth[[sj]]$state$edf + p.state$edf

            ## Update predictor and smooth fit.
            eta[[nx[j]]] <- eta[[nx[j]]] - x[[nx[j]]]$smooth[[sj]]$state$fit + p.state$fit
            x[[nx[j]]]$smooth[[sj]]$state <- p.state
          }
        }

        eps0 <- mean((do.call("cbind", eta) - do.call("cbind", eta0))^2)

        if(any(method == "backfitting") & verbose) {
          IC <- get.ic(family, response, eta, edf, nobs, criterion)

          cat("\r")
          vtxt <- paste(criterion, " ", fmt(IC, width = -8, digits = digits),
            " loglik ", fmt(family$loglik(response, eta), width = -8, digits = digits),
            " edf ", fmt(edf, width = -6, digits = digits + 2),
            " eps ", fmt(eps0, width = -5, digits = digits * 2),
            " iteration ", formatC(iter, width = -1 * nchar(maxit)), sep = "")
          cat(vtxt)

          if(.Platform$OS.type != "unix") flush.console()
        }

        iter <- iter + 1
      }

      IC <- get.ic(family, response, eta, edf, nobs, criterion)

      if(any(method %in% c("backfitting", "backfitting2")) & verbose) {
        cat("\r")
        vtxt <- paste(criterion, " ", fmt(IC, width = -8, digits = digits),
          " loglik ", fmt(family$loglik(response, eta), width = -8, digits = digits),
          " edf ", fmt(edf, width = -6, digits = digits + 2),
          " eps ", fmt(eps0, width = -5, digits = digits * 2),
          " iteration ", formatC(iter, width = -1 * nchar(maxit)), sep = "")
        cat(vtxt)
        if(.Platform$OS.type != "unix") flush.console()
        cat("\n")
      }

      if(iter == maxit)
        warning("the backfitting algorithm did not converge, please check argument eps and maxit!")

      return(list("x" = x, "eta" = eta, "ic" = IC))
    }

    if(length(method) < 2) {
      verbose2 <- if(method == "MCMC")  FALSE else verbose
    }

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
          p.state <- x[[nx[j]]]$smooth[[sj]]$propose(x[[nx[j]]]$smooth[[sj]], family, response, eta, nx[j], rho = rho)

          ## If accepted, set current state to proposed state.
          accepted <- if(is.na(p.state$alpha)) FALSE else log(runif(1)) <= p.state$alpha

          if(accepted) {
            eta[[nx[j]]] <- eta[[nx[j]]] - x[[nx[j]]]$smooth[[sj]]$state$fit + p.state$fit
            x[[nx[j]]]$smooth[[sj]]$state <- p.state 
          }
          x[[nx[j]]]$smooth[[sj]]$state$accepted <- accepted

          ## Save the samples and acceptance.
          if(save) {
            x[[nx[j]]]$smooth[[sj]]$s.alpha[js] <- exp(p.state$alpha)
            x[[nx[j]]]$smooth[[sj]]$s.samples[js, ] <- unlist(x[[nx[j]]]$smooth[[sj]]$state[x[[nx[j]]]$smooth[[sj]]$p.save])
          }
        }
      }

      if(save) deviance[js] <- -2 * family$loglik(response, eta)

      if(verbose) barfun(ptm, n.iter, i, step, nstep)
    }
    if(verbose) cat("\n")
  }

  save.edf <- save.loglik <- NULL
  if(!any(grepl("MCMC", method))) {
    save.edf <- get_edf(x)
    save.loglik <- family$loglik(response, eta)
    eta2 <- eta
    if(verbose) cat("generating samples\n")
    ptm <- proc.time()
    for(js in seq_along(iterthin)) {
      for(j in 1:np) {
        ## And all terms.
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          ## Get proposed states.
          p.state <- x[[nx[j]]]$smooth[[sj]]$propose(x[[nx[j]]]$smooth[[sj]], family, response, eta, nx[j], rho = rho)

          ## Save the samples.
          x[[nx[j]]]$smooth[[sj]]$s.alpha[js] <- exp(p.state$alpha)
          x[[nx[j]]]$smooth[[sj]]$s.samples[js, ] <- unlist(p.state[x[[nx[j]]]$smooth[[sj]]$p.save])

          eta2[[nx[j]]] <- eta2[[nx[j]]] - x[[nx[j]]]$smooth[[sj]]$state$fit + p.state$fit

          deviance[js] <- -2 * family$loglik(response, eta2)
        }
      }
      if(verbose) barfun(ptm, length(iterthin), js, 1, 20, start = FALSE)
    }
    if(verbose) cat("\n")
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
    if(inherits(samples[[1]], "mcmc.list"))
      samples <- do.call("c", samples)
    else samples <- as.mcmc.list(list(samples))
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
          k <- ncol(obj$smooth[[i]]$X)
          pn <- grep(paste(id, "h1", obj$smooth[[i]]$label, sep = ":"), snames,
            value = TRUE, fixed = TRUE) ## FIXME: hlevels!
          pn <- pn[!grepl("tau2", pn) & !grepl("alpha", pn)]
          psamples <- matrix(samples[[j]][, snames %in% pn], ncol = k)

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
          }

          ## Acceptance probalities.
          asamples <- NULL
          alpha <- paste(id, "h1", paste(obj$smooth[[i]]$label, "alpha", sep = "."), sep = ":")
          if(length(alpha <- grep(alpha, snames, fixed = TRUE))) {
            asamples <- as.numeric(samples[[j]][, alpha])
          }

          ## Prediction matrix.
          get.X <- function(x) {
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
            qu <- t(apply(psamples, 2, quantile, probs = c(0.025, 0.5, 0.975)))
            sd <- drop(apply(psamples, 2, sd))
            me <- drop(apply(psamples, 2, mean))
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

