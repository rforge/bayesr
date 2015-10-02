#####################
## Survival models ##
#####################
###############################
## Cox model fitting engine. ##
###############################
## (1) The family object.
cox.bamlss <- function(...)
{
  require("survival")

  rval <- list(
    "family" = "cox",
    "names" = c("lambda", "gamma"),
    "links" = c(lambda = "log", gamma = "identity"),
    "transform" = function(x, ...) {
      surv.transform(x = x$x, y = x$y, data = model.frame(x), family = x$family, is.cox = TRUE, ...)
    },
    "optimizer" = cox.mode,
    "sampler" = cox.mcmc,
    "loglik" = function(y, eta, ...) {
      n <- attr(y, "subdivisions")
      eeta <- exp(eta_Surv_timegrid)
      int <- attr(y, "width") * (0.5 * (eeta[, 1] + eeta[, n]) + apply(eeta[, 2:(n - 1)], 1, sum))
      ll <- (eta$lambda + eta$gamma) * y[, "status"] - exp(eta$gamma) * int
      sum(ll)
    },
    "predict" = bamlss.surv.prob
  )

  class(rval) <- "family.bamlss"
  rval
}


## Posterior mode estimation.
cox.mode <- function(x, y, weights, offset,
  nu = 0.1, update.nu = TRUE, eps = .Machine$double.eps^0.25, maxit = 400,
  verbose = TRUE, digits = 4, ...)
{
  if(nu < 0 | nu > 1)
    stop("nu must be 0 < nu < 1!")

  ## Names of parameters/predictors.
  nx <- names(x)
  
  ## Compute additive predictors.
  eta <- get.eta(x)

  ## For the time dependent part, compute
  ## predictor based on the time grid.
  eta_timegrid <- 0
  for(sj in seq_along(x$lambda$smooth.construct)) {
    g <- get.state(x$lambda$smooth.construct[[sj]], "b")
    eta_timegrid <- eta_timegrid + x$lambda$smooth.construct[[sj]]$fit.fun_timegrid(g)
  }

  ## Extract y.
  y <- y[[1]]

  ## Number of observations.
  nobs <- nrow(y)

  ## Number of subdivions used for the time grid.
  sub <- attr(y, "subdivisions")

  ## The interval width from subdivisons.
  width <- attr(y, "width")

  ## Start the backfitting algorithm.
  logPost0 <- NA
  eps0 <- eps + 1; iter <- 1
  while(eps0 > eps & iter < maxit) {
    eta0 <- eta

    ########################################
    ## Cycle through time-dependent part. ##
    ########################################
    for(sj in seq_along(x$lambda$smooth.construct)) {
      ## The time-dependent design matrix for the grid.
      X <- x$lambda$smooth.construct[[sj]]$fit.fun_timegrid(NULL)

      ## Timegrid lambda.
      eeta <- exp(eta_timegrid)

      ## Compute gradient and hessian integrals.
      int <- survint(X, eeta, width, exp(eta$gamma))
      xgrad <- drop(t(y[, "status"]) %*% x$lambda$smooth.construct[[sj]]$XT - int$grad)
      xgrad <- xgrad + x$lambda$smooth.construct[[sj]]$grad(score = NULL, x$lambda$smooth.construct[[sj]]$state$parameters, full = FALSE)
      xhess <- int$hess + x$lambda$smooth.construct[[sj]]$hess(score = NULL, x$lambda$smooth.construct[[sj]]$state$parameters, full = FALSE)

      ## Compute the inverse of the hessian.
      Sigma <- matrix_inv(xhess)

      ## Update regression coefficients.
      g <- get.state(x$lambda$smooth.construct[[sj]], "b")
      g2 <- drop(g + nu * Sigma %*% xgrad)
      names(g2) <- names(g)
      x$lambda$smooth.construct[[sj]]$state$parameters <- set.par(x$lambda$smooth.construct[[sj]]$state$parameters, g2, "b")

      ## Update additive predictors.
      fit_timegrid <- x$lambda$smooth.construct[[sj]]$fit.fun_timegrid(g2)
      eta_timegrid <- eta_timegrid - x$lambda$smooth.construct[[sj]]$state$fitted_timegrid + fit_timegrid
      x$lambda$smooth.construct[[sj]]$state$fitted_timegrid <- fit_timegrid

      fit <- x$lambda$smooth.construct[[sj]]$fit.fun(x$lambda$smooth.construct[[sj]]$X, g2)
      eta$lambda <- eta$lambda - fitted(x$lambda$smooth.construct[[sj]]$state) + fit
      x$lambda$smooth.construct[[sj]]$state$fitted.values <- fit
      x$lambda$smooth.construct[[sj]]$state$edf <- sum.diag(int$hess %*% Sigma)
    }

    ###########################################
    ## Actual integral of survivor function. ##
    ###########################################
    eeta <- exp(eta_timegrid)
    int <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))

    ##########################################
    ## Cycle through time-independent part. ##
    ##########################################
    for(sj in seq_along(x$gamma$smooth.construct)) {
      ## Compute weights.
      weights <- exp(eta$gamma) * int

      ## Compute score.
      score <- y[, "status"] - exp(eta$gamma) * int

      ## Compute working observations.
      z <- eta$gamma + 1 / weights * score

      ## Compute partial predictor.
      eta$gamma <- eta$gamma - fitted(x$gamma$smooth.construct[[sj]]$state)

      ## Compute reduced residuals.
      e <- z - eta$gamma
      xbin.fun(x$gamma$smooth.construct[[sj]]$binning$sorted.index, weights, e,
        x$gamma$smooth.construct[[sj]]$weights, x$gamma$smooth.construct[[sj]]$rres,
        x$gamma$smooth.construct[[sj]]$binning$order)

      ## Compute mean and precision.
      XWX <- do.XWX(x$gamma$smooth.construct[[sj]]$X,
        1 / x$gamma$smooth.construct[[sj]]$weights,
        x$gamma$smooth.construct[[sj]]$imat)
      if(x$gamma$smooth.construct[[sj]]$fixed) {
        P <- matrix_inv(XWX)
      } else {
        S <- 0
        tau2 <- get.state(x$gamma$smooth.construct[[sj]], "tau2")
        for(j in seq_along(x$gamma$smooth.construct[[sj]]$S))
          S <- S + 1 / tau2[j] * x$gamma$smooth.construct[[sj]]$S[[j]]
        P <- matrix_inv(XWX + S)
      }
      g <- drop(P %*% crossprod(x$gamma$smooth.construct[[sj]]$X, x$gamma$smooth.construct[[sj]]$rres))
      x$gamma$smooth.construct[[sj]]$state$parameters <- set.par(x$gamma$smooth.construct[[sj]]$state$parameters, g, "b")
      x$gamma$smooth.construct[[sj]]$state$edf <- sum.diag(XWX %*% P)

      ## Compute fitted values.
      if(any(is.na(g)) | any(g %in% c(-Inf, Inf))) {
        x$gamma$smooth.construct[[sj]]$state$parameters <- set.par(x$gamma$smooth.construct[[sj]]$state$parameters,
          rep(0, length(x$gamma$smooth.construct[[sj]]$state$g)), "b")
      }
      x$gamma$smooth.construct[[sj]]$state$fitted.values <- x$gamma$smooth.construct[[sj]]$fit.fun(x$gamma$smooth.construct[[sj]]$X,
        get.state(x$gamma$smooth.construct[[sj]], "b"))

      ## Update additive predictor.
      eta$gamma <- eta$gamma + fitted(x$gamma$smooth.construct[[sj]]$state)
    }

    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

    logLik <- sum((eta$lambda + eta$gamma) * y[, "status"] - exp(eta$gamma) * int, na.rm = TRUE)
    logPost <- as.numeric(logLik + get.log.prior(x))

    if(iter > 1 & update.nu) {
      if(logPost < logPost0) {
        nu <- nu * 0.9
        if(verbose)
          cat("\nupdated nu to:", nu, "\n")
      }
    }

    if(verbose) {
      cat("\r")
      vtxt <- paste(
        "logPost ", fmt(logPost, width = 8, digits = digits),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)), sep = ""
      )
      cat(vtxt)
      if(.Platform$OS.type != "unix") flush.console()
    }

    logPost0 <- logPost

    iter <- iter + 1
  }

  if(iter == maxit)
    warning("the backfitting algorithm did not converge, please check argument eps and maxit!")

  if(verbose) cat("\n")

  logLik <- sum((eta$lambda + eta$gamma) * y[, "status"] - exp(eta$gamma) * int, na.rm = TRUE)
  logPost <- as.numeric(logLik + get.log.prior(x))

  return(list("fitted.values" = eta, "parameters" = get.all.par(x),
    "edf" = get.edf(x, type = 2), "logLik" = logLik, "logPost" = logPost))

  return(x)
}


## The MCMC sampling engine.
## Posterior mode estimation.
cox.mcmc <- function(x, y, family, start, weights, offset,
  n.iter = 1200, burnin = 200, thin = 1,
  verbose = TRUE, digits = 4, step = 20, ...)
{
  require("mvtnorm")

  nu <- 1

  if(!is.null(start))
    x <- set.starting.values(x, start)

  ## Names of parameters/predictors.
  nx <- names(x)
  
  ## Compute additive predictors.
  eta <- get.eta(x)

  ## For the time dependent part, compute
  ## predictor based on the time grid.
  eta_timegrid <- 0
  for(sj in seq_along(x$lambda$smooth.construct)) {
    g <- get.state(x$lambda$smooth.construct[[sj]], "b")
    fit_timegrid <- x$lambda$smooth.construct[[sj]]$fit.fun_timegrid(g)
    x$lambda$smooth.construct[[sj]]$state$fitted_timegrid <- fit_timegrid
    eta_timegrid <- eta_timegrid + x$lambda$smooth.construct[[sj]]$fit.fun_timegrid(g)
  }

  ## Extract y.
  y <- y[[1]]

  ## Number of observations.
  nobs <- nrow(y)

  ## Number of subdivions used for the time grid.
  sub <- attr(y, "subdivisions")

  ## The interval width from subdivisons.
  width <- attr(y, "width")

  ## Porcess iterations.
  if(burnin < 1) burnin <- 1
  if(burnin > n.iter) burnin <- floor(n.iter * 0.1)
  if(thin < 1) thin <- 1
  iterthin <- as.integer(seq(burnin, n.iter, by = thin))

  ## Samples.
  samps <- list()
  for(i in nx) {
    samps[[i]] <- list()
    for(j in names(x[[i]]$smooth.construct)) {
      samps[[i]][[j]] <- list(
        "samples" = matrix(NA, nrow = length(iterthin), ncol = length(x[[i]]$smooth.construct[[j]]$state$parameters)),
        "edf" = rep(NA, length = length(iterthin)),
        "alpha" = rep(NA, length = length(iterthin)),
        "accepted" = rep(NA, length = length(iterthin))
      )
      colnames(samps[[i]][[j]]$samples) <- names(x[[i]]$smooth.construct[[j]]$state$parameters)
    }
  }
  logLik.samps <- logPost.samps <- rep(NA, length = length(iterthin))

  foo <- function(x) {
    x <- exp(x)
    if(x < 0)
      x <- 0
    if(x > 1)
      x <- 1
    x
  }

  ## Start sampling.
  cat2("Starting the sampler...")

  nstep <- step
  step <- floor(n.iter / step)

  ptm <- proc.time()
  for(iter in 1:n.iter) {
    if(save <- iter %in% iterthin)
      js <- which(iterthin == iter)
    ########################################
    ## Cycle through time-dependent part. ##
    ########################################
    for(sj in names(x$lambda$smooth.construct)) {
      p.state <- propose_surv_td(x$lambda$smooth.construct[[sj]], y, eta, eta_timegrid, width, sub, nu)

      ## If accepted, set current state to proposed state.
      accepted <- if(is.na(p.state$alpha)) FALSE else log(runif(1)) <= p.state$alpha

      if(accepted) {
        eta_timegrid <- eta_timegrid - x$lambda$smooth.construct[[sj]]$state$fitted_timegrid + p.state$fitted_timegrid
        eta$lambda <- eta$lambda - fitted(x$lambda$smooth.construct[[sj]]$state) + fitted(p.state)
        x$lambda$smooth.construct[[sj]]$state <- p.state 
      }

      ## Save the samples and acceptance.
      if(save) {
        samps$lambda[[sj]]$samples[js, ] <- x$lambda$smooth.construct[[sj]]$state$parameters
        samps$lambda[[sj]]$edf[js] <- x$lambda$smooth.construct[[sj]]$state$edf
        samps$lambda[[sj]]$alpha[js] <- foo(p.state$alpha)
        samps$lambda[[sj]]$accepted[js] <- accepted
      }
    }

    ###########################################
    ## Actual integral of survivor function. ##
    ###########################################
    eeta <- exp(eta_timegrid)
    int <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))

    ##########################################
    ## Cycle through time-independent part. ##
    ##########################################
    for(sj in names(x$gamma$smooth.construct)) {
      p.state <- propose_surv_tc(x$gamma$smooth.construct[[sj]], y, eta, int)

      ## If accepted, set current state to proposed state.
      accepted <- if(is.na(p.state$alpha)) FALSE else log(runif(1)) <= p.state$alpha

      if(accepted) {
        eta$gamma <- eta$gamma - fitted(x$gamma$smooth.construct[[sj]]$state) + fitted(p.state)
        x$gamma$smooth.construct[[sj]]$state <- p.state 
      }

      ## Save the samples and acceptance.
      if(save) {
        samps$gamma[[sj]]$samples[js, ] <- x$gamma$smooth.construct[[sj]]$state$parameters
        samps$gamma[[sj]]$edf[js] <- x$gamma$smooth.construct[[sj]]$state$edf
        samps$gamma[[sj]]$alpha[js] <- foo(p.state$alpha)
        samps$gamma[[sj]]$accepted[js] <- accepted
      }
    }

    if(save) {
      logLik.samps[js] <- sum((eta$lambda + eta$gamma) * y[, "status"] - exp(eta$gamma) * int, na.rm = TRUE)
      logPost.samps[js] <- as.numeric(logLik.samps[js] + get.log.prior(x))
    }

    if(verbose) barfun(ptm, n.iter, iter, step, nstep)
  }

  if(verbose) cat("\n")

  for(i in names(samps)){
    for(j in names(samps[[i]])) {
      samps[[i]][[j]] <- do.call("cbind", samps[[i]][[j]])
      cn <- if(j == "model.matrix") {
        paste(i, "p", j, colnames(samps[[i]][[j]]), sep = ".")
      } else {
        paste(i, "s", j, colnames(samps[[i]][[j]]), sep = ".")
      }
      colnames(samps[[i]][[j]]) <- cn
    }
    samps[[i]] <- do.call("cbind", samps[[i]])
  }
  samps$logLik <- logLik.samps
  samps$logPost <- logPost.samps
  samps <- do.call("cbind", samps)

  ## Compute DIC.
  dev <- -2 * logLik.samps
  mpar <- apply(samps, 2, mean, na.rm = TRUE)
  names(mpar) <- colnames(samps)
  ll <- sum(family$p2d(mpar, log = TRUE), na.rm = TRUE)
  mdev <- -2 * ll
  pd <- mean(dev) - mdev
  DIC <- mdev + 2 * pd

  samps <- cbind(samps, "DIC" = DIC, "pd" = pd)

  return(as.mcmc(samps))
}


## MCMC propose functions.
## Time-varying propose function.
propose_surv_td <- function(x, y, eta, eta_timegrid, width, sub, nu)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)

  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)

  ## Old logLik and prior.
  int <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibeta <- sum((eta$lambda + eta$gamma) * y[, "status"] - exp(eta$gamma) * int, na.rm = TRUE)
  p1 <- x$prior(x$state$parameters)

  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta, width, exp(eta$gamma), index = x$imat)
  xgrad <- drop(t(y[, "status"]) %*% x$XT - int$grad)
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)

  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, index = x$imat)

  ## Save old coefficients.
  g0 <- get.state(x, "b")

  ## Get new position.
  mu <- drop(g0 + nu * Sigma %*% xgrad)

  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma))
  names(g) <- names(g0)
  x$state$parameters <- set.par(x$state$parameters, g, "b")

  ## Compute log priors.
  p2 <- x$prior(x$state$parameters)
  qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)

  ## Update additive predictors.
  fit_timegrid <- x$fit.fun_timegrid(g)
  eta_timegrid <- eta_timegrid - x$state$fitted_timegrid + fit_timegrid
  x$state$fitted_timegrid <- fit_timegrid

  fit <- x$fit.fun(x$X, g)
  eta$lambda <- eta$lambda - fitted(x$state) + fit
  x$state$fitted.values <- fit

  ## New logLik.
  eeta <- exp(eta_timegrid)
  int <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibetaprop <- sum((eta$lambda + eta$gamma) * y[, "status"] - exp(eta$gamma) * int, na.rm = TRUE)

  ## Prior prob.
  int <- survint(X, eeta, width, exp(eta$gamma), index = x$imat)
  xgrad <- drop(t(y[, "status"]) %*% x$XT - int$grad)
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)

  Sigma2 <- matrix_inv(xhess, index = x$imat)
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  qbeta <- dmvnorm(g0, mean = mu2, sigma = Sigma2, log = TRUE)

  ## Save edf.
  x$state$edf <- sum.diag(int$hess %*% Sigma2)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp)) {
    if(!x$fixed & is.null(x$sp)) {
      tau2 <- NULL
      for(j in seq_along(x$S)) {
        a <- x$rank[j] / 2 + x$a
        b <- 0.5 * crossprod(g, x$S[[j]]) %*% g + x$b
        tau2 <- c(tau2, 1 / rgamma(1, a, b))
      }
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    }
  }

  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

  return(x$state)
}


## Time-constant propose function.
propose_surv_tc <- function(x, y, eta, int)
{
  ## Compute weights.
  weights <- exp(eta$gamma) * int

  ## Compute score.
  score <- y[, "status"] - exp(eta$gamma) * int

  ## Compute working observations.
  z <- eta$gamma + 1 / weights * score

  ## Compute old log likelihood and old log coefficients prior.
  pibeta <- sum((eta$lambda + eta$gamma) * y[, "status"] - exp(eta$gamma) * int, na.rm = TRUE)
  p1 <- x$prior(x$state$parameters)

  ## Compute partial predictor.
  eta2 <- eta$gamma <- eta$gamma - fitted(x$state)

  ## Compute reduced residuals.
  e <- z - eta2
  xbin.fun(x$binning$sorted.index, weights, e, x$weights, x$rres, x$binning$order)

  ## Compute mean and precision.
  XWX <- do.XWX(x$X, 1 / x$weights, x$imat)
  S <- 0
  P <- if(x$fixed) {
    if((k <- ncol(x$X)) < 2) {
      1 / XWX
    } else chol2inv(L1 <- chol(P01 <- XWX))
  } else {
    tau2 <- get.par(x$state$parameters, "tau2")
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    chol2inv(L1 <- chol(P01 <- XWX + S))
  }

  P[P == Inf] <- 0
  M <- P %*% crossprod(x$X, x$rres)

  ## Degrees of freedom.
  x$state$edf <- sum.diag(XWX %*% P)

  ## Save old coefficients
  g0 <- drop(get.par(x$state$parameters, "b"))

  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = M, sigma = P))

  ## Compute log priors.
  p2 <- x$prior(c("b" = g, get.par(x$state$parameters, "tau2")))
  qbetaprop <- dmvnorm(g, mean = M, sigma = P, log = TRUE)

  ## Compute fitted values.        
  x$state$fitted.values <- x$fit.fun(x$X, g)

  ## Set up new predictor.
  eta$gamma <- eta$gamma + x$state$fitted.values

  ## Compute new log likelihood.
  pibetaprop <- sum((eta$lambda + eta$gamma) * y[, "status"] - exp(eta$gamma) * int, na.rm = TRUE)

  ## Compute weights.
  weights <- exp(eta$gamma) * int

  ## Compute score.
  score <- y[, "status"] - exp(eta$gamma) * int

  ## Compute working observations.
  z <- eta$gamma + 1 / weights * score

  ## Compute reduced residuals.
  e <- z - eta2
  xbin.fun(x$binning$sorted.index, weights, e, x$weights, x$rres, x$binning$order)

  ## Compute mean and precision.
  XWX <- do.XWX(x$X, 1 / x$weights, x$imat)
  P2 <- if(x$fixed) {
    if(k < 2) {
      1 / (XWX)
    } else chol2inv(L2 <- chol(P02 <- XWX))
  } else {
    chol2inv(L2 <- chol(P02 <- XWX + S))
  }
  P2[P2 == Inf] <- 0
  M2 <- P2 %*% crossprod(x$X, x$rres)

  ## Get the log prior.
  qbeta <- dmvnorm(g0, mean = M2, sigma = P2, log = TRUE)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp)) {
    if(!x$fixed & is.null(x$sp)) {
      tau2 <- NULL
      for(j in seq_along(x$S)) {
        a <- x$rank[j] / 2 + x$a
        b <- 0.5 * crossprod(g, x$S[[j]]) %*% g + x$b
        tau2 <- c(tau2, 1 / rgamma(1, a, b))
      }
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    }
  }

  x$state$parameters <- set.par(x$state$parameters, g, "b")

  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

  return(x$state)
}


################################
## Survival helper functions. ##
################################
Surv2 <- function(..., obs = NULL, subdivisions = 100)
{
  require("survival")
  rval <- cbind(as.matrix(Surv(...)), "obs" = obs)
  class(rval) <- c("matrix", "Surv2")
  rval
}


rSurvTime2 <- function (lambda, x, cens_fct, upper = 1000, ..., file = NULL,
  subdivisions = 1000) 
{
  if(!is.matrix(x)) 
    x <- cbind(x)
  time <- rep(NA, nrow(x))
  Lambda <- function(lambda, x, time) {
    integrate(lambda, 0, time, x = x, subdivisions = subdivisions)$value
  }
  InvLambda <- function(Lambda, lambda, x) {
    negLogU <- -log(runif(1, 0, 1))
    rootfct <- function(time) {
      negLogU - Lambda(lambda, x, time)
    }
    return(uniroot(rootfct, interval = c(0, upper))$root)
  }
  for(i in 1:nrow(x)) {
    time[i] = InvLambda(Lambda, lambda, x[i, ])
  }
  time_event = cens_fct(time, ...)
  data = data.frame(time = time_event[, 1], event = time_event[, 2], x = x)
  names(data) <- gsub("x.", "x", names(data), fixed = TRUE)
  if(!is.null(file)) {
    save(data, file = file)
    invisible(data)
  } else {
    return(data)
  }
}


## Survival models transformer function.
surv.transform <- function(x, y, data, family,
  subdivisions = 100, timedependent = "lambda",
  timevar = NULL, idvar = NULL, is.cox = FALSE, alpha = 0.1, ...)
{
  rn <- names(y)
  y <- y[[rn]]
  ntd <- timedependent
  if(!all(ntd %in% names(x)))
    stop("the predictor names are different from family object names!")

  ## The basic setup.
  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, ...)
  
  ## Remove intercept if Cox.
  if(is.cox) {
    if(!is.null(x$lambda$smooth.construct$model.matrix)) {
      cn <- colnames(x$lambda$smooth.construct$model.matrix$X)
      if("(Intercept)" %in% cn)
        x$lambda$smooth.construct$model.matrix$X <- x$lambda$smooth.construct$model.matrix$X[, cn != "(Intercept)", drop = FALSE]
      if(ncol(x$lambda$smooth.construct$model.matrix$X) < 1) {
        x$lambda$smooth.construct$model.matrix <- NULL
        x$lambda$terms <- drop.terms.bamlss(x$lambda$terms, pterms = FALSE, keep.intercept = FALSE)
      } else {
        x$lambda$smooth.construct$model.matrix$term <- gsub("(Intercept)+", "",
          x$lambda$smooth.construct$model.matrix$term, fixed = TRUE)
        x$lambda$smooth.construct$model.matrix$state$parameters <- x$lambda$smooth.construct$model.matrix$state$parameters[-1]
        attr(x$lambda$terms, "intercept") <- 0
      }
    }
  }
  if(any("alpha" %in% names(x))) {
    x$alpha$smooth.construct$model.matrix$state$parameters[1] <- alpha
    x$alpha$smooth.construct$model.matrix$state$fitted.values <- x$alpha$smooth.construct$model.matrix$X %*% x$alpha$smooth.construct$model.matrix$state$parameters
  }

  ## Create the time grid.
  if(!inherits(y, c("Surv", "Surv2")))
    stop("the response is not a 'Surv' object, use function Surv() or Surv2()!")

  if(is.null(dim(y))) {
    y <- cbind("time" = y, "status" = rep(1, length(y)))
    class(y) <- "Surv"
  }
  if(length(i <- which(y[, "time"] == 0)))
    y[i, "time"] <- 1e-08 * abs(diff(range(y[, "time"])))

  grid <- function(upper, length){
    seq(from = 0, to = upper, length = length)
  }

  take <- NULL
  if(!is.null(idvar)) {
    if(!(idvar %in% names(attr(x, "model.frame"))))
      stop(paste("variable", idvar, "not in supplied data set!"))
    y2 <- cbind(y[, "time"], attr(x, "model.frame")[[idvar]])
    colnames(y2) <- c("time", idvar)
    take <- !duplicated(y2)
    y2 <- y2[take, , drop = FALSE]
    nobs <- nrow(y2)
    grid <- lapply(y2[, "time"], grid, length = subdivisions)
  } else {
    nobs <- nrow(y)
    grid <- lapply(y[, "time"], grid, length = subdivisions)
  }
  width <- rep(NA, nobs)
  for(i in 1:nobs)
    width[i] <- grid[[i]][2]
  attr(y, "width") <- width
  attr(y, "subdivisions") <- subdivisions
  attr(y, "grid") <- grid
  attr(y, "nobs") <- nobs
  attr(y, "take") <- take
  yname <- response.name(formula(as.Formula(x[[ntd[1]]]$formula, rhs = FALSE)))[1]
  if(is.null(timevar))
    timevar <- yname

  ## Assign time grid predict functions
  ## and create time dependant predictor.
  for(i in seq_along(ntd)) {
    if(has_pterms(x[[ntd[i]]]$terms)) {
      x[[ntd[i]]]$smooth.construct$model.matrix <- param_time_transform(x[[ntd[i]]]$smooth.construct$model.matrix,
        drop.terms.bamlss(x[[ntd[i]]]$terms, sterms = FALSE, keep.response = FALSE), data, grid, yname, timevar, take)
    }
    if(length(x[[ntd[i]]]$smooth.construct)) {
      for(j in names(x[[ntd[i]]]$smooth.construct)) {
        if(j != "model.matrix") {
          xterm <- x[[ntd[i]]]$smooth.construct[[j]]$term
          by <- if(x[[ntd[i]]]$smooth.construct[[j]]$by != "NA") x[[ntd[i]]]$smooth.construct[[j]]$by else NULL
          x[[ntd[i]]]$smooth.construct[[j]] <- sm_time_transform(x[[ntd[i]]]$smooth.construct[[j]],
            data[, unique(c(xterm, yname, by, timevar, idvar)), drop = FALSE], grid, yname, timevar, take)
        }
      }
    }
  }

  ## Assign index matrices for fast computation of integrals.
  for(j in ntd) {
    for(sj in seq_along(x[[j]]$smooth.construct)) {
      x[[j]]$smooth.construct[[sj]]$imat <- index_mat(x[[j]]$smooth.construct[[sj]]$X)
    }
  }

  y <- data.frame(y)
  names(y) <- rn

  family$p2d <- function(par, log = FALSE, ...) {
    x <- set.starting.values(x, par)
    eta <- get.eta(x)
    eta_timegrid <- 0
    for(sj in seq_along(x$lambda$smooth.construct)) {
      g <- get.state(x$lambda$smooth.construct[[sj]], "b")
      fit_timegrid <- x$lambda$smooth.construct[[sj]]$fit.fun_timegrid(g)
      x$lambda$smooth.construct[[sj]]$state$fitted_timegrid <- fit_timegrid
      eta_timegrid <- eta_timegrid + x$lambda$smooth.construct[[sj]]$fit.fun_timegrid(g)
    }
    eeta <- exp(eta_timegrid)
    int <- width * (0.5 * (eeta[, 1] + eeta[, subdivisions]) + apply(eeta[, 2:(subdivisions - 1)], 1, sum))
    d <- (eta$lambda + eta$gamma) * y[[1]][, "status"] - exp(eta$gamma) * int
    if(!log)
      d <- exp(d)
    return(d)
  }

  return(list("x" = x, "y" = y, "family" = family))
}


param_time_transform <- function(x, formula, data, grid, yname, timevar, take)
{
  if(!is.null(take))
    data <- data[take, , drop = FALSE]
  X <- Xn <- NULL
  for(j in names(data)) {
    if((!grepl("Surv(", j, fixed = TRUE) & !grepl("Surv2(", j, fixed = TRUE)) & (j != yname) & (j != timevar)) {
      df <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      names(df) <- j
      X <- if(is.null(X)) df else cbind(X, df)
      Xn <- c(Xn, j)
    }
  }
  if(!is.null(X))
    colnames(X) <- Xn
  X <- if(is.null(X)) data.frame(unlist(grid)) else cbind(X, unlist(grid))
  colnames(X)[ncol(X)] <- yname
  if(timevar != yname) {
    X <- cbind(X, unlist(grid))
    colnames(X)[ncol(X)] <- timevar
  }
  X <- model.matrix(formula, data = X)

  gdim <- c(length(grid), length(grid[[1]]))

  x$XT <- extract_XT(X, gdim[1], gdim[2])

  x$fit.fun_timegrid <- function(g) {
    if(is.null(g)) return(X)
    g <- get.par(g, "b")
    f <- drop(X %*% g)
    f <- matrix(f, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    f
  }
  x$state$fitted_timegrid <- x$fit.fun_timegrid(get.state(x, "b"))

  x
}

sm_time_transform <- function(x, data, grid, yname, timevar, take)
{
  if(!is.null(take))
    data <- data[take, , drop = FALSE]
  X <- NULL
  for(j in x$term) {
    if((j != yname) & (j != timevar)) {
      df <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      names(df) <- j
      X <- if(is.null(X)) df else cbind(X, df)
    }
  }
  if(!is.null(X))
    colnames(X) <- x$term[!(x$term %in% c(yname, timevar))]
  X <- if(is.null(X)) data.frame(unlist(grid)) else cbind(X, unlist(grid))
  colnames(X)[ncol(X)] <- yname
  if(timevar != yname) {
    X <- cbind(X, unlist(grid))
    colnames(X)[ncol(X)] <- timevar
  }
  if(x$by != "NA" & x$by != yname)
    X[[x$by]] <- rep(data[[x$by]], each = length(grid[[1]]))
  X <- PredictMat(x, X)
  gdim <- c(length(grid), length(grid[[1]]))

  x$XT <- extract_XT(X, gdim[1], gdim[2])

  x$grid.imat <- index_mat(X)
  ff <- make.fit.fun(x, type = 2)

  x$fit.fun_timegrid <- function(g) {
    if(is.null(g)) return(X)
    f <- ff(X, g, expand = FALSE)
    f <- matrix(f, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    f
  }
  x$state$fitted_timegrid <- x$fit.fun_timegrid(get.state(x, "b"))
  x$state$optimize <- FALSE

  x
}


## Extract the XT matrix.
extract_XT <- function(X, tnr, tnc)
{
  .Call("extract_XT", X, as.integer(tnr), as.integer(tnc))
}


## Survival integrals.
survint <- function(X, eta, width, gamma, eta2 = NULL, index = NULL)
{
  if(check <- is.null(eta2))
    eta2 <- as.numeric(0.0)
  int <- if(is.null(index)) {
    .Call("survint", X, eta, width, gamma, eta2, as.integer(check))
  } else {
    .Call("survint_index", X, eta, width, gamma, eta2, as.integer(check), index)
  }
  return(int)
}


## Survival probabilities.
bamlss.surv.prob <- function(object, newdata, type = c("probabilities", "link", "parameter"),
  FUN = function(x) { mean(x, na.rm = TRUE) }, time, subdivisions = 100, ...)
{
  if(is.null(newdata))
    newdata <- model.frame(object)
  if(length(type) > 1)
    type <- type[1]
  type <- match.arg(type)
  if(type != "probabilities") {
    object$family$predict <- NULL
    return(predict.bamlss(object, newdata, type, ...))
  }
  if(object$family$family != "cox")
    stop("object must be a cox-survival model!")
  if(missing(time))
    stop("please specify the time!")
  yname <- response.name(formula(as.Formula(object$x$lambda$formula, rhs = FALSE)))[1]
  timegrid <- rep(list(seq(0, time, length = subdivisions)), length = nrow(newdata))
  gdim <- c(length(timegrid), length(timegrid[[1]]))
  width <- timegrid[[1]][2]

  pred.setup <- predict.bamlss(object, newdata, type = "link", get.bamlss.predict.setup = TRUE, ...)
  enames <- pred.setup$enames

  pred_tc <- with(pred.setup, .predict.bamlss("gamma",
    object$x$gamma, samps, enames$gamma, intercept,
    nsamps, newdata, env))

  pred_td <- with(pred.setup, .predict.bamlss.surv.td("lambda",
    object$x$lambda$smooth.construct, samps, enames$lambda, intercept,
    nsamps, newdata, env, yname, timegrid,
    drop.terms.bamlss(object$x$lambda$terms, sterms = FALSE, keep.response = FALSE)))

  probs <- NULL
  for(i in 1:ncol(pred_td)) {
    eta <- matrix(pred_td[, i], nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    eeta <- exp(eta)
    int <- width * (0.5 * (eeta[, 1] + eeta[, subdivisions]) + apply(eeta[, 2:(subdivisions - 1)], 1, sum))
    probs <- cbind(probs, exp(-1 * exp(pred_tc[, i]) * int))
  }
return(probs)
  if(!is.null(FUN)) {
    if(!is.matrix(probs))
      probs <- matrix(probs, ncol = 1)
    probs <- apply(probs, 1, FUN, ...)
  }

  return(probs)
}


sm_Xtimegrid <- function(x, data, grid, yname)
{
  ff <- sm_time_transform(x, data, grid, yname, timevar = yname, take = NULL)$fit.fun_timegrid
  ff(NULL)
}

param_Xtimegrid <- function(formula, data, grid, yname)
{
  ff <- param_time_transform(list(), formula, data, grid, yname, timevar = yname, take = NULL)
  ff(NULL)
}


.predict.bamlss.surv.td <- function(id, x, samps, enames, intercept, nsamps, newdata, env,
  yname, grid, formula)
{
  snames <- colnames(samps)
  enames <- gsub("p.Intercept", "p.(Intercept)", enames, fixed = TRUE)
  has_intercept <- any(grepl(paste(id, "p", "(Intercept)", sep = "."), snames, fixed = TRUE))
  if(intercept & has_intercept)
    enames <- c("p.(Intercept)", enames)
  enames <- unique(enames)
  ec <- sapply(enames, function(x) {
    paste(strsplit(x, "")[[1]][1:2], collapse = "")
  })
  enames2 <- sapply(enames, function(x) {
    paste(strsplit(x, "")[[1]][-c(1:2)], collapse = "")
  })

  eta <- NULL
  if(length(i <- grep("p.", ec))) {
    for(j in enames2[i]) {
      if(j != "(Intercept)") {
        f <- as.formula(paste("~", if(has_intercept) "1" else "-1", "+", j), env = env)
        X <- param_Xtimegrid(f, newdata, grid, yname)
        if(has_intercept)
          X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
        sn <- snames[grep2(paste(id, "p", j, sep = "."), snames, fixed = TRUE)]
        eta <- if(is.null(eta)) {
          fitted_matrix(X, samps[, sn, drop = FALSE])
        } else {
          eta + fitted_matrix(X, samps[, sn, drop = FALSE])
        }
      } else {
        if(has_intercept) {
          sn <- snames[grep2(paste(id, "p", j, sep = "."), snames, fixed = TRUE)]
          eta <- eta + as.numeric(fitted_matrix(matrix(1, nrow = length(grid[[1]]) * nrow(newdata), ncol = 1),
            samps[, sn, drop = FALSE]))
        }
      }
    }
  }
  if(length(i <- grep("s.", ec))) {
    for(j in enames2[i]) {
      if(!inherits(x[[j]], "no.mgcv") & !inherits(x[[j]], "special")) {
        X <- sm_Xtimegrid(x[[j]], newdata, grid, yname)
        sn <- snames[grep2(paste(id, "s", j, sep = "."), snames, fixed = TRUE)]
        eta <- if(is.null(eta)) {
          fitted_matrix(X, samps[, sn, drop = FALSE])
        } else {
          eta + fitted_matrix(X, samps[, sn, drop = FALSE])
        }
      } else {
        stop("no predictions for special terms available yet!")
      }
    }
  }
 
  eta
}

