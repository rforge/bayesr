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
    "names" = c("lambda", "mu"),
    "links" = c(lambda = "log", mu = "identity"),
    "transform" = function(x, ...) { surv.transform(x, globalgrid = FALSE, ...) },
    "engine" = cox.engine,
    "loglik" = function(y, eta, ...) {
      n <- attr(y, "subdivisions")
      eeta <- exp(eta_Surv_timegrid)
      int <- attr(y, "width") * (0.5 * (eeta[, 1] + eeta[, n]) + apply(eeta[, 2:(n - 1)], 1, sum))
      ll <- (eta$lambda + eta$mu) * y[, "status"] - exp(eta$mu) * int
      sum(ll)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}

## (2) New model fitting engines.
cox.engine <- function(x, ...)
{
  stacker(x, optimizer = cox.mode, sampler = cox.mcmc, ...)
}

## (3) Posterior mode estimation.
cox.mode <- function(x, nu = 1, eps = .Machine$double.eps^0.25, maxit = 400,
  verbose = TRUE, digits = 4, ...)
{
  ## Names of parameters/predictors.
  nx <- names(x)
  
  ## Compute additive predictors.
  eta <- get.eta(x)

  ## For the time dependent part, compute
  ## predictor based on the time grid.
  eta_timegrid <- attr(x, "eta_Surv_timegrid")

  ## Compute current equivalent degrees of freedom.
  edf <- get.edf(x)

  ## The reponse, a 'Surv' object.
  response <- attr(x, "response.vec")

  ## Number of observations.
  nobs <- nrow(response)

  ## Number of subdivions used for the time grid.
  sub <- attr(response, "subdivisions")

  ## The interval width from subdivisons.
  width <- attr(response, "width")

  ## Save the hessian.
  hessian <- list()

  ## Start the backfitting algorithm.
  eps0 <- eps + 1; iter <- 1
  while(eps0 > eps & iter < maxit) {
    eta0 <- eta

    ########################################
    ## Cycle through time-dependent part. ##
    ########################################
    for(sj in seq_along(x$lambda$smooth)) {
      ## The time-dependent design matrix for the grid.
      X <- x$lambda$smooth[[sj]]$get.mu_timegrid(NULL)

      ## Timegrid lambda.
      eeta <- exp(eta_timegrid)

      ## Dimension of the design matrix.
      xdim <- dim(x$lambda$smooth[[sj]]$X)

      ## Compute the gradient.
      tint <- vector("list", xdim[2])
      for(i in 1:ncol(X)) {
        tint[[i]] <- matrix(X[, i], nrow = nrow(eeta), ncol = ncol(eeta), byrow = TRUE)
        tint[[i]] <- tint[[i]] * eeta
        tint[[i]] <- width * (0.5 * (tint[[i]][, 1] + tint[[i]][, sub]) + apply(tint[[i]][, 2:(sub - 1)], 1, sum))
      }
      tint <- sapply(tint, cbind)
      tint <- tint * exp(eta$mu)
      int <- apply(tint, 2, sum)
      xgrad <- drop(t(response[, "status"]) %*% x$lambda$smooth[[sj]]$X - int)
      xgrad <- xgrad + x$lambda$smooth[[sj]]$grad(score = NULL, x$lambda$smooth[[sj]]$state$parameters, full = FALSE)
      
      ## Compute the hessian.
      tint <- vector("list", nobs)
      xhess <- matrix(0, nrow = xdim[2], ncol = xdim[2])
      for(i in 1:nobs) {
        forward <- sub * (i - 1)
        tint[[i]] <- matrix(0, nrow = xdim[2], ncol = xdim[2])
        for(j in 1:sub) {
          MAT <- X[j + forward,] %o% X[j + forward,] * eeta[i, j]
          if(j == 1 || j == sub) {
            tint[[i]] <- tint[[i]] + 0.5 * MAT
          } else {
            tint[[i]] <- tint[[i]] + MAT
          }
        }
        tint[[i]] <- tint[[i]] * width[i]
        xhess <- xhess + exp(eta$mu[i]) * tint[[i]]
      }
      xhess <- xhess + x$lambda$smooth[[sj]]$hess(score = NULL, x$lambda$smooth[[sj]]$state$parameters, full = FALSE)

      ## Compute the inverse of the hessian.
      hessian[[paste("p", 1, ".t", sj, ".", sep = "")]] <- xhess
      Sigma <- matrix_inv(xhess)

      ## Update regression coefficients.
      g <- get.state(x$lambda$smooth[[sj]], "gamma")
      g2 <- drop(g + nu * Sigma %*% xgrad)
      names(g2) <- names(g)
      x$lambda$smooth[[sj]]$state$parameters <- set.par(x$lambda$smooth[[sj]]$state$parameters, g2, "g")

      ## Update additive predictors.
      fit_timegrid <- x$lambda$smooth[[sj]]$get.mu_timegrid(g2)
      eta_timegrid <- eta_timegrid - x$lambda$smooth[[sj]]$state$fitted_timegrid + fit_timegrid
      x$lambda$smooth[[sj]]$state$fitted_timegrid <- fit_timegrid

      fit <- x$lambda$smooth[[sj]]$get.mu(x$lambda$smooth[[sj]]$X, g2)
      eta$lambda <- eta$lambda - fitted(x$lambda$smooth[[sj]]$state) + fit
      x$lambda$smooth[[sj]]$state$fitted.values <- fit

      ## Save Sigma.
      x$lambda$smooth[[sj]]$state$hessian <- Sigma
    }

    ###########################################
    ## Actual integral of survivor function. ##
    ###########################################
    eeta <- exp(eta_timegrid)
    int <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))

    ##########################################
    ## Cycle through time-independent part. ##
    ##########################################
    for(sj in seq_along(x$mu$smooth)) {
      ## Compute weights.
      weights <- exp(eta$mu) * int

      ## Compute score.
      score <- response[, "status"] - exp(eta$mu) * int

      ## Compute working observations.
      z <- eta$mu + 1 / weights * score

      ## Compute partial predictor.
      eta$mu <- eta$mu - fitted(x$mu$smooth[[sj]]$state)

      ## Compute reduced residuals.
      e <- z - eta$mu
      xbin.fun(x$mu$smooth[[sj]]$xbin.sind, weights, e,
        x$mu$smooth[[sj]]$weights, x$mu$smooth[[sj]]$rres,
        x$mu$smooth[[sj]]$xbin.order)

      ## Compute mean and precision.
      XWX <- crossprod(x$mu$smooth[[sj]]$X, x$mu$smooth[[sj]]$X * x$mu$smooth[[sj]]$weights)
      if(x$mu$smooth[[sj]]$fixed) {
        hessian[[paste("p", 2, ".t", sj, ".", sep = "")]] <- XWX
        P <- matrix_inv(XWX)
      } else {
        S <- 0
        tau2 <- get.state(x$mu$smooth[[sj]], "tau2")
        for(j in seq_along(x$mu$smooth[[sj]]$S))
          S <- S + 1 / tau2[j] * x$mu$smooth[[sj]]$S[[j]]
        hessian[[paste("p", 2, ".t", sj, ".", sep = "")]] <- XWX + S
        P <- matrix_inv(XWX + S)
      }
      g <- drop(P %*% crossprod(x$mu$smooth[[sj]]$X, x$mu$smooth[[sj]]$rres))
      x$mu$smooth[[sj]]$state$parameters <- set.par(x$mu$smooth[[sj]]$state$parameters, g, "g")

      ## Compute fitted values.
      if(any(is.na(g)) | any(g %in% c(-Inf, Inf))) {
        x$mu$smooth[[sj]]$state$parameters <- set.par(x$mu$smooth[[sj]]$state$parameters,
          rep(0, length(x$mu$smooth[[sj]]$state$g)), "g")
      }
      x$mu$smooth[[sj]]$state$fitted.values <- x$mu$smooth[[sj]]$get.mu(x$mu$smooth[[sj]]$X,
        get.state(x$mu$smooth[[sj]], "g"))
      x$mu$smooth[[sj]]$state$edf <- sum(diag(P %*% XWX))
      if(!is.null(x$mu$smooth[[sj]]$xt$center)) {
        if(x$mu$smooth[[sj]]$xt$center) x$mu$smooth[[sj]]$state$edf <- x$mu$smooth[[sj]]$state$edf - 1
      }

      ## Update additive predictor.
      eta$mu <- eta$mu + fitted(x$mu$smooth[[sj]]$state)
    }

    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

    if(verbose) {
      ll <- sum((eta$lambda + eta$mu) * response[, "status"] - exp(eta$mu) * int, na.rm = TRUE)
      cat("\r")
      vtxt <- paste(
        "logLik ", fmt(ll, width = 8, digits = digits),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)), sep = ""
      )
      cat(vtxt)
      if(.Platform$OS.type != "unix") flush.console()
    }

    iter <- iter + 1
  }

  ## Assign time-dependant predictor to .GlobalEnv.
  eta_Surv_timegrid  <<- eta_timegrid

  if(iter == maxit)
    warning("the backfitting algorithm did not converge, please check argument eps and maxit!")

  if(verbose) cat("\n")

  nh <- NULL
  for(j in names(hessian))
    nh <- c(nh, paste(j, paste("g", 1:ncol(hessian[[j]]), sep = ""), sep = ""))
  require("Matrix")
  hessian <- -1 * as.matrix(do.call("bdiag", hessian))
  rownames(hessian) <- colnames(hessian) <- nh
  attr(x, "hessian") <- hessian

  return(x)
}

## (4) The MCMC sampling engine.
cox.mcmc <- function(x, n.iter = 1200, burnin = 200, thin = 1, ...)
{
  null.sampler(x, ...)
}


################################
## Survival helper functions. ##
################################
Surv2 <- function(..., obs = NULL, subdivisions = 100)
{
  require("survival")
  rval <- cbind(as.matrix(Surv(...)), "obs" = obs)
  class(rval) <- c("Surv")
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


integrate2 <- function(f, a, b, ...)
{
  nodes <- c(
    -0.999713726773441,
    -0.998491950639596,
    -0.996295134733125,
    -0.993124937037443,
    -0.988984395242991,
    -0.983877540706057,
    -0.977809358486919,
    -0.970785775763707,
    -0.962813654255816,
    -0.953900782925493,
    -0.944055870136256,
    -0.93328853504308,
    -0.921609298145334,
    -0.90902957098253,
    -0.895561644970728,
    -0.881218679385019,
    -0.866014688497165,
    -0.849964527879591,
    -0.833083879888401,
    -0.815389238339177,
    -0.796897892390315,
    -0.777627909649496,
    -0.757598118519707,
    -0.73682808980202,
    -0.715338117573057,
    -0.693149199355802,
    -0.670283015603141,
    -0.64676190851413,
    -0.622608860203708,
    -0.597847470247179,
    -0.572501932621382,
    -0.546597012065094,
    -0.520158019881763,
    -0.493210789208192,
    -0.465781649773359,
    -0.437897402172032,
    -0.409585291678302,
    -0.38087298162463,
    -0.351788526372422,
    -0.32236034390053,
    -0.292617188038472,
    -0.262588120371503,
    -0.232302481844974,
    -0.201789864095736,
    -0.171080080538604,
    -0.140203137236114,
    -0.109189203580061,
    -0.0780685828134366,
    -0.0468716824215914,
    -0.0156289844215428,
    0.0156289844215435,
    0.0468716824215913,
    0.0780685828134369,
    0.109189203580061,
    0.140203137236114,
    0.171080080538603,
    0.201789864095736,
    0.232302481844974,
    0.262588120371504,
    0.292617188038472,
    0.322360343900529,
    0.351788526372422,
    0.38087298162463,
    0.409585291678302,
    0.437897402172032,
    0.465781649773358,
    0.493210789208191,
    0.520158019881763,
    0.546597012065094,
    0.572501932621381,
    0.597847470247179,
    0.622608860203708,
    0.64676190851413,
    0.670283015603142,
    0.693149199355802,
    0.715338117573057,
    0.736828089802021,
    0.757598118519707,
    0.777627909649496,
    0.796897892390315,
    0.815389238339176,
    0.833083879888401,
    0.849964527879592,
    0.866014688497164,
    0.881218679385019,
    0.895561644970727,
    0.90902957098253,
    0.921609298145334,
    0.93328853504308,
    0.944055870136256,
    0.953900782925492,
    0.962813654255816,
    0.970785775763707,
    0.977809358486919,
    0.983877540706057,
    0.988984395242992,
    0.993124937037444,
    0.996295134733125,
    0.998491950639596,
    0.999713726773442
  )

  weights <- c(
    0.000734634490505196,
    0.00170939265351828,
    0.00268392537155396,
    0.00365596120132663,
    0.00462445006342169,
    0.00558842800386553,
    0.00654694845084563,
    0.00749907325546386,
    0.00844387146966807,
    0.00938041965369529,
    0.0103078025748694,
    0.0112251140231855,
    0.0121314576629791,
    0.0130259478929708,
    0.0139077107037196,
    0.0147758845274415,
    0.015629621077546,
    0.0164680861761455,
    0.0172904605683233,
    0.0180959407221274,
    0.0188837396133761,
    0.019653087494435,
    0.0204032326462095,
    0.0211334421125275,
    0.0218430024162464,
    0.0225312202563371,
    0.0231974231852537,
    0.0238409602659681,
    0.0244612027079576,
    0.0250575444815794,
    0.0256294029102085,
    0.0261762192395456,
    0.0266974591835708,
    0.0271926134465774,
    0.0276611982207921,
    0.0281027556591007,
    0.0285168543223954,
    0.0289030896011252,
    0.0292610841106386,
    0.0295904880599137,
    0.0298909795933315,
    0.0301622651051686,
    0.0304040795264544,
    0.0306161865839809,
    0.030798379031153,
    0.0309504788504913,
    0.0310723374275658,
    0.03116383569621,
    0.0312248842548494,
    0.0312554234538637,
    0.0312554234538625,
    0.0312248842548495,
    0.0311638356962099,
    0.0310723374275674,
    0.030950478850491,
    0.0307983790311533,
    0.0306161865839802,
    0.0304040795264556,
    0.0301622651051682,
    0.0298909795933326,
    0.0295904880599136,
    0.0292610841106384,
    0.0289030896011253,
    0.0285168543223953,
    0.0281027556591007,
    0.0276611982207923,
    0.0271926134465769,
    0.0266974591835716,
    0.0261762192395464,
    0.025629402910208,
    0.0250575444815787,
    0.024461202707957,
    0.0238409602659674,
    0.0231974231852531,
    0.0225312202563353,
    0.0218430024162472,
    0.0211334421125276,
    0.0204032326462092,
    0.0196530874944354,
    0.0188837396133751,
    0.0180959407221291,
    0.0172904605683237,
    0.0164680861761457,
    0.0156296210775464,
    0.0147758845274413,
    0.0139077107037179,
    0.0130259478929714,
    0.0121314576629803,
    0.011225114023186,
    0.0103078025748689,
    0.00938041965369477,
    0.00844387146966944,
    0.00749907325546371,
    0.00654694845084634,
    0.00558842800386576,
    0.00462445006342231,
    0.00365596120132682,
    0.00268392537155354,
    0.00170939265351811,
    0.000734634490505862
  )

  fn <- function(x) {
    (b - a) / 2 * f((b - a) / 2 * x + (a + b) / 2, ...)
  }

  sum(fn(nodes) * weights)
}


integrate3 <- function(f, a, b, m = 2, n = 40, nx = 1000, ret.fun = FALSE, plot = FALSE) {
  require("splines")
  xl <- a
  xu <- b
  xr <- xu - xl
  xl <- xl - xr * 0.001
  xu <- xu + xr * 0.001
  dx <- (xu - xl)/(n - 1)
  k <- seq(xl - dx * (m + 1), xu + dx * (m + 1), length = n + 2 * m + 2)
  x <- seq(a, b, length = nx)
  X <- spline.des(k, x, m + 2, rep(0, nx))$design
  X1 <- spline.des(k, x, m + 2, rep(1, nx))$design
  X <- X[, -1]
  X1 <- X1[, -1]
  y <- f(x)
  bf <- lm.fit(X1, y)
  y2 <- drop(X %*% coef(bf))
  f2 <- splinefun(x, y2)
  if(plot) {
    x2 <- abs(f(x))
    i <- which(x2 <= quantile(x2, prob = 0.01))
    par(mfrow = c(2, 1))
    plot(x, f2(x), type = "l", ylab = "int(f(x))")
    abline(v = x[i], lty = 2)
    plot(x, f(x), type = "l", col = "red", lwd = 2, ylab = "f(x)")
    lines(x, drop(X1 %*% coef(bf)), col = "black")
    abline(h = 0, lty = 2)
    abline(v = x[i], lty = 2)
  }
  if(ret.fun) {
    return(f2)
  } else {
    return(f2(b) - f2(a))
  }
}


#####################
## 2nd Cox version ##
#####################
cox2.bamlss <- function(links = c(lambda = "identity", mu = "identity"), ...)
{
  require("survival")
  rval <- list(
    "family" = "cox",
    "names" = c("lambda", "mu"),
    "links" = parse.links(links, c(lambda = "log", mu = "identity"), ...),
    "transform" = surv.transform,
    "loglik" = function(y, eta, ...) {
      n <- attr(y, "subdivisions")
      eeta <- exp(eta_Surv_timegrid)
      int <- attr(y, "width") * (0.5 * (eeta[, 1] + eeta[, n]) + apply(eeta[, 2:(n - 1)], 1, sum))
      ll <- (eta$lambda + eta$mu) * y[, "status"] - exp(eta$mu) * int
      sum(ll)
    },
    "score" = list(
      "mu" = function(y, eta, ...) {
        n <- attr(y, "subdivisions")
        eeta <- exp(eta_Surv_timegrid)
        int <- attr(y, "width") * (0.5 * (eeta[, 1] + eeta[, n]) + apply(eeta[, 2:(n - 1)], 1, sum))
        y[, "status"] - exp(eta$mu) * int
      }
    ),
    "weights" = list(
      "mu" = function(y, eta, ...) {
        n <- attr(y, "subdivisions")
        eeta <- exp(eta_Surv_timegrid)
        int <- attr(y, "width") * (0.5 * (eeta[, 1] + eeta[, n]) + apply(eeta[, 2:(n - 1)], 1, sum))
        exp(eta$mu) * int
      }
    ),
    "gradient" = list(
      "lambda" = function(g, y, eta, x, ...) {
        n <- attr(y, "subdivisions")
        X <- x$get.mu_timegrid(NULL)
        eeta <- eta_Surv_timegrid + x$get.mu_timegrid(g)
        eeta <- exp(eeta)
        dummy <- vector("list", ncol(x$X))
        for(i in 1:ncol(x$X)) {
          dummy[[i]] <- matrix(X[, i], nrow = nrow(eeta), ncol = ncol(eeta), byrow = TRUE)
          dummy[[i]] <- dummy[[i]] * eeta
          dummy[[i]] <- attr(y, "width") * (0.5 * (dummy[[i]][, 1] + dummy[[i]][, n]) + apply(dummy[[i]][, 2:(n - 1)], 1, sum))
        }
        dummy <- sapply(dummy, cbind)
        dummy <- dummy * exp(eta$mu)
        int <- apply(dummy, 2, sum)
        xgrad <- drop(t(y[, "status"]) %*% x$X - int)
        return(xgrad)
      }
    ),
    "hessian" = list(
      "lambda" = function(g, y, eta, x, ...) {
        n <- attr(y, "subdivisions")
        X <- x$get.mu_timegrid(NULL)
        eeta <- eta_Surv_timegrid + x$get.mu_timegrid(g)
        eeta <- exp(eeta)
        nobs <- nrow(y)
        dummy <- vector("list", nobs)
        width <- attr(y, "width")
        xhess <- matrix(0, ncol = ncol(x$X), nrow = ncol(x$X))
        for(i in 1:nobs) {
          forward <- n * (i - 1)
          dummy[[i]] <- matrix(0, ncol = ncol(X), nrow = ncol(X))
          for(j in 1:n) {
            MAT <- X[j + forward,] %o% X[j + forward,] * eeta[i, j]
            if(j == 1 || j == n){
              dummy[[i]] <- dummy[[i]] + 0.5 * MAT
            } else {
              dummy[[i]] <- dummy[[i]] + MAT
            }
          }
          dummy[[i]] <- dummy[[i]] * width[i]
          xhess <- xhess + exp(eta$mu[i]) * dummy[[i]]
        }
        return(xhess)
      }
    )
  )
  class(rval) <- "family.bamlss"
  rval
}

