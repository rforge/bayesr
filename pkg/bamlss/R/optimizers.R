##################################################
## (1) Generic setup function for smooth terms. ##
##################################################
## For each s() term, an additional list() named "state" must be supplied,
## within the list of specifications returned by smooth.construct(). This list
## contains the following objects:
##
## - fitted.values: numeric, vector containing the current fitted values.
## - parameters: numeric, vector containing the current parameters.
##               Can also contain smoothing variances/parameters,
##               which should be named with "tau2", regression coefficients
##               should be named with "g1", "g2", "g3", ..., "g"k.
## - edf: numeric, the current degrees of freedom of the term.
## - optimize: NULL or logical, if NULL then per default the backfitting
##             algorithm is allowed to optimize possible variance
##             parameters within one update step of a term.
## - interval: numeric, if optimization is allowed, specifies the min and max
##             of the search space for optimizing variances. This can also
##             be supplied with the xt argument of s().
## - grid: integer, the grid length for searching variance parameters, if needed.
##         Can also be supplied within the xt argument in s().
##
## 4 additional functions must be supplied using the provided backfitting functions.
##
## - get.mu(X, gamma): computes the fitted values.
## - prior(parameters): computes the log prior using the parameters.
## - edf(x): computes degrees of freedom of the smooth term.
## - grad(score, parameters, ...): function that computes the gradient.
##
## NOTE: parametric effects are also added to the smooth term list with
##       appropriate penalty structure. The name of the object in the
##       list is "parametric", for later identifyng the pure parametric
##       modeled effects.
##
## bamlss.setup() sets up the basic structure, i.e., adds
## possible parametric terms to the smooth term list in x, also
## adds parametric terms of a random effect presentation of smooth
## terms to the "parametric" term. It calls the generic function
## smooth.bamlss(), which adds additional parts to the
## state list, as this could vary for special terms. A default
## method is provided.
bamlss.setup <- function(x, update = "iwls1", do.optim = NULL, criterion = c("AICc", "BIC", "AIC"), ...)
{
  if(!is.null(attr(x, "bamlss.setup"))) return(x)

  criterion <- match.arg(criterion)

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

  foo <- function(x) {
    if(!any(c("formula", "fake.formula", "response") %in% names(x))) {
      nx <- names(x)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(x)
      if(length(unique(nx)) < length(x)) nx <- 1:length(x)
      for(j in nx)
        x[[j]] <- foo(x[[j]])
    } else {
      if(!is.null(dim(x$X))) {
        if(nrow(x$X) > 0 & !is.na(mean(unlist(x$X), na.rm = TRUE))) {
          if(is.null(x$smooth)) x$smooth <- list()
          label <- if(is.null(colnames(x$X))) {
            paste("b", 1:ncol(x$X), sep = "", collapse = "+")
          } else paste(colnames(x$X), collapse = "+")
          x$smooth[["parametric"]] <- list(
            "X" = x$X,
            "S" = list(diag(0, ncol(x$X))),
            "rank" = ncol(x$X),
            "term" = label,
            "label" = label,
            "bs.dim" = ncol(x$X),
            "fixed" = TRUE,
            "is.parametric" = TRUE,
            "by" = "NA",
            "xt" = list("xbin" = x$binning)
          )
          class(x$smooth[["parametric"]]) <- c(class(x$smooth[["parametric"]]), "no.mgcv")
          x$sterms <- c(x$strems, "parametric")
          x$X <- NULL
        }
      }
      if(length(x$smooth)) {
        for(j in seq_along(x$smooth)) {
          x$smooth[[j]] <- smooth.bamlss(x$smooth[[j]])
          if(is.null(x$smooth[[j]]$update))
            x$smooth[[j]]$update <- eval(parse(text = paste("update", update, sep = "_")))
          if(is.null(x$smooth[[j]]$state$optimize)) {
            if(is.null(do.optim))
              x$smooth[[j]]$state$optimize <- TRUE
            else
              x$smooth[[j]]$state$optimize <- do.optim
          }
          if(!is.null(x$smooth[[j]]$rank))
            x$smooth[[j]]$rank <- as.numeric(x$smooth[[j]]$rank)
          x$smooth[[j]]$criterion <- criterion
          if(!is.null(x$smooth[[j]]$Xf)) {
            x$smooth[[j]]$Xfcn <- paste(paste(paste(x$smooth[[j]]$term, collapse = "."),
              "Xf", sep = "."), 1:ncol(x$smooth[[j]]$Xf), sep = ".")
            colnames(x$smooth[[j]]$Xf) <- x$smooth[[j]]$Xfcn
            if(is.null(x$smooth[["parametric"]])) {
              label <- paste(x$smooth[[j]]$Xfcn, collapse = "+")
              x$smooth[["parametric"]] <- list(
                "X" = x$smooth[[j]]$Xf,
                "S" = list(diag(0, ncol(x$Xf))),
                "rank" = ncol(x$smooth[[j]]$Xf),
                "term" = label,
                "label" = label,
                "bs.dim" = ncol(x$smooth[[j]]$Xf),
                "fixed" = TRUE,
                "is.parametric" = TRUE,
                "by" = "NA"
              )
              x$sterms <- c(x$strems, "parametric")
            } else {
              cn <- colnames(x$smooth[["parametric"]]$X)
              x$smooth[["parametric"]]$X <- cbind(x$smooth[["parametric"]]$X, x$smooth[[j]]$Xf)
              x$smooth[["parametric"]]$S <- list(diag(0, ncol(x$smooth[["parametric"]]$X)))
              x$smooth[["parametric"]]$bs.dim <- list(diag(0, ncol(x$smooth[["parametric"]]$X)))
              cn <- gsub("Intercept.", "Intercept", gsub("X.", "", c(cn , x$smooth[[j]]$Xfcn), fixed = TRUE))
              x$smooth[["parametric"]]$s.colnames <- colnames(x$smooth[["parametric"]]$X) <- cn
            }
          }
        }
      }
    }
    x
  }

  x <- foo(x)

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

  attr(x, "bamlss.setup") <- TRUE

  x
}


## Generic additional setup function for smooth terms.
smooth.bamlss <- function(x, ...) {
  UseMethod("smooth.bamlss")
}

## Simple extractor function.
get.state <- function(x, what = NULL) {
  if(is.null(what)) return(x$state)
  if(what %in% c("tau2", "tau", "lambda")) {
    p <- x$state$parameters
    return(p[grep("tau", names(p))])
  } else {
    if(what %in% c("g", "gamma")) {
      p <- x$state$parameters
      return(p[grep("g", names(p))])
    } else return(x$state[[what]])
  }
}

get.par <- function(x, what = NULL) {
  if(is.null(what)) return(x)
  if(what %in% c("tau2", "tau", "lambda")) {
    return(x[grep("tau", names(x))])
  } else {
    if(what %in% c("g", "gamma")) {
      return(x[grep("g", names(x))])
    } else return(x[what])
  }
}

set.par <- function(x, replacement, what) {
  if(what %in% c("tau2", "tau", "lambda")) {
    x[grep("tau", names(x))] <- replacement
  } else {
    if(what %in% c("g", "gamma")) {
      x[grep("g", names(x))] <- replacement
    } else x[what] <- replacement
  }
  x
}

## The default method.
smooth.bamlss.default <- function(x, ...)
{
  if(is.null(x$xbin.ind) & !is.null(x$xt$xbin)) {
    ind <- as.vector(apply(x$X, 1, function(x2) {
      if(!is.character(x2) & !is.factor(x2)) {
        if(!is.logical(x$xt$xbin))
          rval <- format(x2, digits = x$xt$xbin, nsmall = x$xt$xbin)
        else rval <- sprintf("%.48f", x2)
      } else rval <- x2
      paste(rval, collapse = ",", sep = "")
    }))
    uind <- unique(ind)
    x$xbin.take <- !duplicated(ind)
    x$xbin.ind <- rep(NA, nrow(x$X))
    xbin.uind <- seq_along(uind)
    for(ii in xbin.uind)
      x$xbin.ind[ind == uind[ii]] <- ii
    x$xbin.order <- order(x$xbin.ind)
    x$xbin.k <- length(xbin.uind)
    x$xbin.sind <- x$xbin.ind[x$xbin.order]
    x$X <- x$X[x$xbin.take, , drop = FALSE]
    x$before <- TRUE
  }
  if(is.null(x$before)) x$before <- FALSE
  if(!is.null(x$xbin.ind) & !x$before) {
    x$X <- x$X[x$xbin.take, , drop = FALSE]
  }
  if(is.null(x$xbin.ind)) {
    x$xbin.k <- nrow(x$X)
    x$xbin.ind <- 1:x$xbin.k
    x$xbin.sind <- 1:x$xbin.k
    x$xbin.order <- 1:x$xbin.k
  }
  x$nobs <- length(x$xbin.ind)
  x$weights <- rep(0, length = x$xbin.k)
  x$rres <- rep(0, length = x$xbin.k)
  state <- if(is.null(x$xt$state)) list() else x$xt$state
  if(is.null(x$fixed))
    x$fixed <- if(!is.null(x$fx)) x$fx[1] else FALSE
  if(!x$fixed & is.null(state$interval))
    state$interval <- if(is.null(x$xt$interval)) tau2interval(x) else x$xt$interval
  state$grid <- if(is.null(x$xt$grid)) 40 else x$xt$grid
  ntau2 <- length(x$S)
  if(length(ntau2) < 1) {
    if(x$fixed) {
      x$sp <- 1e+20
      ntau2 <- 1
      x$S <- list(diag(ncol(x$X)))
    } else {
      x$sp <- NULL
    }
  }
  if(!is.null(x$sp)) {
    x$sp <- rep(x$sp, length.out = ntau2)
    for(j in seq_along(x$sp))
      if(x$sp[j] == 0) x$sp[j] <- .Machine$double.eps^0.5
    x$fxsp <- TRUE
  } else x$fxsp <- FALSE
  if(is.null(state$parameters)) {
    state$parameters <- rep(0, ncol(x$X))
    names(state$parameters) <- paste("g", 1:length(state$parameters), sep = "")
    if(is.null(x$is.parametric)) {
      if(ntau2 > 0) {
        tau2 <- if(is.null(x$sp)) {
          if(x$fixed) {
            rep(1e+20, length.out = ntau2)
          } else {
            rep(if(!is.null(x$xt$tau2)) {
              x$xt$tau2
            } else {
              if(!is.null(x$xt$lambda)) 1 / x$xt$lambda else 100
            }, length.out = ntau2)
          }
        } else rep(x$sp, length.out = ntau2)
        names(tau2) <- paste("tau2", 1:ntau2, sep = "")
        state$parameters <- c(state$parameters, tau2)
        x$a <- if(is.null(x$xt$a)) 1e-04 else x$xt$a
        x$b <- if(is.null(x$xt$b)) 1e-04 else x$xt$b
      }
    }
  }
  if(is.null(x$get.mu) | !is.function(x$get.mu)) {
    x$get.mu <- function(X, b, expand = TRUE) {
      if(!is.null(names(b)))
        b <- get.par(b, "gamma")
      f <- drop(X %*% b)
      if(!is.null(x$xbin.ind) & expand)
        f <- f[x$xbin.ind]
      return(f)
    }
  }
  if(!is.null(x$xt$prior))
    x$prior <- x$xt$prior
  if(is.null(x$prior) | !is.function(x$prior)) {
    x$prior <- function(parameters) {
      gamma <- parameters[grep("g", names(parameters))]
      tau2 <-  parameters[grep("tau", names(parameters))]
      if(x$fixed | !length(tau2)) {
        lp <- sum(dnorm(gamma, sd = 1000, log = TRUE))
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
  if(is.null(x$edf)) {
    x$edf <- function(x) {
      tau2 <- get.state(x, "tau2")
      if(x$fixed | !length(tau2)) return(ncol(x$X))
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
  if(is.null(x$grad)) {
    x$grad <- function(score, parameters, full = TRUE) {
      gamma <- get.par(parameters, "g")
      tau2 <-  get.par(parameters, "tau2")
      grad2 <- NULL
      if(x$fixed | !length(tau2)) {
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
      grad <- if(!is.null(x$xbin.ind)) {
        drop(crossprod(x$X[x$xbin.ind, , drop = FALSE], score)) + c(grad, grad2)
      } else drop(crossprod(x$X, score)) + c(grad, grad2)
      return(grad)
    }
  } else {
    if(!is.function(x$grad))
      x$grad <- NULL
  }
  ng <- length(get.par(state$parameters, "g"))
  x$lower <- c(rep(-Inf, ng),
    if(is.list(state$interval)) {
      unlist(sapply(state$interval, function(x) { x[1] }))
    } else state$interval[1])
  x$upper <- c(rep(Inf, ng),
    if(is.list(state$interval)) {
      unlist(sapply(state$interval, function(x) { x[2] }))
    } else state$interval[2])
  names(x$lower) <- names(x$upper) <- names(state$parameters)[1:length(x$upper)]
  if(!is.null(x$sp)) {
    if(length(x$sp) < 1)
      x$sp <- NULL
    if(is.logical(x$sp))
      x[["sp"]] <- NULL
  }
  state$fitted.values <- x$get.mu(x$X, get.par(state$parameters, "g"))
  x$state <- state
  if(!is.null(x$xt$optimize))
    x$state$optimize <- x$xt$optimize
  x$state$edf <- x$edf(x)

  x
}


## Function to find tau2 interval according to the
## effective degrees of freedom
tau2interval <- function(x, lower = .Machine$double.eps^0.25, upper = 1e+10) {
  XX <- crossprod(x$X)
  if(length(x$S) < 2) {
    objfun <- function(tau2, value) {
      df <- sum(diag(matrix_inv(XX + if(x$fixed) 0 else 1 / tau2 * x$S[[1]]) %*% XX))
      return((value - df)^2)
    }
    le <- try(optimize(objfun, c(lower, upper), value = 1)$minimum, silent = TRUE)
    ri <- try(optimize(objfun, c(lower, upper), value = ncol(x$X))$minimum, silent = TRUE)
    if(inherits(le, "try-error")) le <- 0.1
    if(inherits(ri, "try-error")) ri <- 1e+04
    return(c(le, ri))
  } else {
    return(rep(list(c(lower, upper)), length.out = length(x$S)))
  }
}


get.eta <- function(x)
{
  nx <- names(x)
  np <- length(nx)
  eta <- vector(mode = "list", length = np)
  names(eta) <- nx
  for(j in 1:np) {
    eta[[j]] <- 0
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      g <- get.state(x[[nx[j]]]$smooth[[sj]], "g")
      eta[[j]] <- eta[[j]] + as.numeric(x[[nx[j]]]$smooth[[sj]]$get.mu(x[[nx[j]]]$smooth[[sj]]$X, g))
    }
  }
  eta
}

get.edf <- function(x)
{
  nx <- names(x)
  np <- length(nx)
  edf <- 0
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      edf <- edf + x[[nx[j]]]$smooth[[sj]]$edf(x[[nx[j]]]$smooth[[sj]])
    }
  }
  edf
}


## Formatting for printing.
fmt <- function(x, width = 8, digits = 2) {
  txt <- formatC(round(x, digits), format = "f", digits = digits , width = width)
  if(nchar(txt) > width) {
    txt <- strsplit(txt, "")[[1]]
    txt <- paste(txt[1:width], collapse = "", sep = "")
  }
  txt
}

bfit0 <- function(x, criterion = c("AICc", "BIC", "AIC"),
  eps = .Machine$double.eps^0.25, maxit = 400, outer = FALSE, inner = FALSE,
  verbose = TRUE, digits = 4, ...)
{
  criterion <- match.arg(criterion)

  x <- bamlss.setup(x, criterion = criterion, ...)

  family <- attr(x, "family")
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  criterion <- match.arg(criterion)

  np <- length(nx)
  response <- attr(x, "response.vec")
  nobs <- length(response)
  eta <- get.eta(x)
  
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
        eta[[id]] <- eta[[id]] - fitted(x[[sj]]$state) + fitted(p.state)
        x[[sj]]$state <- p.state
      }
      eps0 <- do.call("cbind", eta)
      eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
      if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1
      iter <- iter + 1
    }
    return(list("x" = x, "eta" = eta, "edf" = edf))
  }

  ## Backfitting main function.
  backfit <- function(x, eta, verbose = TRUE) {
    edf <- get.edf(x)
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
            eta[[nx[j]]] <- eta[[nx[j]]] - fitted(x[[nx[j]]]$smooth[[sj]]$state) + fitted(p.state)

            x[[nx[j]]]$smooth[[sj]]$state <- p.state
          }
        }
      }

      eps0 <- do.call("cbind", eta)
      eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
      if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

      peta <- family$map2par(eta)

      if(verbose) {
        IC <- get.ic(family, response, peta, edf, nobs, criterion)
        cat("\r")
        vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
          " logLik ", fmt(family$loglik(response, peta), width = 8, digits = digits),
          " edf ", fmt(edf, width = 6, digits = digits),
          " eps ", fmt(eps0, width = 6, digits = digits + 2),
          " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
        cat(vtxt)

        if(.Platform$OS.type != "unix") flush.console()
      }

      iter <- iter + 1
    }

    IC <- get.ic(family, response, peta, edf, nobs, criterion)

    if(verbose) {
      cat("\r")
      vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
        " logLik ", fmt(family$loglik(response, peta), width = 8, digits = digits),
        " edf ", fmt(edf, width = 6, digits = digits),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
      cat(vtxt)
      if(.Platform$OS.type != "unix") flush.console()
      cat("\n")
    }

    if(iter == maxit)
      warning("the backfitting algorithm did not converge, please check argument eps and maxit!")

    return(list("x" = x, "eta" = eta, "ic" = IC))
  }

  bf <- backfit(x, eta, verbose = verbose)
  x <- bf$x; eta <- bf$eta
  rm(bf)

  return(x)
}


## Extract information criteria.
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


## Function to create full parameter vector.
make_par <- function(x, type = 1) {
  family <- attr(x, "family")
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  np <- length(nx)
  par <- lower <- upper <- NULL
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      tpar <- x[[nx[j]]]$smooth[[sj]]$state$parameters
      tlower <- x[[nx[j]]]$smooth[[sj]]$lower
      tupper <- x[[nx[j]]]$smooth[[sj]]$upper
      g <- get.par(tpar, "g")
      npar <- paste(paste(nx[j], "h1", x[[nx[j]]]$smooth[[sj]]$label, sep = ":"), 1:length(g), sep = ".")
      if(length(tau2 <- get.par(tpar, "tau2"))) {
        npar <- c(npar, paste(nx[j], "h1", paste(x[[nx[j]]]$smooth[[sj]]$label,
          paste("tau2", 1:length(tau2), sep = ""), sep = "."), sep = ":"))
      }
      names(tpar) <- names(tlower) <- names(tupper) <- if(type < 2) {
        paste("p", j, ".t", sj, ".", names(tpar), sep = "")
      } else npar
      par <- c(par, tpar)
      lower <- c(lower, tlower)
      upper <- c(upper, tupper)
    }
  }
  return(list("par" = par, "lower" = lower, "upper" = upper))
}


## Backfitting updating functions.
update_iwls1 <- function(x, family, response, eta, id, ...)
{
  args <- list(...)

  peta <- family$map2par(eta)

  if(is.null(args$weights)) {
    ## Compute weights.
    weights <- family$weights[[id]](response, peta)
  } else weights <- args$weights

  if(is.null(args$z)) {
    ## Score.
    score <- family$score[[id]](response, peta)

    ## Compute working observations.
    z <- eta[[id]] + 1 / weights * score
  } else z <- args$z

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - fitted(x$state)

  ## Compute reduced residuals.
  e <- z - eta[[id]]
  xbin.fun(x$xbin.sind, weights, e, x$weights, x$rres, x$xbin.order)

  ## Compute mean and precision.
  XWX <- crossprod(x$X, x$X * x$weights)
  if(!x$state$optimize | x$fixed | !is.null(x$sp)) {
    if(x$fixed) {
      P <- matrix_inv(XWX)
    } else {
      S <- 0
      tau2 <- get.state(x, "tau2")
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      P <- matrix_inv(XWX + S)
    }
    x$state$parameters <- set.par(x$state$parameters, drop(P %*% crossprod(x$X, x$rres)), "g")
  } else {
    args <- list(...)
    edf0 <- args$edf - x$state$edf
    eta2 <- eta
    objfun <- function(tau2, ...) {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      P <- matrix_inv(XWX + S)
      if(inherits(P, "try-error")) return(NA)
      g <- drop(P %*% crossprod(x$X, x$rres))
      if(any(is.na(g)) | any(g %in% c(-Inf, Inf))) g <- rep(0, length(g))
      fit <- x$get.mu(x$X, g)
      edf <- sum(diag(P %*% XWX))
      if(!is.null(x$xt$center)) {
        if(x$xt$center) edf <- edf - 1
      }
      eta2[[id]] <- eta2[[id]] + fit
      IC <- get.ic(family, response, family$map2par(eta2), edf0 + edf, length(z), x$criterion)
      return(IC)
    }
    if(length(get.state(x, "tau2")) < 2) {
      tau2 <- try(optimize(objfun, interval = x$state$interval)$minimum, silent = TRUE)
      if(inherits(tau2, "try-error"))
        tau2 <- optimize2(objfun, interval = x$state$interval, grid = x$state$grid)$minimum
      x$state$parameters <- set.par(x$state$parameters, if(!length(tau2)) x$interval[1] else tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$lower))
      opt <- try(optim(get.state(x, "tau2"), fn = objfun, method = "L-BFGS-B",
        lower = x$lower[i], upper = x$upper[i]), silent = TRUE)
      if(!inherits(opt, "try-error"))
        x$state$parameters <- set.par(x$state$parameters, opt$par, "tau2")
    }
    S <- 0
    tau2 <- get.state(x, "tau2")
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    P <- matrix_inv(XWX + S)
    x$state$parameters <- set.par(x$state$parameters, drop(P %*% crossprod(x$X, x$rres)), "g")
  }

  ## Compute fitted values.
  g <- get.state(x, "g")
  if(any(is.na(g)) | any(g %in% c(-Inf, Inf)))
    x$state$parameters <- set.par(x$state$parameters, rep(0, length(x$state$g)), "g")
  x$state$fitted.values <- x$get.mu(x$X, get.state(x, "g"))
  x$state$edf <- sum(diag(P %*% XWX))
  if(!is.null(x$xt$center)) {
    if(x$xt$center) x$state$edf <- x$state$edf - 1
  }

  return(x$state)
}


## Updating based on optim.
update_optim1 <- function(x, family, response, eta, id, ...)
{
  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - fitted(x$state)
  eta2 <- eta

  if(is.null(x$state$XX)) x$state$XX <- crossprod(x$X)
  if(!x$fixed) {
    args <- list(...)
    edf0 <- args$edf - x$state$edf
  }

  tpar <- x$state$parameters

  ## Objective for regression coefficients.
  objfun <- function(gamma, tau2 = NULL) {
    tpar <- set.par(tpar, gamma, "g")
    if(!is.null(tau2) & !x$fixed)
      tpar <- set.par(tpar, tau2, "tau2")
    eta2[[id]] <- eta[[id]] + x$get.mu(x$X, tpar)
    ll <- family$loglik(response, family$map2par(eta2))
    lp <- x$prior(tpar)
    -1 * (ll + lp)
  }

  ## Gradient function.
  grad <- if(!is.null(family$score[[id]]) & is.function(x$grad)) {
    function(gamma, tau2 = NULL) {
      tpar <- set.par(tpar, gamma, "g")
      if(!is.null(tau2) & !x$fixed)
        tpar <- set.par(tpar, tau2, "tau2")
      eta2[[id]] <- eta[[id]] + x$get.mu(x$X, tpar)
      peta <- family$map2par(eta2)
      score <- drop(family$score[[id]](response, peta))
      grad <- x$grad(score, tpar, full = FALSE)
      return(drop(-1 * grad))
    }
  } else NULL

  if(!x$fixed & x$state$optimize & is.null(x$sp)) {
    objfun2 <- function(tau2) {
      tpar <- set.par(tpar, tau2, "tau2")
      suppressWarnings(opt <- try(optim(get.par(tpar, "g"), fn = objfun, gr = grad,
        method = "BFGS", control = list(), tau2 = tau2), silent = TRUE))
      if(!inherits(opt, "try-error")) {
        tpar <- set.par(tpar, opt$par, "g")
        x$state$fitted.values <- x$get.mu(x$X, tpar)
      }
      x$state$parameters <- tpar
      edf <- x$edf(x)
      eta2[[id]] <- eta[[id]] + fitted(x$state)
      IC <- get.ic(family, response, family$map2par(eta2), edf0 + edf, length(eta2[[id]]), x$criterion)
      IC
    }
    if(length(get.state(x, "tau2")) < 2) {
      tau2 <- try(optimize(objfun2, interval = x$state$interval)$minimum, silent = TRUE)
      if(inherits(tau2, "try-error"))
        tau2 <- optimize2(objfun2, interval = x$state$interval, grid = x$state$grid)$minimum
      tpar <- set.par(tpar, if(!length(tau2)) x$interval[1] else tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$lower))
      opt <- try(optim(get.state(x, "tau2"), fn = objfun2, method = "L-BFGS-B",
        lower = x$lower[i], upper = x$upper[i]), silent = TRUE)
      if(!inherits(opt, "try-error"))
        tpar <- set.par(tpar, opt$par, "tau2")
    }
  }

  suppressWarnings(opt <- try(optim(get.par(tpar, "g"), fn = objfun, gr = grad,
    method = "BFGS", control = list(), tau2 = get.par(tpar, "tau2")), silent = TRUE))

  if(!inherits(opt, "try-error")) {
    tpar <- set.par(tpar, opt$par, "g")
    x$state$fitted.values <- x$get.mu(x$X, tpar)
    x$state$parameters <- tpar
  }

  if(!x$fixed)
    x$state$edf <- x$edf(x)

  return(x$state)
}


set.all.par <- function(par, x)
{
  nx <- names(x)
  np <- length(x)
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      tpar <- par[grep(paste("p", j, ".t", sj, ".", sep = ""), names(par), fixed = TRUE)]
      x[[nx[j]]]$smooth[[sj]]$state$parameters <- tpar
      x[[nx[j]]]$smooth[[sj]]$state$fitted.values <- x[[nx[j]]]$smooth[[sj]]$get.mu(x[[nx[j]]]$smooth[[sj]]$X,
        get.par(tpar, "g"))
    }
  }
  return(x)
}


log_posterior <- function(par, x, verbose = TRUE, criterion = "AICc", digits = 3, scale = 1)
{
  eta <- get.eta(x)
  nx <- names(eta)
  np <- length(eta)
  family <- attr(x, "family")
  lprior <- 0
  edf <- 0
  start <- 1
  for(j in 1:np) {
    eta[[nx[j]]] <- 0
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      end <- length(x[[nx[j]]]$smooth[[sj]]$state$parameters) + start - 1
      i <- start:end
      tpar <- par[i]
      names(tpar) <- names(x[[nx[j]]]$smooth[[sj]]$state$parameters)
      start <- max(i) + 1
      if(any((tau2 <- get.par(tpar, "tau2")) < 0)) {
        tau2[tau2 < 0] <- rep(.Machine$double.eps^0.25, sum(tau2 < 0))
        tpar <- set.par(tpar, tau2, "tau2")
      }
      x[[nx[j]]]$smooth[[sj]]$state$parameters <- tpar
      x[[nx[j]]]$smooth[[sj]]$state$fitted.values <- x[[nx[j]]]$smooth[[sj]]$get.mu(x[[nx[j]]]$smooth[[sj]]$X,
        get.par(tpar, "g"))
      eta[[nx[j]]] <- eta[[nx[j]]] + fitted(x[[nx[j]]]$smooth[[sj]]$state)
      lprior <- lprior + x[[nx[j]]]$smooth[[sj]]$prior(tpar)
      edf <- edf + x[[nx[j]]]$smooth[[sj]]$edf(x[[nx[j]]]$smooth[[sj]])
    }
  }
  ll <- family$loglik(attr(x, "response.vec"), family$map2par(eta))
  lp <- as.numeric(ll + lprior)

  if(verbose) {
    nobs <- if(is.null(dim(attr(x, "response.vec")))) {
      length(attr(x, "response.vec"))
    } else nrow(attr(x, "response.vec"))
    IC <- get.ic(family, attr(x, "response.vec"), family$map2par(eta), edf, nobs, criterion)
    cat("\r")
    vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
      " logPost ", fmt(lp, width = 8, digits = digits),
      " edf ", fmt(edf, width = 6, digits = digits), sep = "")
    cat(vtxt)
    if(.Platform$OS.type != "unix") flush.console()
  }

  return(lp * scale)
}

grad_posterior <- function(par, x, ...)
{
  eta <- get.eta(x)
  nx <- names(eta)
  np <- length(eta)
  family <- attr(x, "family")
  grad <- NULL
  for(j in 1:np) {
    eta[[nx[j]]] <- 0
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      tpar <- par[grep(paste("p", j, ".t", sj, ".", sep = ""), names(par), fixed = TRUE)]
      if(any((tau2 <- get.par(tpar, "tau2")) < 0)) {
        tau2[tau2 < 0] <- rep(.Machine$double.eps^0.25, sum(tau2 < 0))
        tpar <- set.par(tpar, tau2, "tau2")
      }
      x[[nx[j]]]$smooth[[sj]]$state$parameters <- tpar
      x[[nx[j]]]$smooth[[sj]]$state$fitted.values <- x[[nx[j]]]$smooth[[sj]]$get.mu(x[[nx[j]]]$smooth[[sj]]$X,
        get.par(tpar, "g"))
      eta[[nx[j]]] <- eta[[nx[j]]] + fitted(x[[nx[j]]]$smooth[[sj]]$state)
    }
  }
  for(j in 1:np) {
    score <- family$score[[nx[j]]](attr(x, "response.vec"), family$map2par(eta))
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      tpar <- par[grep(paste("p", j, ".t", sj, ".", sep = ""), names(par), fixed = TRUE)]
      tgrad <- x[[nx[j]]]$smooth[[sj]]$grad(score, tpar)
      grad <- c(grad, tgrad)
    }
  }
  return(grad)
}


opt0 <- function(x, criterion = c("AICc", "BIC", "AIC"), verbose = TRUE, digits = 3, ...)
{
  x <- bamlss.setup(x, ...)

  criterion = match.arg(criterion)

  par <- make_par(x)
  family <- attr(x, "family")

  opt <- optim(par$par, fn = log_posterior, gr = if(!is.null(family$score)) grad_posterior else NULL,
    x = x, method = "L-BFGS-B", lower = par$lower, upper = par$upper, verbose = verbose,
    criterion = criterion, digits = digits, control = list(fnscale = -1), hessian = TRUE)
 
  if(verbose) cat("\n")

  x <- set.all.par(opt$par, x)
  attr(x, "hessian") <- opt$hessian
  attr(x, "converged") <- opt$convergence < 1

  x
}


xbin.fun <- function(ind, weights, e, xweights, xrres, oind)
{
  .Call("xbin_fun", as.integer(ind), as.numeric(weights), 
    as.numeric(e), as.numeric(xweights), as.numeric(xrres),
    as.integer(oind))
  invisible(NULL)
}

