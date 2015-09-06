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
##               should be named with "b1", "b2", "b3", ..., "b"k.
## - edf: numeric, the current degrees of freedom of the term.
## - do.optim: NULL or logical, if NULL then per default the backfitting
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
## - edf(x, weights): computes degrees of freedom of the smooth term.
## - grad(score, parameters, ...): function that computes the gradient.
##
## NOTE: model.matrix effects are also added to the smooth term list with
##       appropriate penalty structure. The name of the object in the
##       list is "model.matrix", for later identifyng the pure model.matrix
##       modeled effects.
##
## bamlss.engine.setup() sets up the basic structure, i.e., adds
## possible model.matrix terms to the smooth term list in x, also
## adds model.matrix terms of a random effect presentation of smooth
## terms to the "model.matrix" term. It calls the generic function
## bamlss.engine.setup.smooth(), which adds additional parts to the
## state list, as this could vary for special terms. A default
## method is provided.
bamlss.engine.setup <- function(x, update = "iwls",
  do.optim = NULL, criterion = c("AICc", "BIC", "AIC"),
  nu = 0.1, start = NULL, df = NULL, ...)
{
  if(!is.null(attr(x, "bamlss.engine.setup"))) return(x)

  criterion <- match.arg(criterion)

  foo <- function(x, id = NULL) {
    if(!any(c("formula", "fake.formula") %in% names(x))) {
      for(j in names(x))
        x[[j]] <- foo(x[[j]], id = c(id, j))
    } else {
      if(is.null(id)) id <- ""
      if(!is.null(dim(x$model.matrix))) {
        if(nrow(x$model.matrix) > 0 & !is.na(mean(unlist(x$model.matrix), na.rm = TRUE))) {
          if(is.null(x$smooth.construct)) x$smooth.construct <- list()
          label <- if(is.null(colnames(x$model.matrix))) {
            paste("b", 1:ncol(x$model.matrix), sep = "", collapse = "+")
          } else paste(colnames(x$model.matrix), collapse = "+")
          x$smooth.construct <- c(list("model.matrix" = list(
            "X" = x$model.matrix,
            "S" = list(diag(0, ncol(x$model.matrix))),
            "rank" = ncol(x$model.matrix),
            "term" = label,
            "label" = label,
            "bs.dim" = ncol(x$model.matrix),
            "fixed" = TRUE,
            "is.model.matrix" = TRUE,
            "by" = "NA",
            "xt" = list("binning" = x$binning)
          )), x$smooth.construct)
          class(x$smooth.construct[["model.matrix"]]) <- c(class(x$smooth.construct[["model.matrix"]]),
            "no.mgcv", "model.matrix")
          x$model.matrix <- NULL
        }
      }
      if(length(x$smooth.construct)) {
        for(j in seq_along(x$smooth.construct)) {
          x$smooth.construct[[j]] <- bamlss.engine.setup.smooth(x$smooth.construct[[j]])
          x$smooth.construct[[j]] <- assign.df(x$smooth.construct[[j]], df)
          if(!is.null(x$smooth.construct[[j]]$xt$update))
            x$smooth.construct[[j]]$update <- x$smooth.construct[[j]]$xt$update
          if(is.null(x$smooth.construct[[j]]$update)) {
            if(is.character(update)) {
              if(!grepl("bfit_", update))
                update <- paste("bfit", update, sep = "_")
              update <- eval(parse(text = update))
            }
            x$smooth.construct[[j]]$update <- update
          }
          if(is.null(x$smooth.construct[[j]]$state$do.optim)) {
            if(is.null(do.optim))
              x$smooth.construct[[j]]$state$do.optim <- TRUE
            else
              x$smooth.construct[[j]]$state$do.optim <- do.optim
          }
          if(!is.null(x$smooth.construct[[j]]$rank))
            x$smooth.construct[[j]]$rank <- as.numeric(x$smooth.construct[[j]]$rank)
          x$smooth.construct[[j]]$criterion <- criterion
          x$smooth.construct[[j]]$nu <- nu
          if(!is.null(x$smooth.construct[[j]]$Xf)) {
            x$smooth.construct[[j]]$Xfcn <- paste(paste(paste(x$smooth.construct[[j]]$term, collapse = "."),
              "Xf", sep = "."), 1:ncol(x$smooth.construct[[j]]$Xf), sep = ".")
            colnames(x$smooth.construct[[j]]$Xf) <- x$smooth.construct[[j]]$Xfcn
            if(is.null(x$smooth.construct[["model.matrix"]])) {
              label <- paste(x$smooth.construct[[j]]$Xfcn, collapse = "+")
              x$smooth.construct[["model.matrix"]] <- list(
                "X" = x$smooth.construct[[j]]$Xf,
                "S" = list(diag(0, ncol(x$Xf))),
                "rank" = ncol(x$smooth.construct[[j]]$Xf),
                "term" = label,
                "label" = label,
                "bs.dim" = ncol(x$smooth.construct[[j]]$Xf),
                "fixed" = TRUE,
                "is.model.matrix" = TRUE,
                "by" = "NA",
                "nu" = nu
              )
              x$smooth.construct <- c(x$smooth.construct, "model.matrix")
            } else {
              x$smooth.construct[["model.matrix"]]$X <- cbind(x$smooth.construct[["model.matrix"]]$X, x$smooth.construct[[j]]$Xf)
              x$smooth.construct[["model.matrix"]]$S <- list(diag(0, ncol(x$smooth.construct[["model.matrix"]]$X)))
              x$smooth.construct[["model.matrix"]]$bs.dim <- list(diag(0, ncol(x$smooth.construct[["model.matrix"]]$X)))
            }
          }
        }
      }
    }
    x
  }

  x <- foo(x)
  attr(x, "bamlss.engine.setup") <- TRUE

  x
}


## Generic additional setup function for smooth terms.
bamlss.engine.setup.smooth <- function(x, ...) {
  UseMethod("bamlss.engine.setup.smooth")
}

## Simple extractor function.
get.state <- function(x, what = NULL) {
  if(is.null(what)) return(x$state)
  if(what %in% c("par", "parameters")) {
    return(x$state$parameters)
  } else {
    if(what %in% c("tau2", "tau", "lambda")) {
      p <- x$state$parameters
      return(p[grep("tau", names(p))])
    } else {
      if(what %in% "b") {
        p <- x$state$parameters
        return(p[!grepl("tau", names(p))])
      } else return(x$state[[what]])
    }
  }
}

get.par <- function(x, what = NULL) {
  if(is.null(what)) return(x)
  if(what %in% c("tau2", "tau", "lambda")) {
    return(x[grep("tau", names(x))])
  } else {
    if(what %in% "b") {
      return(x[!grepl("tau", names(x))])
    } else return(x[what])
  }
}

set.par <- function(x, replacement, what) {
  if(what %in% c("tau2", "tau", "lambda")) {
    x[grep("tau", names(x))] <- replacement
  } else {
    if(what %in% "b") {
      x[!grepl("tau", names(x))] <- replacement
    } else x[what] <- replacement
  }
  x
}

## The default method.
bamlss.engine.setup.smooth.default <- function(x, ...)
{
  if(inherits(x, "special"))
    return(x)
  before <- NULL
  if(is.null(x$binning) & !is.null(x$xt$binning)) {
    x$binning <- match.index(x$X)
    x$binning$order <- order(x$binning$match.index)
    x$binning$sorted.index <- x$binning$match.index[x$binning$order]
    x$X <- x$X[x$binning$nodups, , drop = FALSE]
    before <- TRUE
  }
  if(is.null(before)) before <- FALSE
  if(!is.null(x$binning) & !before) {
    x$X <- x$X[x$binning$nodups, , drop = FALSE]
  }
  if(is.null(x$binning)) {
    nr <- nrow(x$X)
    x$binning <- list(
      "match.index" = 1:nr,
      "nodups" = 1:nr,
      "order" = 1:nr,
      "sorted.index" = 1:nr
    )
  }
  x$nobs <- length(x$binning$match.index)
  k <- length(x$binning$nodups)
  x$weights <- rep(0, length = k)
  x$rres <- rep(0, length = k)
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
    names(state$parameters) <- if(is.null(colnames(x$X))) {
      paste("b", 1:length(state$parameters), sep = "")
    } else colnames(x$X)
    if(is.null(x$is.model.matrix)) {
      if(ntau2 > 0) {
        tau2 <- if(is.null(x$sp)) {
          if(x$fixed) {
            rep(1e+20, length.out = ntau2)
          } else {
            rep(if(!is.null(x$xt$tau2)) {
              x$xt$tau2
            } else {
              if(!is.null(x$xt$lambda)) 1 / x$xt$lambda else 1000
            }, length.out = ntau2)
          }
        } else rep(x$sp, length.out = ntau2)
        names(tau2) <- paste("tau2", 1:ntau2, sep = "")
        state$parameters <- c(state$parameters, tau2)
      }
    }
  }
  if((ntau2 > 0) & !any(grepl("tau2", names(state$parameters))) & is.null(x$is.model.matrix)) {
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
  }
  x$a <- if(is.null(x$xt$a)) 1e-04 else x$xt$a
  x$b <- if(is.null(x$xt$b)) 1e-04 else x$xt$b
  if(!is.null(x$xt$prior))
    x$prior <- x$xt$prior
  if(is.null(x$prior) | !is.function(x$prior)) {
    x$prior <- function(parameters) {
      gamma <- parameters[!grepl("tau", names(parameters))]
      tau2 <-  parameters[grepl("tau", names(parameters))]
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
    x$edf <- function(x, type = 1) {
      if(type > 1) {
        if(!is.null(x$state$edf))
          return(x$state$edf)
      }
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
      return(edf)
    }
  }
  if(is.null(x$grad)) {
    x$grad <- function(score = NULL, parameters, full = TRUE) {
      gamma <- get.par(parameters, "b")
      tau2 <-  get.par(parameters, "tau2")
      grad2 <- NULL
      if(x$fixed | !length(tau2)) {
        grad <- 0
      } else {
        grad <- 0; grad2 <- NULL
        for(j in seq_along(tau2)) {
          gS <- x$S[[j]] %*% gamma
          grad <- grad + drop(-0.5 / tau2[j] * gS)
          if(full & !is.null(tau2[j])) {
            grad2 <- c(grad2, drop(-x$rank[j] / (2 * tau2[j]) - 1 / (2 * tau2[j]^2) * gS %*% gamma + (-x$a - 1) / tau2[j] + x$b / (tau2[j]^2)))
            x$X <- cbind(x$X, 0)
          }
        }
      }
      if(!is.null(score)) {
        grad <- if(!is.null(x$xbin.ind)) {
          drop(crossprod(x$X[x$xbin.ind, , drop = FALSE], score)) + c(grad, grad2)
        } else drop(crossprod(x$X, score)) + c(grad, grad2)
      } else grad <- c(grad, grad2)
      return(grad)
    }
  } else {
    if(!is.function(x$grad))
      x$grad <- NULL
  }
  if(is.null(x$hess)) {
    x$hess <- function(score = NULL, parameters, full = FALSE) {
      tau2 <- get.par(parameters, "tau2")
      if(x$fixed | !length(tau2)) {
        hx <- 0
      } else {
        hx <- 0
        for(j in seq_along(tau2)) {
          hx <- hx + (1 / tau2[j]) * x$S[[j]]
        }
      }
      return(hx)
    }
  } else {
    if(!is.function(x$hess))
      x$hess <- NULL
  }
  ng <- length(get.par(state$parameters, "b"))
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
  if(is.null(x$fit.fun))
    x$fit.fun <- make.fit.fun(x)
  state$fitted.values <- x$fit.fun(x$X, get.par(state$parameters, "b"))
  x$state <- state
  if(!is.null(x$xt$do.optim))
    x$state$do.optim <- x$xt$do.optim
  x$state$edf <- x$edf(x)
  if(!inherits(x, "model.matrix"))
    colnames(x$X) <- NULL
  x$added <- c("nobs", "weights", "rres", "state", "grid", "a", "b", "prior", "edf",
    "grad", "hess", "lower", "upper")

  x
}


## Function to find tau2 interval according to the
## effective degrees of freedom
tau2interval <- function(x, lower = .Machine$double.eps^0.25, upper = 1e+10) {
  XX <- crossprod(x$X)
  if(length(x$S) < 2) {
    objfun <- function(tau2, value) {
      df <- sum(diag(XX %*% matrix_inv(XX + if(x$fixed) 0 else 1 / tau2 * x$S[[1]])))
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


## Assign degrees of freedom.
assign.df <- function(x, df) {
  tau2 <- get.par(x$state$parameters, "tau2")
  if(x$fixed | !length(tau2))
    return(x)
  df <- if(is.null(x$xt$df)) df else x$xt$df
  if(is.null(df))
    df <- ceiling(length(get.par(x$state$parameters, "b")) / 2)
  if(length(tau2) > 1)
    return(x)
  if(df > ncol(x$X))
    df <- ncol(x$X)
  if(df < 1)
    df <- 1
  XX <- crossprod(x$X)
  objfun <- function(tau2) {
    edf <- sum(diag(XX %*% matrix_inv(XX + 1 / tau2 * x$S[[1]])))
    return((df - edf)^2)
  }
  tau2 <- try(optimize(objfun, c(.Machine$double.eps^0.25, 1e+10))$minimum, silent = TRUE)
  if(inherits(tau2, "try-error"))
    return(x)
  x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  x$state$edf <- df
  return(x)
}


get.eta <- function(x, expand = TRUE)
{
  nx <- names(x)
  np <- length(nx)
  eta <- vector(mode = "list", length = np)
  names(eta) <- nx
  for(j in 1:np) {
    eta[[j]] <- 0
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      fit <- x[[nx[j]]]$smooth.construct[[sj]]$fit.fun(x[[nx[j]]]$smooth.construct[[sj]]$X,
        x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, expand)
      eta[[j]] <- eta[[j]] + fit
    }
  }
  eta
}

get.edf <- function(x, type = 1)
{
  nx <- names(x)
  np <- length(nx)
  edf <- 0
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      edf <- edf + if(type < 2) {
        x[[nx[j]]]$smooth.construct[[sj]]$edf(x[[nx[j]]]$smooth.construct[[sj]])
      } else x[[nx[j]]]$smooth.construct[[sj]]$state$edf
    }
  }
  edf
}

get.log.prior <- function(x, type = 1)
{
  nx <- names(x)
  np <- length(nx)
  lp <- 0
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      lp <- lp + if(type < 2) {
        x[[nx[j]]]$smooth.construct[[sj]]$prior(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters)
      } else x[[nx[j]]]$smooth.construct[[sj]]$state$log.prior
    }
  }
  lp
}

get.all.par <- function(x)
{
  nx <- names(x)
  np <- length(nx)
  par <- list()
  for(i in nx) {
    par[[i]] <- list()
    if(!all(c("formula", "fake.formula") %in% names(x[[i]]))) {
      for(k in names(x[[i]])) {
        if(!is.null(x[[i]][[k]]$smooth.construct)) {
          par[[i]][[k]]$s <- list()
          for(j in names(x[[i]][[k]]$smooth.construct)) {
            if(j == "model.matrix") {
              par[[i]][[k]]$p <- x[[i]][[k]]$smooth.construct[[j]]$state$parameters
            } else {
              if(is.null(par[[i]][[k]]$s))
                par[[i]][[k]]$s <- list()
              par[[i]][[k]]$s[[j]] <- x[[i]][[k]]$smooth.construct[[j]]$state$parameters
            }
          }
        }
      }
    } else {
      if(!is.null(x[[i]]$smooth.construct)) {
        for(j in names(x[[i]]$smooth.construct)) {
          if(j == "model.matrix") {
            par[[i]]$p <- x[[i]]$smooth.construct[[j]]$state$parameters
          } else {
            if(is.null(par[[i]]$s))
              par[[i]]$s <- list()
            par[[i]]$s[[j]] <- x[[i]]$smooth.construct[[j]]$state$parameters
          }
        }
      }
    }
  }
  par
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

bfit <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
  criterion = c("AICc", "BIC", "AIC"), eps = .Machine$double.eps^0.25,
  maxit = 400, outer = FALSE, inner = FALSE,
  verbose = TRUE, digits = 4, ...)
{
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("design construct names mismatch with family names!")

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, ...)

  if(!is.null(start))
    x <- add.starting.values(x, start)

  criterion <- match.arg(criterion)
  np <- length(nx)
  nobs <- nrow(y)
  eta <- get.eta(x)
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  if(!is.null(weights))
    weights <- as.data.frame(weights)
  if(!is.null(offset)) {
    offset <- as.data.frame(offset)
    for(j in nx) {
      if(!is.null(offset[[j]]))
        eta[[j]] <- eta[[j]] + offset[[j]]
    }
  }

  inner_bf <- function(x, y, eta, family, edf, id, ...) {
    eps0 <- eps + 1; iter <- 1
    while(eps0 > eps & iter < maxit) {
      eta0 <- eta
      for(sj in seq_along(x)) {
        ## Get updated parameters.
        p.state <- x[[sj]]$update(x[[sj]], family, y, eta, id, edf = edf, ...)

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
  backfit <- function(verbose = TRUE) {
    eps0 <- eps + 1; iter <- 1
    edf <- get.edf(x, type = 2)
    while(eps0 > eps & iter < maxit) {
      eta0 <- eta
      ## Cycle through all parameters
      for(j in 1:np) {
        if(outer) {
          peta <- family$map2par(eta)

          ## Compute weights.
          hess <- family$hess[[nx[j]]](y, peta)

          ## Score.
          score <- family$score[[nx[j]]](y, peta)

          ## Compute working observations.
          z <- eta[[nx[j]]] + 1 / hess * score
        } else z <- hess <- NULL

        ## And all terms.
        if(inner) {
          tbf <- inner_bf(x[[nx[j]]]$smooth.construct, y, eta, family,
            edf = edf, id = nx[j], z = z, hess = hess, weights = weights[[nx[j]]])
          x[[nx[j]]]$smooth.construct <- tbf$x
          edf <- tbf$edf
          eta <- tbf$eta
          rm(tbf)
        } else {
          for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
            ## Get updated parameters.
            p.state <- x[[nx[j]]]$smooth.construct[[sj]]$update(x[[nx[j]]]$smooth.construct[[sj]],
              family, y, eta, nx[j], edf = edf, z = z, hess = hess, weights = weights[[nx[j]]])

            ## Compute equivalent degrees of freedom.
            edf <- edf - x[[nx[j]]]$smooth.construct[[sj]]$state$edf + p.state$edf

            ## Update predictor and smooth fit.
            eta[[nx[j]]] <- eta[[nx[j]]] - fitted(x[[nx[j]]]$smooth.construct[[sj]]$state) + fitted(p.state)

            x[[nx[j]]]$smooth.construct[[sj]]$state <- p.state
          }
        }
      }

      eps0 <- do.call("cbind", eta)
      eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
      if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

      peta <- family$map2par(eta)

      if(verbose) {
        IC <- get.ic(family, y, peta, edf, nobs, criterion)
        cat("\r")
        vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
          " logPost ", fmt(family$loglik(y, peta) + get.log.prior(x), width = 8, digits = digits),
          " edf ", fmt(edf, width = 6, digits = digits),
          " eps ", fmt(eps0, width = 6, digits = digits + 2),
          " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
        cat(vtxt)

        if(.Platform$OS.type != "unix") flush.console()
      }

      iter <- iter + 1
    }

    IC <- get.ic(family, y, peta, edf, nobs, criterion)
    logLik <- family$loglik(y, peta)
    logPost <- as.numeric(logLik + get.log.prior(x))

    if(verbose) {
      cat("\r")
      vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
        " logPost ", fmt(logPost, width = 8, digits = digits),
        " edf ", fmt(edf, width = 6, digits = digits),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
      cat(vtxt)
      if(.Platform$OS.type != "unix") flush.console()
      cat("\n")
    }

    if(iter == maxit)
      warning("the backfitting algorithm did not converge!")

    names(IC) <- criterion

    return(list("fitted.values" = eta, "parameters" = get.all.par(x), "edf" = edf,
      "logLik" = logLik, "logPost" = logPost, "IC" = IC))
  }

  backfit(verbose = verbose)
}


## Extract information criteria.
get.ic <- function(family, y, par, edf, n, type = c("AIC", "BIC", "AICc", "MP"), ...)
{
  type <- match.arg(type)
  ll <- family$loglik(y, par)
  pen <- switch(type,
    "AIC" = -2 * ll + 2 * edf,
    "BIC" = -2 * ll + edf * log(n),
    "AICc" = -2 * ll + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
    "MP" = -1 * (ll + edf)
  )
  return(pen)
}


## Function to create full parameter vector.
make_par <- function(x, type = 1, add.tau2 = FALSE) {
  family <- attr(x, "family")
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  np <- length(nx)
  par <- lower <- upper <- NULL
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      tpar <- x[[nx[j]]]$smooth.construct[[sj]]$state$parameters
      tlower <- x[[nx[j]]]$smooth.construct[[sj]]$lower
      tupper <- x[[nx[j]]]$smooth.construct[[sj]]$upper
      if(!add.tau2) {
        tlower <- tlower[!grepl("tau2", names(tlower))]
        tupper <- tupper[!grepl("tau2", names(tupper))]
        tpar <- tpar[!grepl("tau2", names(tpar))]
      }
      g <- get.par(tpar, "b")
      npar <- paste(paste(nx[j], "h1", x[[nx[j]]]$smooth.construct[[sj]]$label, sep = ":"), 1:length(g), sep = ".")
      if(length(tau2 <- get.par(tpar, "tau2"))) {
        npar <- c(npar, paste(nx[j], "h1", paste(x[[nx[j]]]$smooth.construct[[sj]]$label,
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
bfit_newton <- function(x, family, y, eta, id, ...)
{
  args <- list(...)

  eta[[id]] <- eta[[id]] - fitted(x$state)

  tau2 <- if(!x$fixed) get.par(x$state$parameters, "tau2") else NULL

  lp <- function(g) {
    eta[[id]] <- eta[[id]] + x$fit.fun(x$X, g)
    family$loglik(y, family$map2par(eta)) + x$prior(c(g, tau2))
  }

  if(is.null(family$gradient[[id]])) {
    gfun <- NULL
  } else {
    gfun <- list()
    gfun[[id]] <- function(g, y, eta, x, ...) {
      gg <- family$gradient[[id]](g, y, eta, x, ...)
      if(!is.null(x$grad)) {
        gg <- gg + x$grad(score = NULL, c(g, tau2), full = FALSE)
      }
      drop(gg)
    }
  }

  if(is.null(family$hessian[[id]])) {
    hfun <- NULL
  } else {
    hfun <- list()
    hfun[[id]] <- function(g, y, eta, x, ...) {
      hg <- family$hessian[[id]](g, y, eta, x, ...)
      if(!is.null(x$hess)) {
        hg <- hg + x$hess(score = NULL, c(g, tau2), full = FALSE)
      }
      hg
    }
  }

  g <- get.par(x$state$parameters, "b")
  nu <- if(is.null(x$nu)) 0.1 else x$nu

  g.grad <- grad(fun = lp, theta = g, id = id, prior = NULL,
    args = list("gradient" = gfun, "x" = x, "y" = y, "eta" = eta))

  g.hess <- hess(fun = lp, theta = g, id = id, prior = NULL,
    args = list("gradient" = gfun, "hessian" = hfun, "x" = x, "y" = y, "eta" = eta))

  Sigma <- matrix_inv(g.hess)

  g <- drop(g + nu * Sigma %*% g.grad)

  x$state$parameters <- set.par(x$state$parameters, g, "b")
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
  x$state$hessian <- Sigma

  return(x$state)
}


bfit_iwls <- function(x, family, y, eta, id, weights, ...)
{
  args <- list(...)

  peta <- family$map2par(eta)
  if(is.null(args$hess)) {
    ## Compute weights.
    hess <- family$hess[[id]](y, peta, id = id, ...)
  } else hess <- args$hess

  if(!is.null(weights))
    hess <- hess * weights

  if(is.null(args$z)) {
    ## Score.
    score <- family$score[[id]](y, peta, id = id, ...)

    ## Compute working observations.
    z <- eta[[id]] + 1 / hess * score
  } else z <- args$z

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - fitted(x$state)

  ## Compute reduced residuals.
  e <- z - eta[[id]]
  xbin.fun(x$binning$sorted.index, hess, e, x$weights, x$rres, x$binning$order)

  ## Compute mean and precision.
  XWX <- crossprod(x$X, x$X * x$weights)
  if(!x$state$do.optim | x$fixed | !is.null(x$sp)) {
    if(x$fixed) {
      P <- matrix_inv(XWX)
    } else {
      S <- 0
      tau2 <- get.state(x, "tau2")
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      P <- matrix_inv(XWX + S)
    }
    x$state$parameters <- set.par(x$state$parameters, drop(P %*% crossprod(x$X, x$rres)), "b")
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
      fit <- x$fit.fun(x$X, g)
      edf <- sum(diag(P %*% XWX))
      eta2[[id]] <- eta2[[id]] + fit
      IC <- get.ic(family, y, family$map2par(eta2), edf0 + edf, length(z), x$criterion, ...)
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
    x$state$parameters <- set.par(x$state$parameters, drop(P %*% crossprod(x$X, x$rres)), "b")
  }

  ## Compute fitted values.
  g <- get.state(x, "b")
  if(any(is.na(g)) | any(g %in% c(-Inf, Inf)))
    x$state$parameters <- set.par(x$state$parameters, rep(0, length(x$state$g)), "b")
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
  x$state$edf <- sum(diag(XWX %*% P))
  x$state$log.prior <- x$prior(x$state$parameters)

  return(x$state)
}


## Updating based on optim.
bfit_optim <- function(x, family, y, eta, id, ...)
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
    tpar <- set.par(tpar, gamma, "b")
    if(!is.null(tau2) & !x$fixed)
      tpar <- set.par(tpar, tau2, "tau2")
    eta2[[id]] <- eta[[id]] + x$fit.fun(x$X, tpar)
    ll <- family$loglik(y, family$map2par(eta2))
    lp <- x$prior(tpar)
    val <- -1 * (ll + lp)
    if(!is.finite(val)) val <- NA
    val
  }

  ## Gradient function.
  grad <- if(!is.null(family$score[[id]]) & is.function(x$grad)) {
    function(gamma, tau2 = NULL) {
      tpar <- set.par(tpar, gamma, "b")
      if(!is.null(tau2) & !x$fixed)
        tpar <- set.par(tpar, tau2, "tau2")
      eta2[[id]] <- eta[[id]] + x$fit.fun(x$X, tpar)
      peta <- family$map2par(eta2)
      score <- drop(family$score[[id]](y, peta))
      grad <- x$grad(score, tpar, full = FALSE)
      return(drop(-1 * grad))
    }
  } else NULL

  if(!x$fixed & x$state$do.optim & is.null(x$sp)) {
    objfun2 <- function(tau2) {
      tpar <- set.par(tpar, tau2, "tau2")
      suppressWarnings(opt <- try(optim(get.par(tpar, "b"), fn = objfun, gr = grad,
        method = "BFGS", control = list(), tau2 = tau2), silent = TRUE))
      if(!inherits(opt, "try-error")) {
        tpar <- set.par(tpar, opt$par, "b")
        x$state$fitted.values <- x$fit.fun(x$X, tpar)
      }
      x$state$parameters <- tpar
      edf <- x$edf(x)
      eta2[[id]] <- eta[[id]] + fitted(x$state)
      IC <- get.ic(family, y, family$map2par(eta2), edf0 + edf, length(eta2[[id]]), type = x$criterion, ...)
      IC
    }
    if(length(get.state(x, "tau2")) < 2) {
      tau2 <- try(optimize(objfun2, interval = c(0.1, 1000))$minimum, silent = TRUE)
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

  suppressWarnings(opt <- try(optim(get.par(tpar, "b"), fn = objfun, gr = grad,
    method = "BFGS", control = list(), tau2 = get.par(tpar, "tau2")), silent = TRUE))

  if(!inherits(opt, "try-error")) {
    tpar <- set.par(tpar, opt$par, "b")
    x$state$fitted.values <- x$fit.fun(x$X, tpar)
    x$state$parameters <- tpar
  }

  if(!x$fixed)
    x$state$edf <- x$edf(x)

  return(x$state)
}


set.all.par <- function(par, x, hessian = NULL)
{
  if(!is.null(hessian))
    hessian <- matrix_inv(-1 * hessian)
  nx <- names(x)
  np <- length(x)
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      tpar <- par[grep(paste("p", j, ".t", sj, ".", sep = ""), names(par), fixed = TRUE)]
      x[[nx[j]]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, get.par(tpar, "b"), "b")
      if(any(grepl("tau2", tpar))) {
        x[[nx[j]]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, get.par(tpar, "tau2"), "tau2")
      }
      x[[nx[j]]]$smooth.construct[[sj]]$state$fitted.values <- x[[nx[j]]]$smooth.construct[[sj]]$fit.fun(x[[nx[j]]]$smooth.construct[[sj]]$X,
        get.par(tpar, "b"))
      if(!is.null(hessian)) {
        ntpar <- names(tpar)[!grepl("tau2", names(tpar))]
        sigma <- hessian[ntpar, ntpar]
        if(any(eigen(sigma, symmetric = TRUE)$values <= 0)) {
          require("Matrix")
          sigma2 <- try(nearPD(sigma)$mat, silent = TRUE)
          if(inherits(sigma2, "try-error")) {
            sigma2 <- if(length(sigma) < 2) as.numeric(sigma) else diag(sigma)
            sigma2 <- if(length(sigma2) < 2) matrix(sigma2, 1, 1) else diag(sigma2)
          }
          sigma <- as.matrix(sigma2)
        }
        if(length(sigma) < 2) sigma <- matrix(sigma, 1, 1)
        x[[nx[j]]]$smooth.construct[[sj]]$state$hessian <- sigma
      }
    }
  }
  return(x)
}


log_posterior <- function(par, x, verbose = TRUE, show.edf = TRUE, digits = 3, scale = NULL)
{
  nx <- names(x)
  np <- length(nx)
  eta <- vector(mode = "list", length = np)
  names(eta) <- nx
  family <- attr(x, "family")
  lprior <- 0.0
  edf <- 0.0
  start <- 1
  for(j in 1:np) {
    eta[[nx[j]]] <- 0.0
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      g <- get.state(x[[nx[j]]]$smooth.construct[[sj]], "b")
      end <- length(g) + start - 1
      i <- start:end
      tpar <- par[i]
      names(tpar) <- names(g)
      start <- max(i) + 1
      x[[nx[j]]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, tpar, "b")
      x[[nx[j]]]$smooth.construct[[sj]]$state$fitted.values <- x[[nx[j]]]$smooth.construct[[sj]]$fit.fun(x[[nx[j]]]$smooth.construct[[sj]]$X,
        get.par(tpar, "b"))
      eta[[nx[j]]] <- eta[[nx[j]]] + fitted(x[[nx[j]]]$smooth.construct[[sj]]$state)
      lprior <- lprior + x[[nx[j]]]$smooth.construct[[sj]]$prior(c(tpar, get.state(x[[nx[j]]]$smooth.construct[[sj]], "tau2")))
      edf <- edf + x[[nx[j]]]$smooth.construct[[sj]]$edf(x[[nx[j]]]$smooth.construct[[sj]])
    }
  }
  ll <- family$loglik(attr(x, "y.vec"), family$map2par(eta))
  lp <- as.numeric(ll + lprior)

  if(verbose) {
    cat("\r")
    vtxt <- paste("logLik ", fmt(ll, width = 8, digits = digits),
      " logPost ", fmt(lp, width = 8, digits = digits),
      if(show.edf) paste(" edf ", fmt(edf, width = 6, digits = digits), sep = "") else NULL,
      " iteration ", formatC(bamlss_log_posterior_iteration, width = 4), sep = "")
    cat(vtxt)
    if(.Platform$OS.type != "unix") flush.console()
    bamlss_log_posterior_iteration <<- bamlss_log_posterior_iteration + 1
  }

  if(!is.null(scale))
    lp <- lp * scale

  return(lp)
}

grad_posterior <- function(par, x, ...)
{
  nx <- names(x)
  np <- length(nx)
  eta <- vector(mode = "list", length = np)
  names(eta) <- nx
  family <- attr(x, "family")
  grad <- NULL
  for(j in 1:np) {
    eta[[nx[j]]] <- 0
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      tpar <- par[grep(paste("p", j, ".t", sj, ".", sep = ""), names(par), fixed = TRUE)]
      x[[nx[j]]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, tpar, "b")
      x[[nx[j]]]$smooth.construct[[sj]]$state$fitted.values <- x[[nx[j]]]$smooth.construct[[sj]]$fit.fun(x[[nx[j]]]$smooth.construct[[sj]]$X,
        get.par(tpar, "b"))
      eta[[nx[j]]] <- eta[[nx[j]]] + fitted(x[[nx[j]]]$smooth.construct[[sj]]$state)
    }
  }
  for(j in 1:np) {
    if(!is.null(family$gradient[[nx[j]]])) {
      for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
        tp <- get.state(x[[nx[j]]]$smooth.construct[[sj]], "par")
        tgrad <- drop(family$gradient[[nx[j]]](get.par(tp, "b"), attr(x, "y.vec"), eta, x[[nx[j]]]$smooth.construct[[sj]]))
        if(!is.null(x[[nx[j]]]$smooth.construct[[sj]]$grad)) {
          tgrad <- tgrad + drop(x[[nx[j]]]$smooth.construct[[sj]]$grad(score = NULL, tp, full = FALSE))
        }
        grad <- c(grad, tgrad)
      }
    } else {
      score <- family$score[[nx[j]]](attr(x, "y.vec"), family$map2par(eta), id = nx[j])
      for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
        tgrad <- x[[nx[j]]]$smooth.construct[[sj]]$grad(score, x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, full = FALSE)
        grad <- c(grad, tgrad)
      }
    }
  }
  return(grad)
}


opt0 <- function(x, verbose = TRUE, digits = 3, gradient = TRUE, hessian = FALSE,
  eps = .Machine$double.eps^0.5, maxit = 100, ...)
{
  x <- bamlss.engine.setup(x, ...)

  par <- make_par(x, add.tau2 = FALSE)
  family <- attr(x, "family")

  if(!hessian) {
    if(verbose)
      bamlss_log_posterior_iteration <<- 1

    opt <- optim(par$par, fn = log_posterior,
      gr = if(!is.null(family$score) & gradient) grad_posterior else NULL,
      x = x, method = "BFGS", verbose = verbose, show.edf = FALSE,
      digits = digits, control = list(fnscale = -1, reltol = eps, maxit = maxit),
      hessian = TRUE)
 
    if(verbose) {
      cat("\n")
      rm(bamlss_log_posterior_iteration, envir = .GlobalEnv)
    }

    x <- set.all.par(opt$par, x, opt$hessian)
    attr(x, "hessian") <- opt$hessian
    attr(x, "converged") <- opt$convergence < 1

    return(x)
  } else {
    opt <- optimHess(par$par, fn = log_posterior,
      gr = if(!is.null(family$score) & gradient) grad_posterior else NULL,
      x = x, verbose = verbose, digits = digits,
      control = list(fnscale = -1, reltol = eps, maxit = maxit))
    return(opt)
  }
}


xbin.fun <- function(ind, weights, e, xweights, xrres, oind)
{
  .Call("xbin_fun", as.integer(ind), as.numeric(weights), 
    as.numeric(e), as.numeric(xweights), as.numeric(xrres),
    as.integer(oind))
  invisible(NULL)
}


## Censored normal backfitting.
bfit_cnorm <- function(x, criterion = c("AICc", "BIC", "AIC"),
  eps = .Machine$double.eps^0.25, maxit = 400, outer = FALSE, inner = FALSE,
  verbose = TRUE, digits = 4, ...)
{
  criterion <- match.arg(criterion)

  x <- bamlss.engine.setup(x, criterion = criterion, ...)

  family <- attr(x, "family")
  nx <- family$names[1:2]
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  criterion <- match.arg(criterion)

  np <- length(nx)
  y <- attr(x, "y.vec")
  nobs <- if(is.null(dim(y))) length(y) else nrow(y)
  eta <- get.eta(x)
  edf <- get.edf(x)

  ## Start the backfitting algorithm.
  eps0 <- eps + 1; iter <- 1
  while(eps0 > eps & iter < maxit) {
    eta0 <- eta
    ## Cycle through all parameters
    for(j in 1:np) {
      for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
        ## Get updated parameters.
        p.state <- x[[nx[j]]]$smooth.construct[[sj]]$update(x[[nx[j]]]$smooth.construct[[sj]],
          family, y, eta, nx[j], edf = edf)

        ## Compute equivalent degrees of freedom.
        edf <- edf - x[[nx[j]]]$smooth.construct[[sj]]$state$edf + p.state$edf

        ## Update predictor and smooth fit.
        eta[[nx[j]]] <- eta[[nx[j]]] - fitted(x[[nx[j]]]$smooth.construct[[sj]]$state) + fitted(p.state)

        x[[nx[j]]]$smooth.construct[[sj]]$state <- p.state
      }
    }

    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

    peta <- family$map2par(eta)

    if(verbose) {
      IC <- get.ic(family, y, peta, edf, nobs, criterion)
      cat("\r")
      vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
        " logLik ", fmt(family$loglik(y, peta), width = 8, digits = digits),
        " edf ", fmt(edf, width = 6, digits = digits),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
      cat(vtxt)

      if(.Platform$OS.type != "unix") flush.console()
    }

    iter <- iter + 1
  }

  if(verbose) cat("\n")

  return(x)
}


## Likelihood based boosting.
boost0 <- function(x, criterion = c("AICc", "BIC", "AIC"),
  nu = 1, df = 4, maxit = 100, mstop = NULL, best = TRUE,
  verbose = TRUE, digits = 4,
  eps = .Machine$double.eps^0.25, plot = TRUE, ...)
{
  if(!is.null(mstop))
    maxit <- mstop

  criterion <- match.arg(criterion)

  x <- bamlss.engine.setup(x, criterion = criterion, ...)

  family <- attr(x, "family")
  nx <- family$names[1:2]
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  criterion <- match.arg(criterion)

  np <- length(nx)
  y <- attr(x, "y.vec")
  nobs <- if(is.null(dim(y))) length(y) else nrow(y)

  ## Initialize select indicator and intercepts.
  states <- save.model.matrix <- list()
  eta <- get.eta(x)
  for(j in 1:np) {
    states[[j]] <- list()
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      x[[nx[j]]]$smooth.construct[[sj]] <- assign.df(x[[nx[j]]]$smooth.construct[[sj]], df)
      if(is.null(x[[nx[j]]]$smooth.construct[[sj]]$sp) & !x[[nx[j]]]$smooth.construct[[sj]]$fixed) {
        x[[nx[j]]]$smooth.construct[[sj]]$old.optimize <- x[[nx[j]]]$smooth.construct[[sj]]$state$do.optim
        x[[nx[j]]]$smooth.construct[[sj]]$state$do.optim <- FALSE
        x[[nx[j]]]$smooth.construct[[sj]]$do.optim <- FALSE
      }
      states[[j]][[sj]] <- list()
    }
    ii <- sapply(x[[nx[j]]]$smooth.construct, function(x) { inherits(x, "model.matrix") })
    if(any(ii)) {
      ii <- which(ii)[1]
      model.matrix <- list()
      cn <- colnames(x[[nx[j]]]$smooth.construct[[ii]]$X)
      g0 <- get.par(x[[nx[j]]]$smooth.construct[[ii]]$state$parameters, "b")
      for(pj in 1:ncol(x[[nx[j]]]$smooth.construct[[ii]]$X)) {
        model.matrix[[pj]] <- list()
        model.matrix[[pj]]$label <- cn[pj]
        model.matrix[[pj]]$term <- cn[pj]
        model.matrix[[pj]]$X <- x[[nx[j]]]$smooth.construct[[ii]]$X[, pj, drop = FALSE]
        model.matrix[[pj]]$xbin.ind <- x[[nx[j]]]$smooth.construct[[ii]]$xbin.ind
        model.matrix[[pj]]$xbin.sind <- x[[nx[j]]]$smooth.construct[[ii]]$xbin.sind
        model.matrix[[pj]]$xbin.k <- x[[nx[j]]]$smooth.construct[[ii]]$xbin.k
        model.matrix[[pj]]$xbin.order <- x[[nx[j]]]$smooth.construct[[ii]]$xbin.order
        model.matrix[[pj]]$nobs <- x[[nx[j]]]$smooth.construct[[ii]]$nobs
        model.matrix[[pj]]$fixed <- TRUE
        model.matrix[[pj]]$weights <- x[[nx[j]]]$smooth.construct[[ii]]$weights
        model.matrix[[pj]]$rres <- x[[nx[j]]]$smooth.construct[[ii]]$rres
        model.matrix[[pj]]$fit.fun <- x[[nx[j]]]$smooth.construct[[ii]]$fit.fun
        model.matrix[[pj]]$state <- list("parameters" = g0[pj])
        names(model.matrix[[pj]]$state$parameters) <- "g1"
        model.matrix[[pj]]$state$fitted.values <- drop(model.matrix[[pj]]$X %*% g0[pj])
        model.matrix[[pj]]$state$edf <- 0
        model.matrix[[pj]]$state$do.optim <- FALSE
        model.matrix[[pj]]$fixed <- TRUE
        model.matrix[[pj]]$is.model.matrix <- TRUE
        model.matrix[[pj]]$selected <- rep(0, length = maxit)
        model.matrix[[pj]]$upper <- Inf
        model.matrix[[pj]]$lower <- -Inf
        class(model.matrix[[pj]]) <- class(x[[nx[j]]]$smooth.construct[[ii]])
      }
      save.model.matrix[[j]] <- x[[nx[j]]]$smooth.construct[[ii]]
      x[[nx[j]]]$smooth.construct[[ii]] <- NULL
      x[[nx[j]]]$smooth.construct <- c(model.matrix, x[[nx[j]]]$smooth.construct)
    }
  }

  ## Selector help vector.
  sn <- NULL
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      sn <- c(sn, paste(j, sj, sep = "."))
      x[[nx[j]]]$smooth.construct[[sj]]$selected <- rep(0, length = maxit)
      x[[nx[j]]]$smooth.construct[[sj]]$loglik <- rep(0, length = maxit)
      x[[nx[j]]]$smooth.construct[[sj]]$state$edf <- 0
    }
  }
  select <- rep(0, length(sn))
  names(select) <- sn

  ## Initialize intercepts.
  eps0 <- eps + 1; iter <- 1
  while(eps0 > eps & iter < maxit) {
    eta0 <- eta
    ## Cycle through all parameters
    for(j in 1:np) {
      for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
        if(inherits(x[[nx[j]]]$smooth.construct[[sj]], "model.matrix") & x[[nx[j]]]$smooth.construct[[sj]]$label == "(Intercept)") {
          ## Get updated parameters.
          p.state <- bfit_iwls(x[[nx[j]]]$smooth.construct[[sj]],
            family, y, eta, nx[j], edf = edf)

          ## Update predictor and smooth fit.
          eta[[nx[j]]] <- eta[[nx[j]]] - fitted(x[[nx[j]]]$smooth.construct[[sj]]$state) + fitted(p.state)

          x[[nx[j]]]$smooth.construct[[sj]]$state <- p.state
        }
      }
    }

    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1
    iter <- iter + 1
  }

  ## Initial parameters.
  parameters <- get.all.par(x)

  ## Start boosting.
  eps0 <- 1; iter <- 1
  save.ic <- save.ll <- NULL
  ll <- family$loglik(y, family$map2par(eta))
  while(iter <= maxit) {
    eta0 <- eta

    ## Cycle through all parameters
    for(j in 1:np) {
      peta <- family$map2par(eta)

      ## Compute weights.
      weights <- family$hess[[nx[j]]](y, peta)

      ## Score.
      score <- family$score[[nx[j]]](y, peta)

      ## Compute working observations.
      z <- eta[[nx[j]]] + 1 / weights * score

      ## Residuals.
      resids <- z - eta[[nx[j]]]

      for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
        ## Get updated parameters.
        states[[j]][[sj]] <- boost0_iwls(x[[nx[j]]]$smooth.construct[[sj]],
          family, y, weights, resids, nu)

        ## Compute likelihood contribution.
        eta[[nx[j]]] <- eta[[nx[j]]] + fitted(states[[j]][[sj]])
        select[paste(j, sj, sep = ".")] <- -1 * (ll - family$loglik(y, family$map2par(eta)))
        eta[[nx[j]]] <- eta0[[nx[j]]]
      }
    }

    ## Which term to update.
    take <- as.integer(strsplit(sn[which.max(select)], ".", fixed = TRUE)[[1]])

    ## Update selected base learner.
    eta[[take[1]]] <- eta[[take[1]]] + fitted(states[[take[1]]][[take[2]]])

    ## Write to x.
    x[[take[1]]]$smooth.construct[[take[2]]]$state <- increase(x[[take[1]]]$smooth.construct[[take[2]]]$state, states[[take[1]]][[take[2]]])
    x[[take[1]]]$smooth.construct[[take[2]]]$selected[iter] <- 1
    x[[take[1]]]$smooth.construct[[take[2]]]$loglik[iter] <- max(select)

    edf <- get.edf(x, type = 2)

    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

    peta <- family$map2par(eta)
    IC <- get.ic(family, y, peta, edf, nobs, criterion)
    if(!is.null(save.ic) & best) {
      if(all(IC < save.ic))
        parameters <- get.all.par(x)
    }
    ll <- family$loglik(y, peta)

    save.ic <- c(save.ic, IC)
    save.ll <- c(save.ll, ll)

    if(verbose) {
      cat("\r")
      vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
        " logLik ", fmt(ll, width = 8, digits = digits),
        " edf ", fmt(edf, width = 6, digits = digits),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
      cat(vtxt)

      if(.Platform$OS.type != "unix") flush.console()
    }

    iter <- iter + 1
  }

  if(verbose) cat("\n")

  mstop <- which.min(save.ic)
  if(!best)
    parameters <- get.all.par(x)

  if(verbose) {
    cat("\n")
    cat(criterion, "=", save.ic[mstop], "-> at mstop =", mstop, "\n---\n")
  }
  labels <- NULL
  ll.contrib <- NULL
  bsum <- lmat <- vector(mode = "list", length = np)
  for(j in 1:np) {
    rn <- NULL
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      x[[nx[j]]]$smooth.construct[[sj]]$state$parameters <- parameters[[j]][[sj]]
      g <- get.par(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, "b")
      if(!is.null(tau2 <- attr(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, "true.tau2"))) {
        x[[nx[j]]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters,
          tau2, "tau2")
      }
      x[[nx[j]]]$smooth.construct[[sj]]$state$edf <- attr(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, "edf")
      if(is.null(x[[nx[j]]]$smooth.construct[[sj]]$state$edf))
        x[[nx[j]]]$smooth.construct[[sj]]$state$edf <- 0
      labels <- c(labels, paste(x[[nx[j]]]$smooth.construct[[sj]]$label, nx[j], sep = ":"))
      rn <- c(rn, x[[nx[j]]]$smooth.construct[[sj]]$label)
      bsum[[j]] <- rbind(bsum[[j]], sum(x[[nx[j]]]$smooth.construct[[sj]]$selected[1:mstop]) / mstop * 100)
      lmat[[j]] <- rbind(lmat[[j]], sum(x[[nx[j]]]$smooth.construct[[sj]]$loglik[1:mstop]))
      ll.contrib <- cbind(ll.contrib, cumsum(x[[nx[j]]]$smooth.construct[[sj]]$loglik))
    }
    if(!is.matrix(bsum[[j]])) bsum[[j]] <- matrix(bsum[[j]], nrow = 1)
    bsum[[j]] <- cbind(bsum[[j]], lmat[[j]])
    if(!is.matrix(bsum[[j]])) bsum[[j]] <- matrix(bsum[[j]], nrow = 1)
    colnames(bsum[[j]]) <- c(paste(nx[j], "% selected"), "LogLik contrib.")
    rownames(bsum[[j]]) <- rownames(lmat[[j]]) <- rn
    bsum[[j]] <- bsum[[j]][order(bsum[[j]][, 2], decreasing = TRUE), , drop = FALSE]
    if(verbose) {
      if(length(bsum[[j]]) < 2) print(round(bsum[[j]], digits = 4)) else printCoefmat(bsum[[j]], digits = 4)
      if(j != np)
        cat("---\n")
    }
  }
  if(verbose) cat("\n")

  colnames(ll.contrib) <- labels
  names(bsum) <- nx
  bsum <- list("summary" = bsum, "mstop" = mstop, "criterion" = criterion,
    "ic" = save.ic, "loglik" = ll.contrib)
  assign("boost.summary", bsum, envir = attr(x, "environment"))

  if(plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(1, 2), mar = c(5.1, 4.1, 2.1, 2.1))
    plot(save.ic, type = "l", xlab = "Iteration", ylab = criterion)
    abline(v = mstop, lwd = 3, col = "lightgray")
    axis(3, at = mstop, labels = paste("mstop =", mstop))
    par(mar = c(5.1, 4.1, 2.1, 10.1))
    matplot(1:nrow(ll.contrib), ll.contrib, type = "l", lty = 1,
      xlab = "Iteration", ylab = "LogLik contribution", col = "black")
    abline(v = mstop, lwd = 3, col = "lightgray")
    axis(4, at = ll.contrib[nrow(ll.contrib), ], labels = labels, las = 1)
    axis(3, at = mstop, labels = paste("mstop =", mstop))
  }

  ## Collect model.matrix effects.
  for(j in 1:np) {
    g <- drop <- NULL
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      if(!is.null(x[[nx[j]]]$smooth.construct[[sj]]$is.model.matrix)) {
        g <- c(g, x[[nx[j]]]$smooth.construct[[sj]]$state$parameters)
        drop <- c(drop, sj)
      }
    }
    if(!is.null(g)) {
      x[[nx[j]]]$smooth.construct[drop] <- NULL
      names(g) <- paste("b", 1:length(g), sep = "")
      save.model.matrix[[j]]$state$fitted.values <- drop(save.model.matrix[[j]]$X %*% g)
      save.model.matrix[[j]]$state$parameters <- g
      x[[nx[j]]]$smooth.construct <- c(list(save.model.matrix[[j]]), x[[nx[j]]]$smooth.construct)
    }
  }

  return(x)
}


## Boosting iwls.
boost0_iwls <- function(x, family, y, weights, resids, nu, ...)
{
  ## Initial parameters and fit.
  g0 <- get.par(x$state$parameters, "b")
  fit0 <- fitted(x$state)

  ## Compute reduced residuals.
  xbin.fun(x$xbin.sind, weights, resids, x$weights, x$rres, x$xbin.order)

  ## Compute mean and precision.
  XW <- x$X * x$weights
  XWX <- crossprod(x$X, XW)
  if(x$fixed) {
    P <- matrix_inv(XWX)
  } else {
    S <- 0
    tau2 <- get.state(x, "tau2")
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    P <- matrix_inv(XWX + S)
  }

  ## New parameters.
  g <- nu * drop(P %*% crossprod(x$X, x$rres))

  ## Finalize.
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))

  ## Find edf.
  xbin.fun(x$xbin.sind, weights, resids + fit0 + fitted(x$state), x$weights, x$rres, x$xbin.order)

  XW <- x$X * x$weights
  XWX <- crossprod(x$X, XW)
  if(x$fixed) {
    P <- matrix_inv(XWX)
  } else {
    g0 <- g0 + g

    objfun <- function(tau2) {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      P <- matrix_inv(XWX + S)
      g1 <- drop(P %*% crossprod(x$X, x$rres))
      sum((g1 - g0)^2)
    }

    if(length(get.state(x, "tau2")) < 2) {
      tau2 <- try(optimize(objfun, interval = x$state$interval)$minimum, silent = TRUE)
      if(inherits(tau2, "try-error"))
        tau2 <- optimize2(objfun, interval = x$state$interval, grid = x$state$grid)$minimum
    } else {
      i <- grep("tau2", names(x$lower))
      tau2 <- if(!is.null(x$state$true.tau2)) x$state$true.tau2 else get.state(x, "tau2")
      opt <- try(optim(tau2, fn = objfun, method = "L-BFGS-B",
        lower = x$lower[i], upper = x$upper[i]), silent = TRUE)
      if(!inherits(opt, "try-error"))
        tau2 <- opt$par
    }
    if(inherits(tau2, "try-error"))
      stop(paste("problem in finding optimum smoothing parameter for term ", x$label, "!", sep = ""))

    attr(x$state$parameters, "true.tau2") <- tau2

    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    P <- matrix_inv(XWX + S)
  }

  ## Assign degrees of freedom.
  x$state$edf <- sum(diag(XWX %*% P))
  attr(x$state$parameters, "edf") <- x$state$edf

  return(x$state)
}


## Increase coeffiecients.
increase <- function(state0, state1)
{
  g <- get.par(state0$parameters, "b") + get.par(state1$parameters, "b")
  state0$fitted.values <- fitted(state0) + fitted(state1)
  state0$parameters <- set.par(state0$parameters, g, "b")
  state0$edf <- state1$edf
  attr(state0$parameters, "true.tau2") <- attr(state1$parameters, "true.tau2")
  attr(state0$parameters, "edf") <- attr(state1$parameters, "edf")
  state0
}


## Smallish summary function.
sboost <- function(object, plot = TRUE, ...)
{
  bs <- get("boost.summary", envir = attr(object, "environment"))
  np <- length(bs$summary)

  cat("\n")
  cat(bs$criterion, "=", bs$ic[bs$mstop], "-> at mstop =", bs$mstop, "\n---\n")
  for(j in 1:np) {
    if(length(bs$summary[[j]]) < 2) {
      print(round(bs$summary[[j]], digits = 4))
    } else printCoefmat(bs$summary[[j]], digits = 4)
    if(j != np)
      cat("---\n")
  }
  cat("\n")

  if(plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(1, 2), mar = c(5.1, 4.1, 2.1, 2.1))
    plot(bs$ic, type = "l", xlab = "Iteration", ylab = bs$criterion)
    abline(v = bs$mstop, lwd = 3, col = "lightgray")
    axis(3, at = bs$mstop, labels = paste("mstop =", bs$mstop))
    par(mar = c(5.1, 4.1, 2.1, 10.1))
    matplot(1:nrow(bs$loglik), bs$loglik, type = "l", lty = 1,
      xlab = "Iteration", ylab = "LogLik contribution", col = "black")
    abline(v = bs$mstop, lwd = 3, col = "lightgray")
    axis(4, at = bs$loglik[nrow(bs$loglik), ], labels = colnames(bs$loglik), las = 1)
    axis(3, at = bs$mstop, labels = paste("mstop =", bs$mstop))
  }

  return(invisible(bs))
}


## Assign starting values.
add.starting.values <- function(x, start)
{
  if(!is.null(start)) {
    nx <- names(x)
    for(id in nx) {
      if(!is.null(x[[id]]$smooth.construct)) {
        if(!is.null(x[[id]]$smooth.construct$model.matrix)) {
          if(length(take <- grep(paste(id, "p", sep = "."), names(start), fixed = TRUE, value = TRUE))) {
            cn <- paste(id, "p", colnames(x[[id]]$smooth.construct$model.matrix$X), sep = ".")
            i <- grep2(take, cn, fixed = TRUE)
            if(length(i)) {
              tpar <- start[take[i]]
              names(tpar) <- gsub(paste(id, "p.", sep = "."), "", names(tpar), fixed = TRUE)
              i <- grep2(c("edf", "accepted", "alpha"), names(tpar))
              x[[id]]$smooth.construct$model.matrix$state$parameters <- if(length(i)) tpar[-i] else tpar
              x[[id]]$smooth.construct$model.matrix$state$fitted.values <- x[[id]]$smooth.construct$model.matrix$fit.fun(x[[id]]$smooth.construct$model.matrix$X, x[[id]]$smooth.construct$model.matrix$state$parameters)
            }
          }
        }
        for(j in seq_along(x[[id]]$smooth.construct)) {
          tl <- x[[id]]$smooth.construct[[j]]$label
          if(x[[id]]$smooth.construct[[j]]$by != "NA")
            tl <- paste(tl, x[[id]]$smooth.construct[[j]]$by, sep = ":")
          take <- grep(tl <- paste(id, "s", tl, sep = "."),
            names(start), fixed = TRUE, value = TRUE)
          if(x[[id]]$smooth.construct[[j]]$by == "NA") {
            take <- take[!grepl(paste(tl, ":", sep = ""), take, fixed = TRUE)]
          }
          if(length(take)) {
            tpar <- start[take]
            names(tpar) <- gsub(paste(tl, ".", sep = ""), "", names(tpar), fixed = TRUE)
            i <- grep2(c("edf", "accepted", "alpha"), names(tpar))
            x[[id]]$smooth.construct[[j]]$state$parameters <- if(length(i)) tpar[-i] else tpar
            x[[id]]$smooth.construct[[j]]$state$fitted.values <- x[[id]]$smooth.construct[[j]]$fit.fun(x[[id]]$smooth.construct[[j]]$X, x[[id]]$smooth.construct[[j]]$state$parameters)
          }
        }
      }
    }
  }
  return(x)
}

