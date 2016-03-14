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
  do.optim = NULL, optim.grid = NULL, criterion = c("AICc", "BIC", "AIC"),
  nu = 0.1, start = NULL, df = NULL, ...)
{
  if(!is.null(attr(x, "bamlss.engine.setup"))) return(x)

  criterion <- match.arg(criterion)
  if(!is.null(optim.grid)) {
    if(is.na(optim.grid))
      optim.grid <- NULL
    if(!optim.grid)
      optim.grid <- NULL
  }

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
          x$smooth.construct[[j]] <- bamlss.engine.setup.smooth(x$smooth.construct[[j]], ...)
          tdf <- NULL
          if(!is.null(df)) {
            if(!is.null(names(df))) {
              if((x$smooth.construct[[j]]$label %in% names(df)))
                tdf <- df[x$smooth.construct[[j]]$label]
            } else tdf <- df[1]
          }
          x$smooth.construct[[j]] <- assign.df(x$smooth.construct[[j]], tdf)
          x$smooth.construct[[j]]$optim.grid <- optim.grid
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
        return(p[!grepl("tau", names(p)) & !grepl("edf", names(p))])
      } else return(x$state[[what]])
    }
  }
}

get.par <- function(x, what = NULL) {
  if(is.null(what) | is.null(names(x))) return(x)
  if(what %in% c("tau2", "tau", "lambda")) {
    return(x[grep("tau", names(x))])
  } else {
    if(what %in% "b") {
      return(x[!grepl("tau", names(x)) & !grepl("edf", names(x))])
    } else return(x[what])
  }
}

set.par <- function(x, replacement, what) {
  if(is.null(replacement))
    return(x)
  if(what %in% c("tau2", "tau", "lambda")) {
    x[grep("tau", names(x))] <- replacement
  } else {
    if(what %in% "b") {
      x[!grepl("tau", names(x)) & !grepl("edf", names(x))] <- replacement
    } else x[what] <- replacement
  }
  x
}

## The default method.
bamlss.engine.setup.smooth.default <- function(x, spam = FALSE, ...)
{
  if(inherits(x, "special"))
    return(x)
  if(is.null(x$binning) & !is.null(x$xt[["binning"]])) {
    x$binning <- match.index(x$X)
    x$binning$order <- order(x$binning$match.index)
    x$binning$sorted.index <- x$binning$match.index[x$binning$order]
    x$X <- x$X[x$binning$nodups, , drop = FALSE]
  }
  if(!is.null(x$binning)) {
    if(nrow(x$X) != length(x$binning$nodups))
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
  x$fit.reduced <- rep(0, length = k)
  state <- if(is.null(x$xt[["state"]])) list() else x$xt[["state"]]
  if(is.null(x$fixed))
    x$fixed <- if(!is.null(x$fx)) x$fx[1] else FALSE
  if(!x$fixed & is.null(state$interval))
    state$interval <- if(is.null(x$xt[["interval"]])) tau2interval(x) else x$xt[["interval"]]
  state$grid <- if(is.null(x$xt[["grid"]])) 10 else x$xt[["grid"]]
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
            rep(if(!is.null(x$xt[["tau2"]])) {
              x$xt[["tau2"]]
            } else {
              if(!is.null(x$xt[["lambda"]])) 1 / x$xt[["lambda"]] else 1000
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
        rep(if(!is.null(x$xt[["tau2"]])) {
          x$xt[["tau2"]]
        } else {
          if(!is.null(x$xt[["lambda"]])) 1 / x$xt[["lambda"]] else 100
        }, length.out = ntau2)
      }
    } else rep(x$sp, length.out = ntau2)
    names(tau2) <- paste("tau2", 1:ntau2, sep = "")
    state$parameters <- c(state$parameters, tau2)
  }
  x$a <- if(is.null(x$xt[["a"]])) 1e-04 else x$xt[["a"]]
  x$b <- if(is.null(x$xt[["b"]])) 1e-04 else x$xt[["b"]]
  if(!is.null(x$xt[["prior"]]))
    x$prior <- x$xt[["prior"]]
  if(is.null(x$prior) | !is.function(x$prior))
    x$prior <- make.prior(x)
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
      S <- 0
      for(j in seq_along(tau2))
        S <- S + 1 / tau2[j] * x$S[[j]]
      P <- matrix_inv(x$state$XX + S, index = x$sparse.setup)
      edf <- sum.diag(x$state$XX %*% P)
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
          tauS <- -1 / tau2[j] * x$S[[j]]
          grad <- grad + tauS %*% gamma
          if(full & !is.null(tau2[j])) {
            grad2 <- c(grad2, drop(-x$rank[j] / (2 * tau2[j]) - 1 / (2 * tau2[j]^2) * x$S[[j]] %*% gamma + (-x$a - 1) / tau2[j] + x$b / (tau2[j]^2)))
            x$X <- cbind(x$X, 0)
          }
          grad <- drop(grad)
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
  if(!is.null(x$xt[["do.optim"]]))
    x$state$do.optim <- x$xt[["do.optim"]]
  x$state$edf <- x$edf(x)
  x$sparse.setup <- sparse.setup(x$X, S = x$S)
  x$fit.fun <- make.fit.fun(x)
  x$state$fitted.values <- x$fit.fun(x$X, get.par(x$state$parameters, "b"))
  x$added <- c("nobs", "weights", "rres", "state", "grid", "a", "b", "prior", "edf",
    "grad", "hess", "lower", "upper")

  if(spam) {
    require("spam")
    x$X <- as.spam(x$X)
    xx <- crossprod.spam(x$X)
    for(j in seq_along(x$S)) {
      x$S[[j]] <- as.spam(x$S[[j]])
      xx <- xx + x$S[[j]]
    }
    x$sparse.setup$spam.cholFactor <- chol.spam(xx)
  }

  if(ntau2 > 0) {
    tau2 <- NULL
    if(length(x$margin)) {
      for(j in seq_along(x$margin)) {
        if(!is.null(x$margin[[j]]$xt$tau2))
          tau2 <- c(tau2, x$margin[[j]]$xt$tau2)
      }
    } else {
      if(!is.null(x$xt$tau2))
        tau2 <- x$xt$tau2
    }
    if(!is.null(tau2)) {
      tau2 <- rep(tau2, length.out = ntau2)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    }
  }

  x
}


## Function to find tau2 interval according to the
## effective degrees of freedom
tau2interval <- function(x, lower = .Machine$double.eps^0.8, upper = 1e+10)
{
  if(length(x$S) < 2) {
    return(c(lower, upper))
  } else {
    return(rep(list(c(lower, upper)), length.out = length(x$S)))
  }
}


## Assign degrees of freedom.
assign.df <- function(x, df)
{
  tau2 <- get.par(x$state$parameters, "tau2")
  if(x$fixed | !length(tau2))
    return(x)
  df <- if(is.null(x$xt$df)) df else x$xt$df
  if(is.null(df)) {
    nc <- ncol(x$X)
    df <- ceiling(nc * if(nc < 11) 0.9 else 0.5)
  }
  if(df > ncol(x$X))
    df <- ncol(x$X)
  if(df < 1)
    df <- 1
  int <- c(.Machine$double.eps^0.25, 1e+10)
  XX <- if(inherits(x$X, "spam")) crossprod.spam(x$X) else crossprod(x$X)
  if(length(tau2) > 1) {
    if(FALSE) {
      df.part <- df / length(tau2)
      for(j in seq_along(tau2)) {
        objfun <- function(val) {
          tau2[j] <- val
          S <- 0
          for(i in seq_along(x$S))
            S <- S + 1 / tau2[i] * x$S[[i]]
          edf <- sum.diag(XX %*% matrix_inv(XX + S, index = x$sparse.setup))
          return((df - edf)^2)
        }
        opt <- try(optimize(objfun, int)$minimum, silent = TRUE)
        if(!inherits(opt, "try-error"))
          tau2[j] <- opt
      }
      tau2 <- rep(1000, length(tau2))
    }
  } else {
    objfun <- function(tau2) {
      edf <- sum.diag(XX %*% matrix_inv(XX + 1 / tau2 * x$S[[1]], index = x$sparse.setup))
      return((df - edf)^2)
    }
    tau2 <- try(optimize(objfun, int)$minimum, silent = TRUE)
    if(inherits(tau2, "try-error"))
      return(x)
  }
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

## Initialze.
init.eta <- function(eta, y, family, nobs)
{
  if(is.null(family$initialize))
    return(eta)
  for(j in family$names) {
    if(!is.null(family$initialize[[j]])) {
      linkfun <- make.link2(family$links[j])$linkfun
      eta[[j]] <- linkfun(family$initialize[[j]](y))
    }
  }
  return(eta)
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

get.all.par <- function(x, drop = FALSE, list = TRUE)
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
              if(!is.null(edf <- x[[i]][[k]]$smooth.construct[[j]]$state$edf))
                par[[i]][[k]]$s[[j]] <- c(par[[i]][[k]]$s[[j]], "edf" = edf)
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
            if(!is.null(edf <- x[[i]]$smooth.construct[[j]]$state$edf))
              par[[i]]$s[[j]] <- c(par[[i]]$s[[j]], "edf" = edf)
          }
        }
      }
    }
  }
  if(!list) {
    par <- unlist(par)
    if(drop) {
      for(j in c(".edf", ".tau2", ".alpha"))
        par <- par[!grepl(j, names(par), fixed = TRUE)]
    }
  }
  par
}


get.hessian <- function(x)
{
  npar <- names(get.all.par(x, list = FALSE, drop = TRUE))
  hessian <- list(); nh <- NULL
  for(i in names(x)) {
    for(j in names(x[[i]]$smooth.construct)) {
      pn <- if(j == "model.matrix") paste(i, "p", sep = ".") else paste(i, "s", j, sep = ".")
      hessian[[pn]] <- x[[i]]$smooth.construct[[j]]$state$hessian
      cn <- colnames(x[[i]]$smooth.construct[[j]]$X)
      if(is.null(cn))
        cn <- paste("b", 1:ncol(x[[i]]$smooth.construct[[j]]$X), sep = "")
      pn <- paste(pn, cn, sep = ".")
      nh <- c(nh, pn)
    }
  }
  require("Matrix")
  hessian <- -1 * as.matrix(do.call("bdiag", hessian))
  rownames(hessian) <- colnames(hessian) <- nh
  hessian <- hessian[npar, npar]
  return(hessian)
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
  maxit = 400, outer = FALSE, inner = FALSE, mgcv = FALSE,
  verbose = TRUE, digits = 4, ...)
{
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("design construct names mismatch with family names!")

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, ...)

  criterion <- match.arg(criterion)
  np <- length(nx)
  nobs <- nrow(y)
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  if(!is.null(start))
    x <- set.starting.values(x, start)

  eta <- init.eta(get.eta(x), y, family, nobs)
  ##eta$bd <- rep(1, nobs)
  
  if(!is.null(weights))
    weights <- as.data.frame(weights)
  if(!is.null(offset)) {
    offset <- as.data.frame(offset)
    for(j in nx) {
      if(!is.null(offset[[j]]))
        eta[[j]] <- eta[[j]] + offset[[j]]
    }
  }

  ia <- interactive()

  if(mgcv) {
    outer <- TRUE
    inner <- TRUE
  }

  if(!mgcv) {
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
  } else {
    inner_bf <- function(x, y, eta, family, edf, id, z, hess, weights, ...) {
      X <- lapply(x, function(x) { x$X })
      S <- lapply(x, function(x) { x$S })
      nt <- nt0 <- names(X)
      nt <- rmf(nt)
      names(X) <- names(S) <- nt
      if("modelmatrix" %in% nt)
        S <- S[!(nt %in% "modelmatrix")]
      X$z <- z
      f <- paste("z", paste(c("-1", nt), collapse = " + "), sep = " ~ ")
      f <- as.formula(f)
      if(!is.null(weights))
        hess <- hess * weights
      b <- gam(f, data = X, weights = hess, paraPen = S)
      cb <- coef(b)
      ncb <- names(cb)
      tau2 <- if(length(b$sp)) 1 / b$sp else NULL
      fitted <- 0
      for(sj in seq_along(x)) {
        tn <- rmf(nt0[sj])
        par <- cb[grep(tn, ncb, fixed = TRUE)]
        tedf <- sum(b$edf[grep(tn, ncb, fixed = TRUE)])
        names(par) <- paste("b", 1:length(par), sep = "")
        if(!is.null(tau2) & (tn != "modelmatrix")) {
          ttau2 <- tau2[grep(tn, names(tau2), fixed = TRUE)]
          names(ttau2) <- paste("tau2", 1:length(ttau2), sep = "")
          lo <- x[[sj]]$lower[grep("tau2", names(x[[sj]]$lower), fixed = TRUE)]
          up <- x[[sj]]$upper[grep("tau2", names(x[[sj]]$upper), fixed = TRUE)]
          if(any(j <- ttau2 < lo))
            ttau2[j] <- lo[j]
          if(any(j <- ttau2 > up))
            ttau2[j] <- up[j]
          par <- c(par, ttau2)
        } else {
          names(par) <- colnames(x[[sj]]$X)
          par <- c(par, "tau21" = 1e+20)
        }
        x[[sj]]$state$parameters <- par
        x[[sj]]$state$fitted.values <- x[[sj]]$fit.fun(x[[sj]]$X, par)
        fitted <- fitted + x[[sj]]$state$fitted.values
        edf <- edf - x[[sj]]$state$edf + tedf
        x[[sj]]$state$edf <- tedf
        x[[sj]]$state$prior <- x[[sj]]$prior(par)
      }
      eta[[id]] <- fitted
      return(list("x" = x, "eta" = eta, "edf" = edf))
    }
  }

  ## Backfitting main function.
  backfit <- function(verbose = TRUE) {
    eps0 <- eps + 1; iter <- 1
    edf <- get.edf(x, type = 2)
    while(eps0 > eps & iter < maxit) {
      eta0 <- eta
      ## Cycle through all parameters
      for(j in 1:np) {
        if(outer | iter < 2) {
          peta <- family$map2par(eta)

          ## Compute weights.
          hess <- family$hess[[nx[j]]](y, peta, id = nx[j])

          ## Score.
          score <- family$score[[nx[j]]](y, peta, id = nx[j])

          ## Compute working observations.
          z <- eta[[nx[j]]] + 1 / hess * score
        } else z <- hess <- NULL

        if(iter < 2)
          eta[[nx[j]]] <- get.eta(x)[[nx[j]]]

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
              family, y, eta, nx[j], edf = edf, z = z, hess = hess, weights = weights[[nx[j]]],
              iteration = iter)

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
        cat(if(ia) "\r" else "\n")
        vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
          " logPost ", fmt(family$loglik(y, peta) + get.log.prior(x), width = 8, digits = digits),
          " logLik ", fmt(family$loglik(y, peta), width = 8, digits = digits),
          " edf ", fmt(edf, width = 6, digits = digits),
          " eps ", fmt(eps0, width = 6, digits = digits + 2),
          " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
        cat(vtxt)

        if(.Platform$OS.type != "unix" & ia) flush.console()
      }

      iter <- iter + 1
    }

    IC <- get.ic(family, y, peta, edf, nobs, criterion)
    logLik <- family$loglik(y, peta)
    logPost <- as.numeric(logLik + get.log.prior(x))

    if(verbose) {
      cat(if(ia) "\r" else "\n")
      vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
        " logPost ", fmt(logPost, width = 8, digits = digits),
        " logLik ", fmt(family$loglik(y, peta), width = 8, digits = digits),
        " edf ", fmt(edf, width = 6, digits = digits),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
      cat(vtxt)
      if(.Platform$OS.type != "unix" & ia) flush.console()
      cat("\n")
    }

    if(iter == maxit)
      warning("the backfitting algorithm did not converge!")

    names(IC) <- criterion

    rval <- list("fitted.values" = eta, "parameters" = get.all.par(x), "edf" = edf,
      "logLik" = logLik, "logPost" = logPost, "converged" = iter < maxit)
    rval[[names(IC)]] <- IC
    rval
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

get.ic2 <- function(logLik, edf, n, type = c("AIC", "BIC", "AICc", "MP"), ...)
{
  type <- match.arg(type)
  pen <- switch(type,
    "AIC" = -2 * logLik + 2 * edf,
    "BIC" = -2 * logLik + edf * log(n),
    "AICc" = -2 * logLik + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
    "MP" = -1 * (logLik + edf)
  )
  return(pen)
}


cround <- function(x, digits = 2)
{
  cdigits <- Vectorize(function(x) {
    if(abs(x) >= 1)
      return(0)
    scipen <- getOption("scipen")
    on.exit(options("scipen" = scipen))
    options("scipen" = 100)
    x <- strsplit(paste(x), "")[[1]]
    x <- x[which(x == "."):length(x)][-1]
    i <- which(x != "0")
    x <- x[1:(i[1] - 1)]
    n <- length(x)
    if(n < 2) {
      if(x != "0")
        return(1)
      else return(n + 1)
    } else return(n + 1)
  })

  round(x, digits = cdigits(x) + digits)
}


## Naive smoothing parameter optimization.
tau2.optim <- function(f, start, ..., scale = 10, eps = 0.0001, maxit = 10)
{
  foo <- function(par, start, k) {
    start[k] <- par
    return(f(start, ...))
  }

  start <- cround(start)

  iter <- 1; eps0 <- eps + 1
  while((eps0 > eps) & (iter < maxit)) {
    start0 <- start
    for(k in seq_along(start)) {
      xr <- c(start[k] / scale, start[k] * scale)
      tpar <- try(optimize(foo, interval = xr, start = start, k = k), silent = TRUE)
      if(!inherits(tpar, "try-error"))
        start[k] <- cround(tpar$minimum)
    }
    if(length(start) < 2)
      break

    eps0 <- mean(abs((start - start0) / start0))
    iter <- iter + 1
  }

  return(start)
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

  Sigma <- matrix_inv(g.hess, index = x$sparse.setup)

  g <- drop(g + nu * Sigma %*% g.grad)

  x$state$parameters <- set.par(x$state$parameters, g, "b")
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
  x$state$hessian <- Sigma

  return(x$state)
}


bfit_lm <- function(x, family, y, eta, id, weights, ...)
{
  args <- list(...)

  peta <- family$map2par(eta)

  hess <- family$hess[[id]](y, peta, id = id, ...)

  ## Score.
  score <- family$score[[id]](y, peta, id = id, ...)

  ## Compute working observations.
  z <- eta[[id]] + 1 / hess * score

  ## Compute reduced residuals.
  e <- z - eta[[id]] + fitted(x$state)

  if(!is.null(weights))
    hess <- hess * weights
  if(x$fixed | x$fxsp) {
    b <- lm.wfit(x$X, e, hess)
  } else {
    tau2 <- get.par(x$state$parameters, "tau2")
    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    n <- nrow(S)
    w <- c(hess, rep(0, n))
    e <- c(e, rep(1, n))
    b <- lm.wfit(rbind(x$X, S), e, w)
  }

  x$state$parameters <- set.par(x$state$parameters, coef(b), "b")
  x$state$fitted.values <- x$X %*% coef(b)

  x$state
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
  XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)
  if(!x$state$do.optim | x$fixed | x$fxsp) {
    if(x$fixed) {
      P <- matrix_inv(XWX, index = x$sparse.setup)
    } else {
      S <- 0
      tau2 <- get.state(x, "tau2")
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      P <- matrix_inv(XWX + S, index = x$sparse.setup)
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
      P <- matrix_inv(XWX + S, index = x$sparse.setup)
      if(inherits(P, "try-error")) return(NA)
      g <- drop(P %*% crossprod(x$X, x$rres))
      if(any(is.na(g)) | any(g %in% c(-Inf, Inf))) g <- rep(0, length(g))
      fit <- x$fit.fun(x$X, g)
      edf <- sum.diag(XWX %*% P)
      eta2[[id]] <- eta2[[id]] + fit
      IC <- get.ic(family, y, family$map2par(eta2), edf0 + edf, length(z), x$criterion, ...)
      return(IC)
    }

    tau2 <- tau2.optim(objfun, start = get.state(x, "tau2"))
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    P <- matrix_inv(XWX + S, index = x$sparse.setup)
    x$state$parameters <- set.par(x$state$parameters, drop(P %*% crossprod(x$X, x$rres)), "b")
  }

  ## Compute fitted values.
  g <- get.state(x, "b")
  if(any(is.na(g)) | any(g %in% c(-Inf, Inf)))
    x$state$parameters <- set.par(x$state$parameters, rep(0, length(x$state$g)), "b")
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
  x$state$edf <- sum.diag(XWX %*% P)
  if(!is.null(x$prior)) {
    if(is.function(x$prior))
      x$state$log.prior <- x$prior(x$state$parameters)
  }

  return(x$state)
}


bfit_iwls_spam <- function(x, family, y, eta, id, weights, ...)
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
  XWX <- (t(diag.spam(x$weights) %*% x$X)) %*% x$X
  Xr <- crossprod.spam(x$X, x$rres)

  if(!x$state$do.optim | x$fixed | x$fxsp) {
    if(!x$fixed) {
      tau2 <- get.state(x, "tau2")
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      U <- update.spam.chol.NgPeyton(x$sparse.setup$spam.cholFactor, XWX + S)
    } else {
      U <- update.spam.chol.NgPeyton(x$sparse.setup$spam.cholFactor, XWX)
    }
    P <- chol2inv.spam(U)
    b <- P %*% Xr
    x$state$parameters <- set.par(x$state$parameters, b, "b")
  } else {
    args <- list(...)
    edf0 <- args$edf - x$state$edf
    eta2 <- eta
    objfun <- function(tau2, ...) {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      U <- update.spam.chol.NgPeyton(x$sparse.setup$spam.cholFactor, XWX + S)
      P <- chol2inv.spam(U)
      b <- P %*% Xr
      fit <- x$fit.fun(x$X, b)
      edf <- sum.diag(XWX %*% P)
      eta2[[id]] <- eta2[[id]] + fit
      IC <- get.ic(family, y, family$map2par(eta2), edf0 + edf, length(z), x$criterion, ...)
      return(IC)
    }
    if(length(get.state(x, "tau2")) < 2) {
      if(is.null(x$optim.grid)) {
        tau2 <- try(optimize(objfun, interval = x$state$interval)$minimum, silent = TRUE)
        if(inherits(tau2, "try-error"))
          tau2 <- optimize2(objfun, interval = x$state$interval, grid = x$state$grid)$minimum
      } else {
        tau2 <- optimize2(objfun, interval = x$state$interval, grid = x$state$grid)$minimum
      }
      x$state$parameters <- set.par(x$state$parameters, if(!length(tau2)) x$interval[1] else tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$lower))
      opt <- try(optim(get.state(x, "tau2"), fn = objfun, method = "L-BFGS-B",
        lower = x$lower[i], upper = x$upper[i]), silent = TRUE)
      if(!inherits(opt, "try-error"))
        x$state$parameters <- set.par(x$state$parameters, opt$par, "tau2")
    }
    tau2 <- get.state(x, "tau2")
    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    U <- update.spam.chol.NgPeyton(x$sparse.setup$spam.cholFactor, XWX + S)
    P <- chol2inv.spam(U)
    b <- P %*% Xr
    x$state$parameters <- set.par(x$state$parameters, b, "b")
  }

  ## Compute fitted values.
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
  x$state$edf <- sum.diag(XWX %*% P)
  if(!is.null(x$prior)) {
    if(is.function(x$prior))
      x$state$log.prior <- x$prior(x$state$parameters)
  }

  return(x$state)
}


## Updating based on optim.
bfit_optim <- function(x, family, y, eta, id, weights, ...)
{
  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - fitted(x$state)
  eta2 <- eta

  tpar <- x$state$parameters

  ## Objective for regression coefficients.
  objfun <- function(gamma, tau2 = NULL) {
    tpar <- set.par(tpar, gamma, "b")
    if(!is.null(tau2) & !x$fixed)
      tpar <- set.par(tpar, tau2, "tau2")
    eta2[[id]] <- eta[[id]] + x$fit.fun(x$X, tpar)
    ll <- if(is.null(weights)) {
      family$loglik(y, family$map2par(eta2))
    } else {
      sum(family$d(y, family$map2par(eta2)) * weights, na.rm = TRUE)
    }
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


## Compute fitted.values from set of parameters.
get.eta.par <- function(par, x)
{
  nx <- names(x)
  eta <- vector(mode = "list", length = length(nx))
  names(eta) <- nx
  for(j in nx) {
    eta[[j]] <- 0.0
    for(sj in names(x[[j]]$smooth.construct)) {
      xl <- if(sj != "model.matrix") {
        paste(j, "s", x[[j]]$smooth.construct[[sj]]$label, sep = ".")
      } else {
        paste(j, "p", strsplit(x[[j]]$smooth.construct[[sj]]$label, "+", fixed = TRUE)[[1]], sep = ".")
      }
      tpar <- par[grep2(xl, names(par), fixed = TRUE)]
      x[[j]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[j]]$smooth.construct[[sj]]$state$parameters, tpar, "b")
      x[[j]]$smooth.construct[[sj]]$state$fitted.values <- x[[j]]$smooth.construct[[sj]]$fit.fun(x[[j]]$smooth.construct[[sj]]$X,
        get.par(tpar, "b"))
      eta[[j]] <- eta[[j]] + fitted(x[[j]]$smooth.construct[[sj]]$state)
    }
    if(!is.null(x[[j]]$model.matrix)) {
      xl <- paste(j, "p", colnames(x[[j]]$model.matrix), sep = ".")
      tpar <- par[grep(xl, names(par), fixed = TRUE)]
      eta[[j]] <- eta[[j]] + drop(x[[j]]$model.matrix %*% tpar)
    }
  }
  return(eta)
}


## The log-posterior.
log_posterior <- function(par, x, y, family, verbose = TRUE, digits = 3, scale = NULL)
{
  nx <- names(x)
  eta <- vector(mode = "list", length = length(nx))
  names(eta) <- nx
  lprior <- 0.0
  for(j in nx) {
    eta[[j]] <- 0.0
    for(sj in names(x[[j]]$smooth.construct)) {
      xl <- if(sj != "model.matrix") {
        paste(j, "s", x[[j]]$smooth.construct[[sj]]$label, sep = ".")
      } else {
        paste(j, "p", strsplit(x[[j]]$smooth.construct[[sj]]$label, "+", fixed = TRUE)[[1]], sep = ".")
      }
      tpar <- par[grep2(xl, names(par), fixed = TRUE)]
      x[[j]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[j]]$smooth.construct[[sj]]$state$parameters, tpar, "b")
      x[[j]]$smooth.construct[[sj]]$state$fitted.values <- x[[j]]$smooth.construct[[sj]]$fit.fun(x[[j]]$smooth.construct[[sj]]$X,
        get.par(tpar, "b"))
      eta[[j]] <- eta[[j]] + fitted(x[[j]]$smooth.construct[[sj]]$state)
      lprior <- lprior + x[[j]]$smooth.construct[[sj]]$prior(c(tpar, get.state(x[[j]]$smooth.construct[[sj]], "tau2")))
    }
  }
  ll <- family$loglik(y, family$map2par(eta))
  lp <- as.numeric(ll + lprior)

  if(verbose) {
    cat(if(interactive()) "\r" else "\n")
    vtxt <- paste("logLik ", fmt(ll, width = 8, digits = digits),
      " logPost ", fmt(lp, width = 8, digits = digits),
      " iteration ", formatC(bamlss_log_posterior_iteration, width = 4), sep = "")
    cat(vtxt)
    if(.Platform$OS.type != "unix" & interactive()) flush.console()
    bamlss_log_posterior_iteration <<- bamlss_log_posterior_iteration + 1
  }

  if(!is.null(scale))
    lp <- lp * scale

  return(lp)
}


## Gradient vecor of the log-posterior.
grad_posterior <- function(par, x, y, family, ...)
{
  nx <- names(x)
  eta <- vector(mode = "list", length = length(nx))
  names(eta) <- nx
  grad <- NULL
  for(j in nx) {
    eta[[j]] <- 0
    for(sj in names(x[[j]]$smooth.construct)) {
      xl <- if(sj != "model.matrix") {
        paste(j, "s", x[[j]]$smooth.construct[[sj]]$label, sep = ".")
      } else {
        paste(j, "p", strsplit(x[[j]]$smooth.construct[[sj]]$label, "+", fixed = TRUE)[[1]], sep = ".")
      }
      tpar <- par[grep2(xl, names(par), fixed = TRUE)]
      x[[j]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[j]]$smooth.construct[[sj]]$state$parameters, tpar, "b")
      x[[j]]$smooth.construct[[sj]]$state$fitted.values <- x[[j]]$smooth.construct[[sj]]$fit.fun(x[[j]]$smooth.construct[[sj]]$X,
        get.par(tpar, "b"))
      eta[[j]] <- eta[[j]] + fitted(x[[j]]$smooth.construct[[sj]]$state)
    }
  }
  for(j in nx) {
    score <- family$score[[j]](y, family$map2par(eta), id = j)
    for(sj in names(x[[j]]$smooth.construct)) {
      tgrad <- x[[j]]$smooth.construct[[sj]]$grad(score, x[[j]]$smooth.construct[[sj]]$state$parameters, full = FALSE)
      grad <- c(grad, tgrad)
    }
  }
  return(grad)
}


## Optimizer based on optim().
opt <- function(x, y, family, start = NULL, verbose = TRUE, digits = 3,
  gradient = TRUE, hessian = FALSE, eps = .Machine$double.eps^0.5, maxit = 100, ...)
{
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("design construct names mismatch with family names!")

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, ...)

  if(!is.null(start))
    x <- set.starting.values(x, start)

  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  par <- get.all.par(x, list = FALSE, drop = TRUE)

  if(!hessian) {
    if(verbose)
      bamlss_log_posterior_iteration <<- 1

    opt <- optim(par, fn = log_posterior,
      gr = if(!is.null(family$score) & gradient) grad_posterior else NULL,
      x = x, y = y, family = family, method = "BFGS", verbose = verbose,
      digits = digits, control = list(fnscale = -1, reltol = eps, maxit = maxit),
      hessian = TRUE)
 
    if(verbose) {
      cat("\n")
      rm(bamlss_log_posterior_iteration, envir = .GlobalEnv)
    }

    eta <- get.eta.par(opt$par, x)

    return(list("parameters" = opt$par, "fitted.values" = eta,
      "logPost" = opt$value, "logLik" = family$loglik(y, family$map2par(eta)),
      "hessian" = opt$hessian, "converged" = opt$convergence < 1))
  } else {
    fn <- if(is.null(family$p2d)) {
      log_posterior
    } else function(par, ...) { sum(family$p2d(par, log = TRUE), na.rm = TRUE) }

    opt <- optimHess(par, fn = fn,
      gr = if(!is.null(family$score) & gradient & is.null(family$p2d)) grad_posterior else NULL,
      x = x, y = y, family = family, verbose = verbose, digits = digits,
      control = list(fnscale = -1, reltol = eps, maxit = maxit))
    return(opt)
  }
}


## Fast computation of weights and residuals when binning.
xbin.fun <- function(ind, weights, e, xweights, xrres, oind)
{
  .Call("xbin_fun", as.integer(ind), as.numeric(weights), 
    as.numeric(e), as.numeric(xweights), as.numeric(xrres),
    as.integer(oind))
  invisible(NULL)
}


## Likelihood based boosting.
boost_logLik <- function(x, y, family, weights = NULL, offset = NULL,
  criterion = c("AICc", "BIC", "AIC"),
  nu = 1, df = 4, maxit = 100, mstop = NULL, best = TRUE,
  verbose = TRUE, digits = 4,
  eps = .Machine$double.eps^0.25, plot = TRUE, ...)
{
  if(!is.null(mstop))
    maxit <- mstop

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, ...)

  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  criterion <- match.arg(criterion)

  np <- length(nx)
  y <- y[[1]]
  nobs <- if(is.null(dim(y))) length(y) else nrow(y)

  ## Setup boosting structure, i.e, all parametric
  ## terms get an entry in $smooth.construct object.
  ## Intercepts are initalized.
  x <- boost.transform(x, y, df, family, weights, offset, maxit, eps, ...)

  ## Create a list() that saves the states for
  ## all parameters and model terms.
  states <- make.state.list(x)

  ## Term selector help vectors.
  select <- unlist(states)
  term.names <- names(select)

  ## Extract actual predictor.
  eta <- get.eta(x)

  ## Initial parameters.
  parameters <- get.all.par(x)

  ## Start boosting.
  eps0 <- 1; iter <- 1
  save.ic <- save.ll <- NULL
  ll <- family$loglik(y, family$map2par(eta))
  while(iter <= maxit) {
    eta0 <- eta

    ## Cycle through all parameters
    for(i in nx) {
      peta <- family$map2par(eta)

      ## Compute weights.
      hess <- family$hess[[i]](y, peta, id = i)

      ## Score.
      score <- family$score[[i]](y, peta, id = i)

      ## Compute working observations.
      z <- eta[[i]] + 1 / hess * score

      ## Residuals.
      resids <- z - eta[[i]]

      for(j in names(x[[i]]$smooth.construct)) {
        ## Get updated parameters.
        states[[i]][[j]] <- boost_iwls(x[[i]]$smooth.construct[[j]], hess, resids, nu)

        ## Compute likelihood contribution.
        eta[[i]] <- eta[[i]] + fitted(states[[i]][[j]])
        select[paste(i, j, sep = ".")] <- -1 * (ll - family$loglik(y, family$map2par(eta)))
        eta[[i]] <- eta0[[i]]
      }
    }

    ## Which term to update.
    take <- strsplit(term.names[which.max(select)], ".", fixed = TRUE)[[1]]

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
      cat(if(interactive()) "\r" else "\n")
      vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
        " logLik ", fmt(ll, width = 8, digits = digits),
        " edf ", fmt(edf, width = 6, digits = digits),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
      cat(vtxt)

      if(.Platform$OS.type != "unix" & interactive()) flush.console()
    }

    iter <- iter + 1
  }

  if(verbose) cat("\n")

  mstop <- which.min(save.ic)

  ## Overwrite parameter state.
  if(best) {
    for(i in nx) {
      rn <- NULL
      for(j in names(x[[i]]$smooth.construct)) {
        x[[i]]$smooth.construct[[j]]$state$parameters <- parameters[[i]]$s[[j]]
        g <- get.par(x[[i]]$smooth.construct[[j]]$state$parameters, "b")
        if(!is.null(tau2 <- attr(parameters[[i]]$s[[j]], "true.tau2"))) {
          x[[i]]$smooth.construct[[j]]$state$parameters <- set.par(x[[i]]$smooth.construct[[j]]$state$parameters,
            tau2, "tau2")
        }
        x[[i]]$smooth.construct[[j]]$state$edf <- attr(parameters[[i]]$s[[j]], "edf")
        if(is.null(x[[i]]$smooth.construct[[j]]$state$edf))
          x[[i]]$smooth.construct[[j]]$state$edf <- 0
      }
    }
  }

  if(verbose) {
    cat("---\n", criterion, "=", save.ic[mstop], "-> at mstop =", mstop, "\n---\n")
  }

  bsum <- boost.summary(x, mstop, criterion, save.ic)
  x <- boost.retransform(x)

  list("parameters" = get.all.par(x), "fitted.values" = get.eta(x), "boost.summary" = bsum)
}


## 2nd booster.
boost <- function(x, y, family, weights = NULL, offset = NULL,
  nu = 0.1, df = 4, maxit = 100, mstop = NULL, best = TRUE,
  verbose = TRUE, digits = 4,
  eps = .Machine$double.eps^0.25, plot = TRUE, ...)
{
  if(is.null(family$score))
    stop("need score functions in family object for boosting!")


  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")

  if(!is.null(mstop))
    maxit <- mstop

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, ...)

  np <- length(nx)
  nobs <- nrow(y)
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  ## Setup boosting structure, i.e, all parametric
  ## terms get an entry in $smooth.construct object.
  ## Intercepts are initalized.
  x <- boost.transform(x, y, df, family, weights, offset, maxit, eps, ...)

  ## Create a list() that saves the states for
  ## all parameters and model terms.
  states <- make.state.list(x)

  ## Term selector help vectors.
  select <- rep(NA, length = length(nx))
  names(select) <- nx
  loglik <- select

  ## Save rss in list().
  rss <- make.state.list(x, type = 2)

  ## Extract actual predictor.
  eta <- get.eta(x)

  ## Start boosting.
  eps0 <- 1; iter <- 1
  save.ll <- NULL
  ll <- family$loglik(y, family$map2par(eta))
  while(iter <= maxit) {
    eta0 <- eta

    ## Cycle through all parameters
    for(i in nx) {
      peta <- family$map2par(eta)

      ## Actual gradient.
      grad <- family$score[[i]](y, peta, id = i)

      ## Fit to gradient.
      for(j in names(x[[i]]$smooth.construct)) {
        ## Get updated parameters.
        states[[i]][[j]] <- boost_fit(x[[i]]$smooth.construct[[j]], grad, nu)

        ## Get rss.
        rss[[i]][j] <- states[[i]][[j]]$rss
      }

      ## Which one is best?
      select[i] <- which.min(rss[[i]])

      ## Compute likelihood contribution.
      eta[[i]] <- eta[[i]] + fitted(states[[i]][[j]])
      loglik[i] <- -1 * (ll - family$loglik(y, family$map2par(eta)))
      eta[[i]] <- eta0[[i]]
    }

    i <- which.max(loglik)

    ## Which term to update.
    take <- c(nx[i], names(rss[[i]])[select[i]])

    ## Update selected base learner.
    eta[[take[1]]] <- eta[[take[1]]] + fitted(states[[take[1]]][[take[2]]])

    ## Write to x.
    x[[take[1]]]$smooth.construct[[take[2]]]$state <- increase(x[[take[1]]]$smooth.construct[[take[2]]]$state, states[[take[1]]][[take[2]]])
    x[[take[1]]]$smooth.construct[[take[2]]]$selected[iter] <- 1
    x[[take[1]]]$smooth.construct[[take[2]]]$loglik[iter] <- loglik[i]

    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

    peta <- family$map2par(eta)
    ll <- family$loglik(y, peta)

    save.ll <- c(save.ll, ll)

    if(verbose) {
      cat(if(interactive()) "\r" else "\n")
      vtxt <- paste(
        " logLik ", fmt(ll, width = 8, digits = digits),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
      cat(vtxt)

      if(.Platform$OS.type != "unix" & interactive()) flush.console()
    }

    iter <- iter + 1
  }

  if(verbose) cat("\n")

  bsum <- boost.summary(x, maxit, "logLik", save.ll)
  x <- boost.retransform(x)

  list("parameters" = get.all.par(x), "fitted.values" = get.eta(x), "boost.summary" = bsum)
}



## Boost setup.
boost.transform <- function(x, y, df, family, weights = NULL, offset = NULL,
  maxit = 100, eps = .Machine$double.eps^0.25, ...)
{
  np <- length(x)
  nx <- names(x)

  ## Initialize select indicator and intercepts.
  eta <- get.eta(x)
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      x[[nx[j]]]$smooth.construct[[sj]] <- assign.df(x[[nx[j]]]$smooth.construct[[sj]], df)
      if(!x[[nx[j]]]$smooth.construct[[sj]]$fxsp & !x[[nx[j]]]$smooth.construct[[sj]]$fixed) {
        x[[nx[j]]]$smooth.construct[[sj]]$old.optimize <- x[[nx[j]]]$smooth.construct[[sj]]$state$do.optim
        x[[nx[j]]]$smooth.construct[[sj]]$state$do.optim <- FALSE
        x[[nx[j]]]$smooth.construct[[sj]]$do.optim <- FALSE
      }
    }
    if(has_pterms(x[[nx[j]]]$terms)) {
      ii <- which(names(x[[nx[j]]]$smooth.construct) == "model.matrix")
      model.matrix <- list()
      cn <- colnames(x[[nx[j]]]$smooth.construct[[ii]]$X)
      g0 <- get.par(x[[nx[j]]]$smooth.construct[[ii]]$state$parameters, "b")
      nm <- NULL
      for(pj in 1:ncol(x[[nx[j]]]$smooth.construct[[ii]]$X)) {
        model.matrix[[pj]] <- list()
        model.matrix[[pj]]$label <- cn[pj]
        model.matrix[[pj]]$term <- cn[pj]
        model.matrix[[pj]]$X <- x[[nx[j]]]$smooth.construct[[ii]]$X[, pj, drop = FALSE]
        if(cn[pj] != "(Intercept)") {
          mx <- mean(model.matrix[[pj]]$X)
          sdx <- sd(model.matrix[[pj]]$X)
          model.matrix[[pj]]$X <- (model.matrix[[pj]]$X - mx) / sdx
          model.matrix[[pj]]$boost.scale <- list("mean" = mx, "sd" = sdx)
        }
        model.matrix[[pj]]$binning <- x[[nx[j]]]$smooth.construct[[ii]]$binning
        model.matrix[[pj]]$nobs <- x[[nx[j]]]$smooth.construct[[ii]]$nobs
        model.matrix[[pj]]$fixed <- TRUE
        model.matrix[[pj]]$fxsp <- FALSE
        model.matrix[[pj]]$weights <- x[[nx[j]]]$smooth.construct[[ii]]$weights
        model.matrix[[pj]]$rres <- x[[nx[j]]]$smooth.construct[[ii]]$rres
        model.matrix[[pj]]$fit.reduced <- x[[nx[j]]]$smooth.construct[[ii]]$fit.reduced
        model.matrix[[pj]]$fit.fun <- x[[nx[j]]]$smooth.construct[[ii]]$fit.fun
        model.matrix[[pj]]$state <- list("parameters" = g0[pj])
        model.matrix[[pj]]$state$fitted.values <- drop(model.matrix[[pj]]$X %*% g0[pj])
        model.matrix[[pj]]$state$edf <- 0
        model.matrix[[pj]]$state$do.optim <- FALSE
        model.matrix[[pj]]$is.model.matrix <- TRUE
        model.matrix[[pj]]$selected <- rep(0, length = maxit)
        model.matrix[[pj]]$upper <- Inf
        model.matrix[[pj]]$lower <- -Inf
        class(model.matrix[[pj]]) <- class(x[[nx[j]]]$smooth.construct[[ii]])
      }
      names(model.matrix) <- cn
      x[[nx[j]]]$smooth.construct[[ii]] <- NULL
      x[[nx[j]]]$smooth.construct <- c(model.matrix, x[[nx[j]]]$smooth.construct)
    }
  }

  ## Save more info.
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      x[[nx[j]]]$smooth.construct[[sj]]$selected <- rep(0, length = maxit)
      x[[nx[j]]]$smooth.construct[[sj]]$loglik <- rep(0, length = maxit)
      x[[nx[j]]]$smooth.construct[[sj]]$state$edf <- 0
    }
  }

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
            family, y, eta, nx[j], weights = weights[[nx[j]]], offset = offset[[nx[j]]], ...)

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

  return(x)
}


## Simple list() generator for
## saving states of model terms.
make.state.list <- function(x, type = 1)
{
  elmts <- c("formula", "fake.formula")
  if(all(elmts %in% names(x))) {
    rval <- list()
    if(!is.null(x$model.matrix))
      rval$model.matrix <- NA
    if(!is.null(x$smooth.construct)) {
      for(j in names(x$smooth.construct))
        rval[[j]] <- NA
    }
    if(type > 1)
      rval <- unlist(rval)
  } else {
    rval <- list()
    for(j in names(x)) {
      rval[[j]] <- make.state.list(x[[j]], type)
    }
  }
  return(rval)
}


#	for(i in seq_along(coefs))
#	{
#		curname<-names(coefs)[i]
#		cpos<-match(curname, colnames(dfr))
#		usedScale<-attr(dfr[[cpos]], "scaled:scale")
#		usedCenter<-attr(dfr[[cpos]], "scaled:center")
#		if(! is.null(usedScale))
#		{
#			catwif(verbosity > 0, "Scaling back for variable", curname)
#			catwif(verbosity >1, "usedScale structure")
#			if(verbosity > 1) str(usedScale)
#			catwif(verbosity >1, "usedCenter structure")
#			if(verbosity > 1) str(usedCenter)
#			oldcoef<-coefs[i]
#			itc<-itc - ((oldcoef * usedCenter)/usedScale)
#			coefs[i]<-oldcoef / usedScale
#		}
#	}
#	coefs<-c(itc, coefs)
#	names(coefs)[1]<-itcname
#	return(coefs)


#          if(!is.null(x[[i]]$smooth.construct[[j]]$boost.scale)) {
#            mx <- x[[i]]$smooth.construct[[j]]$boost.scale$mean
#            sdx <- x[[i]]$smooth.construct[[j]]$boost.scale$sd
#            x[[i]]$smooth.construct[[j]]$X <- x[[i]]$smooth.construct[[j]]$X * sdx + mx
#            if(!is.null(intercept))
#              intercept <- intercept - ((b * mx) / sdx)
#            b <- b / sdx
#          }


## Retransform 'x' to 'bamlss.frame' structure.
boost.retransform <- function(x) {
  for(i in names(x)) {
    if(has_pterms(x[[i]]$terms)) {
      state <- list()
      X <- drop <- xscales <- NULL
      for(j in names(x[[i]]$smooth.construct)) {
        if(inherits(x[[i]]$smooth.construct[[j]], "model.matrix")) {
          drop <- c(drop, j)
          b <- get.par(x[[i]]$smooth.construct[[j]]$state$parameters, "b")
          X <- cbind(X, x[[i]]$smooth.construct[[j]]$X)
          state$parameters <- c(state$parameters, b)
        }
      }
      label <- paste(drop, collapse = "+")
      binning <- x[[i]]$smooth.construct[[drop[1]]]$binning
      state$fitted.values <- drop(X %*% state$parameters)
      x[[i]]$smooth.construct[drop] <- NULL
      x[[i]]$smooth.construct$model.matrix <- list(
        "X" = X,
        "S" = list(diag(0, ncol(X))),
        "rank" = ncol(X),
        "term" = label,
        "label" = label,
        "bs.dim" = ncol(X),
        "fixed" = TRUE,
        "is.model.matrix" = TRUE,
        "by" = "NA",
        "xt" = list("binning" = binning),
        "state" = state
      )
      x[[i]]$smooth.construct$model.matrix$fit.fun <- make.fit.fun(x[[i]]$smooth.construct$model.matrix)
      
    }
  }
  return(x)
}


## Boosting iwls.
boost_iwls <- function(x, hess, resids, nu)
{
  ## Initial parameters and fit.
  g0 <- get.par(x$state$parameters, "b")
  fit0 <- fitted(x$state)

  ## Compute reduced residuals.
  xbin.fun(x$binning$sorted.index, hess, resids, x$weights, x$rres, x$binning$order)

  ## Compute mean and precision.
  XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)
  if(x$fixed) {
    P <- matrix_inv(XWX, index = x$sparse.setup)
  } else {
    S <- 0
    tau2 <- get.state(x, "tau2")
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    P <- matrix_inv(XWX + S, index = x$sparse.setup)
  }

  ## New parameters.
  g <- nu * drop(P %*% crossprod(x$X, x$rres))

  ## Finalize.
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))

  ## Find edf.
  xbin.fun(x$binning$sorted.index, hess, resids + fit0 + fitted(x$state), x$weights, x$rres, x$binning$order)

  XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)
  if(x$fixed) {
    P <- matrix_inv(XWX, index = x$sparse.setup)
  } else {
    g0 <- g0 + g

    objfun <- function(tau2) {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      P <- matrix_inv(XWX + S, index = x$sparse.setup)
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
    P <- matrix_inv(XWX + S, index = x$sparse.setup)
  }

  ## Assign degrees of freedom.
  x$state$edf <- sum.diag(XWX %*% P)
  attr(x$state$parameters, "edf") <- x$state$edf

  return(x$state)
}


## Boosting gradient fit.
boost_fit <- function(x, y, nu)
{
  ## Compute reduced residuals.
  xbin.fun(x$binning$sorted.index, rep(1, length = length(y)), y, x$weights, x$rres, x$binning$order)

  ## Compute mean and precision.
  XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)
  if(x$fixed) {
    P <- matrix_inv(XWX, index = x$sparse.setup)
  } else {
    S <- 0
    tau2 <- get.state(x, "tau2")
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    P <- matrix_inv(XWX + S, index = x$sparse.setup)
  }

  ## New parameters.
  g <- nu * drop(P %*% crossprod(x$X, x$rres))

  ## Finalize.
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
  x$state$rss <- sum((x$state$fitted.values - y)^2)

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


## Extract summary for boosting.
boost.summary <- function(x, mstop, criterion, save.ic)
{
  nx <- names(x)
  labels <- NULL
  ll.contrib <- NULL
  bsum <- lmat <- list()
  for(i in nx) {
    rn <- NULL
    for(j in names(x[[i]]$smooth.construct)) {
      labels <- c(labels, paste(x[[i]]$smooth.construct[[j]]$label, i, sep = "."))
      rn <- c(rn, x[[i]]$smooth.construct[[j]]$label)
      bsum[[i]] <- rbind(bsum[[i]], sum(x[[i]]$smooth.construct[[j]]$selected[1:mstop]) / mstop * 100)
      lmat[[i]] <- rbind(lmat[[i]], sum(x[[i]]$smooth.construct[[j]]$loglik[1:mstop]))
      ll.contrib <- cbind(ll.contrib, cumsum(x[[i]]$smooth.construct[[j]]$loglik))
    }
    if(!is.matrix(bsum[[i]])) bsum[[i]] <- matrix(bsum[[i]], nrow = 1)
    bsum[[i]] <- cbind(bsum[[i]], lmat[[i]])
    if(!is.matrix(bsum[[i]])) bsum[[i]] <- matrix(bsum[[i]], nrow = 1)
    colnames(bsum[[i]]) <- c(paste(i, "% selected"), "LogLik contrib.")
    rownames(bsum[[i]]) <- rownames(lmat[[i]]) <- rn
    bsum[[i]] <- bsum[[i]][order(bsum[[i]][, 2], decreasing = TRUE), , drop = FALSE]
  }
  colnames(ll.contrib) <- labels
  names(bsum) <- nx
  bsum <- list("summary" = bsum, "mstop" = mstop, "criterion" = criterion,
    "ic" = save.ic, "loglik" = ll.contrib)
  class(bsum) <- "boost.summary"
  return(bsum)
}


## Smallish print function for boost summaries.
print.boost.summary <- function(object, summary = TRUE, plot = FALSE, ...)
{
  if(summary) {
    np <- length(object$summary)
    cat("\n")
    cat(object$criterion, "=", object$ic[object$mstop], "-> at mstop =", object$mstop, "\n---\n")
    for(j in 1:np) {
      if(length(object$summary[[j]]) < 2) {
        print(round(object$summary[[j]], digits = 4))
      } else printCoefmat(object$summary[[j]], digits = 4)
      if(j != np)
        cat("---\n")
    }
    cat("\n")
  }

  if(plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(1, 2), mar = c(5.1, 4.1, 2.1, 2.1))
    plot(object$ic, type = "l", xlab = "Iteration", ylab = object$criterion)
    abline(v = object$mstop, lwd = 3, col = "lightgray")
    axis(3, at = object$mstop, labels = paste("mstop =", object$mstop))
    par(mar = c(5.1, 4.1, 2.1, 10.1))
    matplot(1:nrow(object$loglik), object$loglik, type = "l", lty = 1,
      xlab = "Iteration", ylab = "LogLik contribution", col = "black")
    abline(v = object$mstop, lwd = 3, col = "lightgray")
    axis(4, at = object$loglik[nrow(object$loglik), ], labels = colnames(object$loglik), las = 1)
    axis(3, at = object$mstop, labels = paste("mstop =", object$mstop))
  }

  return(invisible(object))
}

plot.boost.summary <- function(x, ...)
{
  print.boost.summary(x, summary = FALSE, plot = TRUE) 
}


## Assign starting values.
set.starting.values <- function(x, start)
{
  if(!is.null(start)) {
    if(is.list(start)) {
      if("parameters" %in% names(start))
        start <- start$parameters
    }
    if(is.list(start))
      start <- unlist(start)
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
              i <- grep2(c(".edf", ".accepted", ".alpha"), names(tpar))
              x[[id]]$smooth.construct$model.matrix$state$parameters <- if(length(i)) tpar[-i] else tpar
              x[[id]]$smooth.construct$model.matrix$state$fitted.values <- x[[id]]$smooth.construct$model.matrix$fit.fun(x[[id]]$smooth.construct$model.matrix$X, x[[id]]$smooth.construct$model.matrix$state$parameters)
            }
          }
        }
        for(j in seq_along(x[[id]]$smooth.construct)) {
          tl <- x[[id]]$smooth.construct[[j]]$label
          take <- grep(tl <- paste(id, "s", tl, sep = "."),
            names(start), fixed = TRUE, value = TRUE)
          if(x[[id]]$smooth.construct[[j]]$by == "NA") {
            take <- take[!grepl(paste(tl, ":", sep = ""), take, fixed = TRUE)]
          }
          if(length(take)) {
            tpar <- start[take]
            names(tpar) <- gsub(paste(tl, ".", sep = ""), "", names(tpar), fixed = TRUE)
            i <- grep2(c(".edf", ".accepted", ".alpha"), names(tpar))
            tpar <- if(length(i)) tpar[-i] else tpar
            spar <- x[[id]]$smooth.construct[[j]]$state$parameters
            spar <- set.par(spar, get.par(tpar, "b"), "b")
            if(any(grepl("tau2", names(tpar)))) {
              spar <- set.par(spar, get.par(tpar, "tau2"), "tau2")
            }
            x[[id]]$smooth.construct[[j]]$state$parameters <- spar
            x[[id]]$smooth.construct[[j]]$state$fitted.values <- x[[id]]$smooth.construct[[j]]$fit.fun(x[[id]]$smooth.construct[[j]]$X, x[[id]]$smooth.construct[[j]]$state$parameters)
          }
        }
      }
    }
  }
  return(x)
}

