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
bamlss.setup <- function(x, update = "iwls", do.optim = NULL, criterion = c("AICc", "BIC", "AIC"),
  nu = 0.1, coefficients = NULL, ...)
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

  foo <- function(x, id = NULL) {
    if(!any(c("formula", "fake.formula", "response") %in% names(x))) {
      nx <- names(x)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(x)
      if(length(unique(nx)) < length(x)) nx <- 1:length(x)
      for(j in nx)
        x[[j]] <- foo(x[[j]], id = j)
    } else {
      if(is.null(id)) id <- ""
      if(!is.null(dim(x$X))) {
        if(nrow(x$X) > 0 & !is.na(mean(unlist(x$X), na.rm = TRUE))) {
          if(is.null(x$smooth)) x$smooth <- list()
          label <- if(is.null(colnames(x$X))) {
            paste("b", 1:ncol(x$X), sep = "", collapse = "+")
          } else paste(colnames(x$X), collapse = "+")
          x$smooth <- c(list("parametric" = list(
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
          )), x$smooth)
          if(!is.null(coefficients)) {
            if(any(grepl(id, names(coefficients)))) {
              label <- strsplit(x$smooth$parametric$label, "+", fixed = TRUE)[[1]]
              label <- paste(label, id, sep = ":")
              if(!all(label %in% names(coefficients))) {
                nlabel <- label[!(label %in% names(coefficients))]
                stop(paste("cannot find the following coefficients: ",
                  paste(nlabel, collapse = ", "), "!", sep = ""))
              }
              x$smooth$parametric$xt$state <- list("parameters" = coefficients[label])
              names(x$smooth$parametric$xt$state$parameters) <- paste("g", 1:length(x$smooth$parametric$xt$state$parameters), sep = "")
            } else stop(paste("coefficients for parameter", id, "are missing!"))
          }
          class(x$smooth[["parametric"]]) <- c(class(x$smooth[["parametric"]]), "no.mgcv", "parametric")
          x$sterms <- c(x$strems, "parametric")
          x$X <- NULL
        }
      }
      if(length(x$smooth)) {
        for(j in seq_along(x$smooth)) {
          if(!is.null(coefficients)) {
            if(!any(grepl(id, names(coefficients))))
              stop(paste("coefficients for parameter", id, "are missing!"))
            if(any(grepl(id, names(coefficients))) & is.null(x$smooth[[j]]$is.parametric)) {
              label <- x$smooth[[j]]$label
              if(any(take <- grepl(label, names(coefficients), fixed = TRUE))) {
                tcoef <- grep(id, names(coefficients)[take], fixed = TRUE, value = TRUE)
                tcoef <- tcoef[!grepl("edf", tcoef, fixed = TRUE)]
                tcoef <- tcoef[!grepl("alpha", tcoef, fixed = TRUE)]
                tau2 <- tcoef[grepl("tau2", tcoef, fixed = TRUE)]
                tcoef <- tcoef[!grepl("tau2", tcoef, fixed = TRUE)]
                x$smooth[[j]]$xt$state <- list("parameters" = coefficients[tcoef])
                names(x$smooth[[j]]$xt$state$parameters) <- paste("g", 1:length(x$smooth[[j]]$xt$state$parameters), sep = "")
                if(length(tau2)) {
                  tau2 <- coefficients[tau2]
                  names(tau2) <- paste("tau2", 1:length(tau2), sep = "")
                  x$smooth[[j]]$xt$state$parameters <- c(x$smooth[[j]]$xt$state$parameters, tau2)
                }
              } else stop(paste("cannot find coefficients for term ", label, " for parameter", id, "!", sep = ""))
            }
          }
          x$smooth[[j]] <- smooth.bamlss(x$smooth[[j]])
          if(!is.null(x$smooth[[j]]$xt$update))
            x$smooth[[j]]$update <- x$smooth[[j]]$xt$update
          if(is.null(x$smooth[[j]]$update)) {
            if(is.character(update)) {
              if(!grepl("bfit0_", update))
                update <- paste("bfit0", update, sep = "_")
              update <- eval(parse(text = update))
            }
            x$smooth[[j]]$update <- update
          }
          if(!is.null(x$smooth[[j]]$do.optim))
            x$smooth[[j]]$state$optimize <- x$smooth[[j]]$do.optim
          if(is.null(x$smooth[[j]]$state$optimize)) {
            if(is.null(do.optim))
              x$smooth[[j]]$state$optimize <- TRUE
            else
              x$smooth[[j]]$state$optimize <- do.optim
          }
          if(!is.null(x$smooth[[j]]$rank))
            x$smooth[[j]]$rank <- as.numeric(x$smooth[[j]]$rank)
          x$smooth[[j]]$criterion <- criterion
          x$smooth[[j]]$nu <- nu
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
                "by" = "NA",
                "nu" = nu
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
  if(what %in% c("par", "parameters")) {
    return(x$state$parameters)
  } else {
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
  if(inherits(x, "special"))
    return(x)
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
      }
    }
  }
  if((ntau2 > 0) & !any(grepl("tau2", names(state$parameters))) & is.null(x$is.parametric)) {
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
      return(edf)
    }
  }
  if(is.null(x$grad)) {
    x$grad <- function(score = NULL, parameters, full = TRUE) {
      gamma <- get.par(parameters, "g")
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
  if(!inherits(x, "parametric"))
    colnames(x$X) <- NULL

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
  nobs <- if(is.null(dim(response))) length(response) else nrow(response)
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
get.ic <- function(family, response, eta, edf, n, type = c("AIC", "BIC", "AICc", "MP"), ...)
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
make_par <- function(x, type = 1, add.tau2 = FALSE) {
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
      if(!add.tau2) {
        tlower <- tlower[!grepl("tau2", names(tlower))]
        tupper <- tupper[!grepl("tau2", names(tupper))]
        tpar <- tpar[!grepl("tau2", names(tpar))]
      }
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
bfit0_newton <- function(x, family, response, eta, id, ...)
{
  args <- list(...)

  eta[[id]] <- eta[[id]] - fitted(x$state)

  tau2 <- if(!x$fixed) get.par(x$state$parameters, "tau2") else NULL

  lp <- function(g) {
    eta[[id]] <- eta[[id]] + x$get.mu(x$X, g)
    family$loglik(response, family$map2par(eta)) + x$prior(c(g, tau2))
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

  g <- get.par(x$state$parameters, "gamma")
  nu <- if(is.null(x$nu)) 0.1 else x$nu

  g.grad <- grad(fun = lp, theta = g, id = id, prior = NULL,
    args = list("gradient" = gfun, "x" = x, "y" = response, "eta" = eta))

  g.hess <- hess(fun = lp, theta = g, id = id, prior = NULL,
    args = list("gradient" = gfun, "hessian" = hfun, "x" = x, "y" = response, "eta" = eta))

  Sigma <- matrix_inv(g.hess)

  g <- drop(g + nu * Sigma %*% g.grad)

  x$state$parameters <- set.par(x$state$parameters, g, "g")
  x$state$fitted.values <- x$get.mu(x$X, get.state(x, "g"))
  x$state$hessian <- Sigma

  return(x$state)
}


bfit0_iwls <- function(x, family, response, eta, id, ...)
{
  args <- list(...)

  peta <- family$map2par(eta)
  if(is.null(args$weights)) {
    ## Compute weights.
    weights <- family$weights[[id]](response, peta, ...)
  } else weights <- args$weights

  if(is.null(args$z)) {
    ## Score.
    score <- family$score[[id]](response, peta, ...)

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
      eta2[[id]] <- eta2[[id]] + fit
      IC <- get.ic(family, response, family$map2par(eta2), edf0 + edf, length(z), x$criterion, ...)
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

  return(x$state)
}


## Updating based on optim.
bfit0_optim <- function(x, family, response, eta, id, ...)
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
    val <- -1 * (ll + lp)
    if(!is.finite(val)) val <- NA
    val
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
      IC <- get.ic(family, response, family$map2par(eta2), edf0 + edf, length(eta2[[id]]), type = x$criterion, ...)
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


set.all.par <- function(par, x, hessian = NULL)
{
  if(!is.null(hessian))
    hessian <- matrix_inv(-1 * hessian)
  nx <- names(x)
  np <- length(x)
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      tpar <- par[grep(paste("p", j, ".t", sj, ".", sep = ""), names(par), fixed = TRUE)]
      x[[nx[j]]]$smooth[[sj]]$state$parameters <- set.par(x[[nx[j]]]$smooth[[sj]]$state$parameters, get.par(tpar, "g"), "g")
      if(any(grepl("tau2", tpar))) {
        x[[nx[j]]]$smooth[[sj]]$state$parameters <- set.par(x[[nx[j]]]$smooth[[sj]]$state$parameters, get.par(tpar, "tau2"), "tau2")
      }
      x[[nx[j]]]$smooth[[sj]]$state$fitted.values <- x[[nx[j]]]$smooth[[sj]]$get.mu(x[[nx[j]]]$smooth[[sj]]$X,
        get.par(tpar, "g"))
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
        x[[nx[j]]]$smooth[[sj]]$state$hessian <- sigma
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
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      g <- get.state(x[[nx[j]]]$smooth[[sj]], "g")
      end <- length(g) + start - 1
      i <- start:end
      tpar <- par[i]
      names(tpar) <- names(g)
      start <- max(i) + 1
      x[[nx[j]]]$smooth[[sj]]$state$parameters <- set.par(x[[nx[j]]]$smooth[[sj]]$state$parameters, tpar, "g")
      x[[nx[j]]]$smooth[[sj]]$state$fitted.values <- x[[nx[j]]]$smooth[[sj]]$get.mu(x[[nx[j]]]$smooth[[sj]]$X,
        get.par(tpar, "g"))
      eta[[nx[j]]] <- eta[[nx[j]]] + fitted(x[[nx[j]]]$smooth[[sj]]$state)
      lprior <- lprior + x[[nx[j]]]$smooth[[sj]]$prior(c(tpar, get.state(x[[nx[j]]]$smooth[[sj]], "tau2")))
      edf <- edf + x[[nx[j]]]$smooth[[sj]]$edf(x[[nx[j]]]$smooth[[sj]])
    }
  }
  ll <- family$loglik(attr(x, "response.vec"), family$map2par(eta))
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
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      tpar <- par[grep(paste("p", j, ".t", sj, ".", sep = ""), names(par), fixed = TRUE)]
      x[[nx[j]]]$smooth[[sj]]$state$parameters <- set.par(x[[nx[j]]]$smooth[[sj]]$state$parameters, tpar, "g")
      x[[nx[j]]]$smooth[[sj]]$state$fitted.values <- x[[nx[j]]]$smooth[[sj]]$get.mu(x[[nx[j]]]$smooth[[sj]]$X,
        get.par(tpar, "g"))
      eta[[nx[j]]] <- eta[[nx[j]]] + fitted(x[[nx[j]]]$smooth[[sj]]$state)
    }
  }
  for(j in 1:np) {
    if(!is.null(family$gradient[[nx[j]]])) {
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        tp <- get.state(x[[nx[j]]]$smooth[[sj]], "par")
        tgrad <- drop(family$gradient[[nx[j]]](get.par(tp, "g"), attr(x, "response.vec"), eta, x[[nx[j]]]$smooth[[sj]]))
        if(!is.null(x[[nx[j]]]$smooth[[sj]]$grad)) {
          tgrad <- tgrad + drop(x[[nx[j]]]$smooth[[sj]]$grad(score = NULL, tp, full = FALSE))
        }
        grad <- c(grad, tgrad)
      }
    } else {
      score <- family$score[[nx[j]]](attr(x, "response.vec"), family$map2par(eta))
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        tgrad <- x[[nx[j]]]$smooth[[sj]]$grad(score, x[[nx[j]]]$smooth[[sj]]$state$parameters, full = FALSE)
        grad <- c(grad, tgrad)
      }
    }
  }
  return(grad)
}


opt0 <- function(x, verbose = TRUE, digits = 3, gradient = TRUE, hessian = FALSE,
  eps = .Machine$double.eps^0.5, maxit = 100, ...)
{
  x <- bamlss.setup(x, ...)

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

  x <- bamlss.setup(x, criterion = criterion, ...)

  family <- attr(x, "family")
  nx <- family$names[1:2]
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  criterion <- match.arg(criterion)

  np <- length(nx)
  response <- attr(x, "response.vec")
  nobs <- if(is.null(dim(response))) length(response) else nrow(response)
  eta <- get.eta(x)
  edf <- get.edf(x)

  ## Start the backfitting algorithm.
  eps0 <- eps + 1; iter <- 1
  while(eps0 > eps & iter < maxit) {
    eta0 <- eta
    ## Cycle through all parameters
    for(j in 1:np) {
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        ## Get updated parameters.
        p.state <- x[[nx[j]]]$smooth[[sj]]$update(x[[nx[j]]]$smooth[[sj]],
          family, response, eta, nx[j], edf = edf)

        ## Compute equivalent degrees of freedom.
        edf <- edf - x[[nx[j]]]$smooth[[sj]]$state$edf + p.state$edf

        ## Update predictor and smooth fit.
        eta[[nx[j]]] <- eta[[nx[j]]] - fitted(x[[nx[j]]]$smooth[[sj]]$state) + fitted(p.state)

        x[[nx[j]]]$smooth[[sj]]$state <- p.state
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

  if(verbose) cat("\n")

  return(x)
}


## Likelihood based boosting.
boost0 <- function(x, criterion = c("AICc", "BIC", "AIC"),
  nu = 1, maxit = 400, mstop = NULL,
  verbose = TRUE, digits = 4, tau2 = 10, ...)
{
  if(!is.null(mstop))
    maxit <- mstop

  criterion <- match.arg(criterion)

  x <- bamlss.setup(x, criterion = criterion, ...)

  family <- attr(x, "family")
  nx <- family$names[1:2]
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  criterion <- match.arg(criterion)

  np <- length(nx)
  response <- attr(x, "response.vec")
  nobs <- if(is.null(dim(response))) length(response) else nrow(response)

  ## Initialize select indicator and intercepts.
  states <- save.parametric <- list()
  eta <- get.eta(x)
  for(j in 1:np) {
    states[[j]] <- list()
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      if(is.null(x[[nx[j]]]$smooth[[sj]]$sp) & !x[[nx[j]]]$smooth[[sj]]$fixed) {
        x[[nx[j]]]$smooth[[sj]]$old.optimize <- x[[nx[j]]]$smooth[[sj]]$state$optimize
        x[[nx[j]]]$smooth[[sj]]$state$optimize <- FALSE
        x[[nx[j]]]$smooth[[sj]]$optimize <- FALSE
      }
      if(inherits(x[[nx[j]]]$smooth[[sj]], "parametric")) {
        parametric <- list()
        cn <- colnames(x[[nx[j]]]$smooth[[sj]]$X)
        g0 <- get.par(x[[nx[j]]]$smooth[[sj]]$state$parameters, "g")
        for(pj in 1:ncol(x[[nx[j]]]$smooth[[sj]]$X)) {
          parametric[[pj]] <- list()
          parametric[[pj]]$label <- cn[pj]
          parametric[[pj]]$term <- cn[pj]
          parametric[[pj]]$X <- x[[nx[j]]]$smooth[[sj]]$X[, pj, drop = FALSE]
          parametric[[pj]]$xbin.ind <- x[[nx[j]]]$smooth[[sj]]$xbin.ind
          parametric[[pj]]$xbin.sind <- x[[nx[j]]]$smooth[[sj]]$xbin.sind
          parametric[[pj]]$xbin.k <- x[[nx[j]]]$smooth[[sj]]$xbin.k
          parametric[[pj]]$xbin.order <- x[[nx[j]]]$smooth[[sj]]$xbin.order
          parametric[[pj]]$nobs <- x[[nx[j]]]$smooth[[sj]]$nobs
          parametric[[pj]]$fixed <- TRUE
          parametric[[pj]]$weights <- x[[nx[j]]]$smooth[[sj]]$weights
          parametric[[pj]]$rres <- x[[nx[j]]]$smooth[[sj]]$rres
          parametric[[pj]]$get.mu <- x[[nx[j]]]$smooth[[sj]]$get.mu
          parametric[[pj]]$edf <- function(X, g) { return(1) }
          parametric[[pj]]$state <- list("parameters" = g0[pj])
          names(parametric[[pj]]$state$parameters) <- "g1"
          parametric[[pj]]$state$fitted.values <- drop(parametric[[pj]]$X %*% g0[pj])
          parametric[[pj]]$state$edf <- 1
          parametric[[pj]]$is.parametric <- TRUE
          parametric[[pj]]$selected <- rep(0, length = maxit)
          parametric[[pj]]$upper <- Inf
          parametric[[pj]]$lower <- -Inf
        }
        save.parametric[[j]] <- x[[nx[j]]]$smooth[[sj]]
        x[[nx[j]]]$smooth[[sj]] <- NULL
        x[[nx[j]]]$smooth <- c(parametric, x[[nx[j]]]$smooth)
      }
      states[[j]][[sj]] <- list()
    }
  }
  sn <- NULL
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      sn <- c(sn, paste(j, sj, sep = "."))
      x[[nx[j]]]$smooth[[sj]]$selected <- rep(0, length = maxit)
      if(length(i <- grep("tau2", names(x[[nx[j]]]$smooth[[sj]]$state$parameters)))) {
        x[[nx[j]]]$smooth[[sj]]$state$parameters[i] <- rep(tau2, length = length(i))
      }
    }
  }
  select <- rep(0, length(sn))
  names(select) <- sn

  ## Start boosting.
  eps <-  .Machine$double.eps^0.25
  eps0 <- 1; iter <- 1
  save.ic <- save.ll <- NULL
  ll <- family$loglik(response, family$map2par(eta))
  while(iter < maxit) {
    eta0 <- eta

    ## Cycle through all parameters
    for(j in 1:np) {
      peta <- family$map2par(eta)

      ## Compute weights.
      weights <- family$weights[[nx[j]]](response, peta)

      ## Score.
      score <- family$score[[nx[j]]](response, peta)

      ## Compute working observations.
      z <- eta[[nx[j]]] + 1 / weights * score

      ## Residuals.
      resids <- z - eta[[nx[j]]]

      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        eta[[nx[j]]] <- eta[[nx[j]]] - fitted(x[[nx[j]]]$smooth[[sj]]$state)

        ## Get updated parameters.
        states[[j]][[sj]] <- boost0_iwls(x[[nx[j]]]$smooth[[sj]],
          family, response, eta, nx[j], weights, resids, nu)

        ## Compute likelihood contribution.
        eta[[nx[j]]] <- eta[[nx[j]]] + fitted(states[[j]][[sj]])
        select[paste(j, sj, sep = ".")] <- -2 * (ll - family$loglik(response, family$map2par(eta)))
        eta[[nx[j]]] <- eta0[[nx[j]]]
      }
    }

    ## Which term to update.
    take <- as.integer(strsplit(sn[which.max(select)], ".", fixed = TRUE)[[1]])

    ## Write to x.
    eta[[take[1]]] <- eta[[take[1]]] - fitted(x[[take[1]]]$smooth[[take[2]]]$state)
    x[[take[1]]]$smooth[[take[2]]]$state <- states[[take[1]]][[take[2]]]
    x[[take[1]]]$smooth[[take[2]]]$selected[iter] <- 1

    ## Update selected base learner.
    eta[[take[1]]] <- eta[[take[1]]] + fitted(x[[take[1]]]$smooth[[take[2]]]$state)

    edf <- get.edf(x)

    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

    peta <- family$map2par(eta)
    IC <- get.ic(family, response, peta, edf, nobs, criterion)
    ll <- family$loglik(response, peta)

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

  cat("\n")
  cat(criterion, "=", save.ic[mstop], "-> at mstop =", mstop, "\n---\n")
  cat("Frequencies\n---\n")
  labels <- NULL
  for(j in 1:np) {
    fmat <- rn <- NULL
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      labels <- c(labels, paste(x[[nx[j]]]$smooth[[sj]]$label, nx[j], sep = ":"))
      rn <- c(rn, x[[nx[j]]]$smooth[[sj]]$label)
      fmat <- rbind(fmat, sum(x[[nx[j]]]$smooth[[sj]]$selected) / maxit * 100)
    }
    rownames(fmat) <- rn
    colnames(fmat) <- paste(nx[j], "% selected")
    if(length(fmat) < 2) print(fmat) else printCoefmat(fmat, digits = 2)
    if(j != np)
      cat("---\n")
  }
  cat("\n")

  par(mfrow = c(2, 1), mar = c(4.1, 4.1, 4.1, 4.1))
  col <- rainbow_hcl(2)
  plot(save.ic, type = "l", xlab = "Iteration", ylab = criterion,
    main = paste("mstop =", mstop), col = col[1], lwd = 1.5)
  par(new = TRUE)
  plot(save.ll, type = "l", xlab = "iteration", ylab = "", col = col[2], lwd = 1.5, axes = FALSE)
  axis(4)
  box()
  mtext("Log. Lik", side = 4, line = 2.5)
  legend("right", c(criterion, "Log. Lik"), col = col, lwd = 1.5, bg = NA, box.col = NA)
  plot(c(0, length(select) + 1), c(0, maxit + 1), type = "n", axes = FALSE,
    xlab = "Term", ylab = "Iteration")
  axis(1, at = 1:length(select), labels = labels)
  i <- 1
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      if(length(yp <- which(x[[nx[j]]]$smooth[[sj]]$selected > 0))) {
        xp <- rep(i, length = length(yp))
        points(xp, yp, pch = 4)
      }
      i <- i + 1
    }
  }
  axis(2)
  box()

  ## Collect parametric effects.
  for(j in 1:np) {
    g <- drop <- NULL
    for(sj in seq_along(x[[nx[j]]]$smooth)) {
      if(!is.null(x[[nx[j]]]$smooth[[sj]]$is.parametric)) {
        g <- c(g, x[[nx[j]]]$smooth[[sj]]$state$parameters)
        drop <- c(drop, sj)
      }
    }
    if(!is.null(g)) {
      x[[nx[j]]]$smooth[drop] <- NULL
      names(g) <- paste("g", 1:length(g), sep = "")
      save.parametric[[j]]$state$fitted.values <- drop(save.parametric[[j]]$X %*% g)
      save.parametric[[j]]$state$parameters <- g
      x[[nx[j]]]$smooth <- c(list(save.parametric[[j]]), x[[nx[j]]]$smooth)
    }
  }

  return(x)
}


## Boosting iwls.
boost0_iwls <- function(x, family, response, eta, id, weights, resids, nu, ...)
{
  args <- list(...)

  ## Old parameters.
  g0 <- get.par(x$state$parameters, "gamma")

  ## Compute reduced residuals.
  xbin.fun(x$xbin.sind, weights, resids, x$weights, x$rres, x$xbin.order)

  ## Compute mean and precision.
  XWX <- crossprod(x$X, x$X * x$weights)
  if(x$fixed) {
    P <- matrix_inv(XWX)
  } else {
    S <- 0
    tau2 <- get.state(x, "tau2")
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    P <- matrix_inv(XWX + S)
  }

  ## New parameters
  g <- g0 + nu * drop(P %*% crossprod(x$X, x$rres))

  ## Finalize.
  if(!(any(is.na(g)) | any(g %in% c(-Inf, Inf)))) {
    x$state$parameters <- set.par(x$state$parameters, g, "g")
    x$state$fitted.values <- x$get.mu(x$X, get.state(x, "g"))
    x$state$edf <- sum(diag(P %*% XWX))
  }

  return(x$state)
}

