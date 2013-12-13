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
          "term" = "parametric",
          "bs.dim" = ncol(obj$X),
          "fixed" = TRUE
        )
        obj$sterms <- c(obj$strems, "parametric")
        obj$X <- NULL
      }
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          obj$smooth[[j]] <- smooth.IWLS(obj$smooth[[j]])
        }
      }
    }
    obj
  }

  x <- tIWLS(x, ...)

  attr(x, "call") <- call
  attr(x, "response.vec") <- attr(x, "model.frame")[, attr(attr(x, "model.frame"), "response.name")]
  attr(x, "model.frame") <- NULL

  x
}


## Function to setup IWLS smooths.
smooth.IWLS <- function(x, ...) {
  UseMethod("smooth.IWLS")
}

smooth.IWLS.default <- function(x, ...)
{
  if(is.null(x$get.mu)) {
    x$get.mu <- if(is.null(x$xt$get.mu)) {
      function(X, g) {
        X %*% as.numeric(g)
      }
    } else x$xt$get.mu
  }

  x$a <- if(is.null(x$xt$a)) 1e-04 else x$xt$a
  x$b <- if(is.null(x$xt$b)) 1e-04 else x$xt$b

  if(is.null(x$state)) {
    x$p.save <- c("g", "tau2")
    x$state <- list()
    x$state$g <- runif(ncol(x$X), 1.001, 1.002)
    x$state$tau2 <- runif(1, 0.99, 1)
    x$s.colnames <- if(is.null(x$s.colnames)) {
      c(paste("g[", 1:length(x$state$g), "]", sep = ""), "tau2[1]")
    } else x$s.colnames
    x$np <- ncol(x$X) + 1
  }

  if(is.null(x$propose)) {
    if(is.null(x$rand)) {
      x$propose <- function(x, family, response, eta, id, ...) {
        ## Compute weights.
        weights <- family$weights[[id]](response, eta)

        ## Score.
        score <- family$score[[id]](response, eta)

        ## Compute working observations.
        z <- eta[[id]] + 1 / weights * score

        ## Save old predictor.
        eta0 <- eta[[id]]

        ## Compute old log likelihood and old log coefficients prior.
        pibeta <- family$loglik(response, eta)
        p1 <- drop(-1 / x$state$tau2 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g)

        ## Compute partial predictor.
        eta[[id]] <- eta0 - x$state$fit

        ## Compute mean and precision.
        XW <- t(x$X * weights)
        P <- chol2inv(chol(P0 <- XW %*% x$X + 1 / x$state$tau2 * x$S[[1]]))
        M <- P %*% (XW %*% (z - eta[[id]]))

        ## Save old coefficients
        g0 <- x$state$g

        ## Sample new parameters.
        x$state$g <- drop(rmvnorm(n = 1, mean = M, sigma = P))

        ## Compute log priors
        p2 <- drop(-1 / x$state$tau2 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g)
        g1 <- x$state$g - M
        qbetaprop <- 0.5 * sum(log((diag(chol(P0, symmetric = TRUE))^2))) -0.5 * crossprod(g1, P0) %*% g1
        ## qbetaprop2 <- dmvnorm(x$state$g, mean = M, sigma = P, log = TRUE)

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

        ## Compute partial predictor.
        eta[[id]] <- eta[[id]] - x$state$fit

        ## Compute mean and precision.
        XW <- t(x$X * weights)
        P2 <- chol2inv(chol(P0 <- XW %*% x$X + 1 / x$state$tau2 * x$S[[1]]))
        M2 <- P2 %*% (XW %*% (z - eta[[id]]))

        ## Get the log prior.
        qbeta <- dmvnorm(g0, mean = M2, sigma = P2, log = TRUE)
        g2 <- g0 - M2
        qbeta <- 0.5 * sum(log((diag(chol(P0, symmetric = TRUE))^2))) -0.5 * crossprod(g2, P0) %*% g1
        ## qbeta2 <- dmvnorm(g0, mean = M2, sigma = P2, log = TRUE)

        ## Sample variance parameter.
        a <- x$rank / 2 + x$a
        b <- 0.5 * crossprod(x$state$g, x$S[[1]]) %*% x$state$g + x$b
        x$state$tau2 <- 1 / rgamma(1, a, b)

        ## Compute acceptance probablity
        x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

        return(x$state)
      }
    }
  }

  x
}


## Sampler based on IWLS proposals.
samplerIWLS <- function(x, n.iter = 12000, thin = 10, burnin = 2000,
  verbose = TRUE, step = 100, tdir = NULL, ...)
{
  require("mvtnorm")

  family <- attr(x, "family")
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  response <- attr(x, "response.vec")

  ## Actual number of samples to save.
  if(n.iter < burnin) stop("argument burnin exceeds n.iter!")
  if(thin > (n.iter - burnin)) stop("argument thin is set too large!")
  iterthin <- as.integer(seq(burnin, n.iter, by = thin))
  n.save <- length(iterthin)
  
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

  ## Start sampling
  for(i in 1:n.iter) {
    if(save <- i %in% iterthin)
      js <- which(iterthin == i)
    
    ## Cycle through all parameters
    for(j in 1:np) {
      ## And all terms.
      for(sj in seq_along(x[[nx[j]]]$smooth)) {
        ## Get proposed states.
        p.state <- x[[nx[j]]]$smooth[[sj]]$propose(x[[nx[j]]]$smooth[[sj]], family, response, eta, nx[j])

        ## If accepted, set current state to proposed state.
        accepted <- if(is.na(p.state$alpha)) FALSE else log(runif(1)) <= p.state$alpha

        if(accepted) {
          eta[[nx[j]]] <- eta[[nx[j]]] - x[[nx[j]]]$smooth[[sj]]$state$fit + p.state$fit
          x[[nx[j]]]$smooth[[sj]]$state <- p.state 
        }

        ## Save the samples and acceptance.
        if(save) {
          x[[nx[j]]]$smooth[[sj]]$s.alpha[js] <- accepted
          x[[nx[j]]]$smooth[[sj]]$s.samples[js, ] <- unlist(x[[nx[j]]]$smooth[[sj]]$state[x[[nx[j]]]$smooth[[sj]]$p.save])
        }
      }
    }
    if(verbose) {
      if(i %% step == 0)
        cat("iteration:", i, "\n")
    }
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
      if(!is.null(dim(obj$X))) {
        fn <- file.path(tdir, paste(id, if(!is.null(id)) ":", "h", obj$hlevel, ":", "param", ".raw", sep = ""))
        colnames(obj$s.samples) <- paste("b[", 1:ncol(obj$X), "]", sep = "")
        write.table(obj$s.samples, file = fn, row.names = FALSE, quote = FALSE)
      }
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          fn <- file.path(tdir, paste(id, if(!is.null(id)) ":", "h",
            obj$hlevel, ":", obj$smooth[[j]]$label, ".raw", sep = ""))
          colnames(obj$smooth[[j]]$s.samples) <- obj$smooth[[j]]$s.colnames
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
    st <- as.matrix(read.table(file.path(tdir, j), header = TRUE))
    colnames(st) <- paste(gsub(".raw", "", j), colnames(st), sep = ":")
    samples <- cbind(samples, st)
  }

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

  }
}

