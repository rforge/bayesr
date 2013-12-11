## Setup for IWLS sampler, handling
## sampling functions.
setupIWLS <- function(x, ...)
{
  call <- x$call; x$call <- NULL

  sIWLS <- function(obj, ...) {
    if(!any(c("formula", "fake.formula", "response") %in% names(obj))) {
      nx <- names(obj)
      nx <- nx[nx != "call"]
      if(is.null(nx)) nx <- 1:length(obj)
      if(length(unique(nx)) < length(obj)) nx <- 1:length(obj)
      for(j in nx)
        obj[[j]] <- sIWLS(obj[[j]], ...)
    } else {
      if(!is.null(dim(obj$X))) {
        obj$p.save <- "b"
        obj$samples <- list("b" = runif(ncol(obj$X), 0.95, 1.05))
        obj$update = function(X, z, eta, weights) {
          XW <- crossprod(X, diag(weights))
          P <- chol2inv(chol(XW %*% X))
          M <- P %*% XW %*% (z - eta)
          drop(mvrnorm(n = 1, mu = M, Sigma = P))
        }
      }
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          obj$smooth[[j]] <- smooth.IWLS(obj$smooth[[j]])
        }
      }
    }
    obj
  }

  x <- sIWLS(x, ...)

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
  x$p.save <- if(is.null(x$p.save)) c("g", "tau2") else x$p.save
  if(is.null(x$samples)) {
    x$samples <- list()
    x$samples$g <- runif(ncol(x$X), 0.001, 0.002)
    x$samples$tau2 <- runif(1, 0.99, 1)
    x$samples$M <- rep(1, ncol(x$X))
    x$samples$P <- solve(crossprod(x$X))
    x$s.colnames <- if(is.null(x$s.colnames)) {
      c(paste("g[", 1:length(x$samples$g), "]", sep = ""), "tau2[1]")
    } else x$s.colnames
    x$samples$log.prior <- drop(-1 / x$samples$tau2 * crossprod(x$samples$g, x$S[[1]]) %*% x$samples$g)
    x$np <- ncol(x$X) + 1
  }

  if(is.null(x$update)) {
    if(is.null(x$rand)) {
      ## Is a list, if length sample > 1 means additional accept check in MCMC step.
      x$update <- function(x, z, eta, weights) {
        XW <- crossprod(x$X, diag(weights))
        P <- chol2inv(chol(XW %*% x$X + 1 / x$samples$tau2 * x$S[[1]]))
        M <- P %*% (XW %*% (z - eta))
        x$samples$g <- drop(mvrnorm(n = 1, mu = M, Sigma = P))
        a <- x$rank / 2 + x$a
        b <- 0.5 * crossprod(x$samples$g, x$S[[1]]) %*% x$samples$g + x$b
        x$samples$tau2 <- 1 / rgamma(1, a, b)
        x$samples$fs <- drop(x$X %*% x$samples$g)
        x$samples$log.prior <- drop(-1 / x$samples$tau2 * crossprod(x$samples$g, x$S[[1]]) %*% x$samples$g)
        x$samples$M <- M
        x$samples$P <- P
        x$samples
      }
    }
  }

  x
}


## Sampler based on IWLS proposals.
samplerIWLS <- function(x, n.iter = 1000, thin = 2, burnin = 200,
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
      if(!is.null(dim(obj$X))) {
        obj$s.alpha <- rep(0, n.save)
        obj$s.samples <- matrix(0, nrow = n.save, ncol = ncol(obj$X))
        obj$fp <- rep(0, nrow(obj$X))
      }
      if(length(obj$smooth)) {
        for(j in seq_along(obj$smooth)) {
          obj$smooth[[j]]$s.alpha <- rep(0, nrow = n.save)
          obj$smooth[[j]]$s.samples <- matrix(0, nrow = n.save, ncol = obj$smooth[[j]]$np)
          obj$smooth[[j]]$samples$fs <- rep(0, nrow(obj$smooth[[j]]$X))
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

  ## IWLS family specs.
  iwls <- family$iwls

  ## Start sampling
  for(i in 1:n.iter) {
    if(save <- i %in% iterthin)
      js <- which(iterthin == i)
      
    for(j in 1:np) {
      ## Compute weights.
      weights <- iwls[[nx[j]]]$weights(response, eta)

      ## Score.
      score <- iwls[[nx[j]]]$score(response, eta)

      ## Compute working observations.
      z <- eta[[nx[j]]] + 1 / weights * score

      ## Sample smooth effects
      if(length(x[[nx[j]]]$smooth)) {
        for(sj in seq_along(x[[nx[j]]]$smooth)) {
          ## Save old predictor.
          eta0 <- eta[[nx[j]]]

          ## Compute old log likelihood + log prior.
          p0 <- iwls$loglik(response, eta)

          ## Compute partial predictor.
          eta[[nx[j]]] <- eta0 - x[[nx[j]]]$smooth[[sj]]$samples$fs

          ## Sample.
          sms <- x[[nx[j]]]$smooth[[sj]]$update(x[[nx[j]]]$smooth[[sj]], z, eta[[nx[j]]], weights)

          ## Set up new predictor.
          eta[[nx[j]]] <- eta[[nx[j]]] + sms$fs

          ## Compute new log likelihood + log prior.
          p1 <- iwls$loglik(response, eta)

          ## Compute acceptance probablity.

          alpha <- (p1 + x[[nx[j]]]$smooth[[sj]]$samples$log.priors[1] + sms$log.priors[2]) -
            (p0 + sms$log.priors[1] + x[[nx[j]]]$smooth[[sj]]$samples$log.priors[2])

          ## If accepted, set current state to proposed state.
          accepted <- if(is.na(alpha)) FALSE else log(runif(1)) <= alpha
print(accepted)
          if(accepted) {
            x[[nx[j]]]$smooth[[sj]]$samples <- sms
          } else eta[[nx[j]]] <- eta0

          ## Save the samples and acceptance.
          if(save) {
            x[[nx[j]]]$smooth[[sj]]$s.alpha[js] <- accepted
            x[[nx[j]]]$smooth[[sj]]$s.samples[js, ] <- unlist(x[[nx[j]]]$smooth[[sj]]$samples[x[[nx[j]]]$smooth[[sj]]$p.save])
          }
        }
      }

      ## Sample parametric effects.
      if(!is.null(x[[nx[j]]]$X)) {
        ## Save old predictor.
        eta0 <- eta[[nx[j]]]

        ## Compute old log likelihood.
        p0 <- iwls$loglik(response, eta)

        ## Compute partial predictor.
        eta[[nx[j]]] <- eta0 - x[[nx[j]]]$fp

        ## Sample.
        ps <- x[[nx[j]]]$update(x[[nx[j]]]$X, z, eta[[nx[j]]], weights)

        ## Set up new predictor.
        eta[[nx[j]]] <- eta[[nx[j]]] + x[[nx[j]]]$X %*% ps

        ## Compute new log likelihood.
        p1 <- iwls$loglik(response, eta)

        ## Compute acceptance probablity.
        alpha <- p0 - p1

        ## If accepted, set current state to proposed state.
        accepted <- if(is.na(alpha)) FALSE else log(runif(1)) <= alpha
        if(accepted) {
          x[[nx[j]]]$samples <- ps
        } else eta[[nx[j]]] <- eta0

        ## Save the samples and acceptance.
        if(save) {
          x[[nx[j]]]$s.alpha[js] <- accepted
          x[[nx[j]]]$s.samples[js, ] <- ps
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

